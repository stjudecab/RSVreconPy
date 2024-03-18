#!/usr/bin/env python3

import os
import json
import re
import sys
import argparse
import subprocess
import glob
import fileinput
from yaml import safe_load


###################################################
##      load packages
###################################################

# Assuming RSV_functions is in the same directory as this script
from RSV_functions import get_sub_folders, elements_not_in_array, pct_sum, determine_subtype, processIGV
from Report_functions import generate_pdf_report

###################################################
##      functions
###################################################

def generate_csv_fasta(report, sequence_file, reference_folder_name, working_folder_name, ref_use):
    # get all reference DB
    all_ref_db = get_sub_folders(reference_folder_name)

    if ref_use != "all":
        ref_use_list = ref_use.split(',')
        not_exist = elements_not_in_array(ref_use_list, all_ref_db)
        if not_exist:
            string = ', '.join(not_exist)
            raise ValueError(f"These references can not be found under your reference folder: {string}")
        else:
            all_ref_db = ref_use_list

    # create output CSV
    CSV_header = 'Sample name,before_filtering_total_reads, before_filtering_q20_rate, before_filtering_q30_rate, after_filtering_total_reads, after_filtering_q20_rate, after_filtering_q30_rate, QC rate,'
    for ref_db in all_ref_db:
        CSV_header += f"Mapping {ref_db} Uniquely mapped reads %,Mapping {ref_db} MULTI-MAPPING READS %,Mapping {ref_db} UNMAPPED READS%,Mapping {ref_db} CHIMERIC READS%,"
    CSV_header += "Subtype suggestion"

    with open(report, 'w') as out:
        out.write(CSV_header + "\n")

    with open(sequence_file, 'w') as fasta:

        # collect information from each sample
        Sample_folders = [f for f in os.listdir(working_folder_name) if os.path.isdir(os.path.join(working_folder_name, f))]
        for cur_folder in sorted(Sample_folders):
            if os.path.isfile(cur_folder):
                continue
            if cur_folder == "Temp":
                continue

            print(f"Start process {cur_folder} ...")

            # Step 1: collect QC info from JSON file
            QC_json_file = os.path.join(working_folder_name, cur_folder, 'fastp.json')

            with open(QC_json_file, 'r') as txt:
                json_str = txt.read()

            decoded_json = json.loads(json_str)

            before = decoded_json['summary']['before_filtering']
            after = decoded_json['summary']['after_filtering']

            # Step 2: collect mapping info for all ref DB
            result_array = []

            for ref_db in all_ref_db:
                cur_map_log_file = os.path.join(working_folder_name, cur_folder, ref_db, 'Log.final.out')

                cur_map_hash = {}
                with open(cur_map_log_file, 'r') as fh:
                    for line in fh:
                        if '|' in line:
                            tmp_line = line
                            tmp_line = re.sub(r'[\s%:]', '', tmp_line)
                            tmp_arr = tmp_line.split('|')
                            cur_map_hash[tmp_arr[0]] = tmp_arr[1]

                result_array.append(cur_map_hash)

            # Step 3: output info to a CSV
            rate = after['total_reads'] / before['total_reads'] * 100
            QC_str = f"{before['total_reads']},{before['q20_rate']},{before['q30_rate']},{after['total_reads']},{after['q20_rate']},{after['q30_rate']},{rate},"
            map_str = ''
            uniquely_mapped_reads_hash = {}

            for i, result in enumerate(result_array):
                map_str += f"{result['Uniquelymappedreads']},{pct_sum(result['ofreadsmappedtomultipleloci'], result['ofreadsmappedtotoomanyloci'])},{pct_sum(result['ofreadsunmappedtoomanymismatches'], result['ofreadsunmappedtooshort'], result['ofreadsunmappedother'])},{result['ofchimericreads']},"
                uniquely_mapped_reads_hash[all_ref_db[i]] = float(result['Uniquelymappedreads'])

            subtype_str = determine_subtype(uniquely_mapped_reads_hash)

            sample = cur_folder
            sample = sample.split("/")[-1]

            with open(report, 'a') as out:
                out.write(f"{sample},{QC_str}{map_str}{subtype_str}\n")

            # Step 4: assemble genome sequence
            if subtype_str in ['SubTypeA','SubTypeB']:
                wig_file = os.path.join(working_folder_name, cur_folder, subtype_str, 'alignments.cov.wig')
                genome_sequence = processIGV(wig_file)

                fasta.write('>' + sample + "\n" + genome_sequence + "\n")

# generate tree 
def generate_phylogenetic_tree(root_file_path, working_folder, addtional_results = None):
    # set up file path
    reference_file = os.path.join(root_file_path, 'Resource','Know_ref.fasta')
    sequence_file = os.path.join(working_folder_name, "Sequence.fasta")
    merge_sequence_file = os.path.join(working_folder, "Sequence_for_tree.fasta")
    merge_alignment_file = os.path.join(working_folder, "Alignment_for_tree.fasta")
    tree_file = os.path.join(working_folder, "RSV.tree")

    # add additional dataset
    if addtional_results is not None:
        files = glob.glob(os.path.join(addtional_results, '*.fa'))
        print(files)
        # open a new file for writing
        with open(sequence_file, 'a') as seq_file:
            for filename in files:
                # open the file for reading
                with open(filename, 'r') as file:
                    # read all lines from the file
                    lines = file.readlines()
                    for line in lines:
                        if line.startswith('>'):
                            line = line.rstrip() + '|NEXT_RSV'
                            seq_file.write(line + "\n")
                        else:
                            seq_file.write(line)
                        

    if not os.path.isfile(tree_file):
        # merge sequencing results with reference sequences
        cmd = f"cat {sequence_file} {reference_file} > {merge_sequence_file}"
        subprocess.run(cmd, shell=True, cwd=working_folder)

        # make multiple sequence alignment
        cmd = f"muscle -in {merge_sequence_file} -out {merge_alignment_file}"
        subprocess.run(cmd, shell=True, cwd=working_folder)

        # make phylogenetic tree
        cmd = f"FastTree -nt -gtr {merge_alignment_file} > {tree_file}"
        subprocess.run(cmd, shell=True, cwd=working_folder)

###################################################
##      main
###################################################

if __name__ == "__main__":
    ###################################################
    ##      load parameters from YAML file
    ###################################################

    # Create command-line argument parser
    parser = argparse.ArgumentParser(description='Run pipeline to map NGS reads against reference genomes.sortedByCoord')
    parser.add_argument('yaml_file', help='Path to the YAML file containing the setting.')

    # Parse command-line arguments
    args = parser.parse_args()

    try:
        with open(args.yaml_file, 'r') as f :
            config = safe_load(f)
    except:
        print('Can not open this file!')
        sys.exit()

    data_folder_name = config['DATA_DIR']
    reference_folder_name = config['REFERENCE_DIR']
    working_folder_name = config['OUTPUT_DIR']

    star_ThreadN = config.get('STAR_THREAD_N', 16)
    ref_use = config.get('REF_USE', 'all')

    ###################################################
    ##      setup working dir
    ###################################################
    working_folder_name = os.path.normpath(working_folder_name)
    if not os.path.exists(working_folder_name):
        os.mkdir(working_folder_name)

    root_file_path = os.path.dirname(os.path.realpath(__file__))

    ###################################################
    ##      start process
    ###################################################

    report = os.path.join(working_folder_name, "Report.csv")
    sequence_file = os.path.join(working_folder_name, "Sequence.fasta")
    print("Aggregating results ...")
    generate_csv_fasta(report, sequence_file, reference_folder_name, working_folder_name, ref_use)

    ###################################################
    ##      generate phylogenetic tree with reference
    ###################################################

    print("Generating tree ...")
    addtional_results_path = os.path.join(working_folder_name, '..', 'consensus')
    generate_phylogenetic_tree(root_file_path, working_folder_name, addtional_results_path)

    ###################################################
    ##      generate pdf report
    ###################################################

    print("Generating Pdf report...")
    generate_pdf_report(report, working_folder_name)

    ###################################################
    ##      program finished
    ###################################################

    print("Report and sequences have been generated!")
