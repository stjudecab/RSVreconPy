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

def generate_csv_fasta(report, sequence_file, reference_folder_name, working_folder_name, sequence_file_a, sequence_file_b):
    # create output CSV
    CSV_header = 'Sample name,before_filtering_total_reads, before_filtering_q20_rate, before_filtering_q30_rate, after_filtering_total_reads, after_filtering_q20_rate, after_filtering_q30_rate, QC rate,'
    CSV_header += "Uniquely mapped reads %,MULTI-MAPPING READS %,UNMAPPED READS%,CHIMERIC READS%,"
    CSV_header += "Subtype suggestion, reference_accession, ref_subtype"

    subtype_a_names = []
    subtype_b_names = []

    with open(report, 'w') as out:
        out.write(CSV_header + "\n")

    with open(sequence_file, 'w') as fasta, open(sequence_file_a, 'w') as fasta_a, open(sequence_file_b, 'w') as fasta_b, open(report, 'a') as out:

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

            # Step 2: collect mapping info
            cur_map_log_file = os.path.join(working_folder_name, cur_folder, 'mapping', 'Log.final.out')

            cur_map_hash = {}
            with open(cur_map_log_file, 'r') as fh:
                for line in fh:
                    if '|' in line:
                        tmp_line = line
                        tmp_line = re.sub(r'[\s%:]', '', tmp_line)
                        tmp_arr = tmp_line.split('|')
                        cur_map_hash[tmp_arr[0]] = tmp_arr[1]

            # Step 3: output info to a CSV
            rate = after['total_reads'] / before['total_reads'] * 100
            QC_str = f"{before['total_reads']},{before['q20_rate']},{before['q30_rate']},{after['total_reads']},{after['q20_rate']},{after['q30_rate']},{rate},"

            map_str = f"{cur_map_hash['Uniquelymappedreads']},{pct_sum(cur_map_hash['ofreadsmappedtomultipleloci'], cur_map_hash['ofreadsmappedtotoomanyloci'])},{pct_sum(cur_map_hash['ofreadsunmappedtoomanymismatches'], cur_map_hash['ofreadsunmappedtooshort'], cur_map_hash['ofreadsunmappedother'])},{cur_map_hash['ofchimericreads']},"

            # determine subtype
            kma_out_file = os.path.join(working_folder_name, cur_folder, 'KMA', cur_folder + '.res')
            subtype_str = determine_subtype(kma_out_file, reference_folder_name)

            sample = cur_folder
            sample = sample.split("/")[-1]

            out.write(f"{sample},{QC_str}{map_str}{subtype_str}\n")

            # Step 4: assemble genome sequence
            wig_file = os.path.join(working_folder_name, cur_folder, 'mapping', 'alignments.cov.wig')
            genome_sequence = processIGV(wig_file)
            my_subtype = subtype_str.split(',')[0]
            if my_subtype in ['SubtypeA','SubtypeB']:
                fasta.write('>' + sample + "\n" + genome_sequence + "\n")
                if my_subtype == 'SubtypeA':
                    fasta_a.write('>' + sample + "\n" + genome_sequence + "\n")
                    subtype_a_names.append(sample)
                else:
                    fasta_b.write('>' + sample + "\n" + genome_sequence + "\n")
                    subtype_b_names.append(sample)
            else:
                sequence_file_error = sequence_file + '.err'
                with open(sequence_file_error, 'a') as fasta_err:
                    fasta_err.write('>' + sample + "\n" + genome_sequence + "\n")

    return subtype_a_names, subtype_b_names

# generate tree 
def generate_phylogenetic_tree(root_file_path, working_folder, subtype_a_names, subtype_b_names, Temp_folder_name, addtional_results = None):
    # set up file path
    reference_file_a = os.path.join(root_file_path, 'Resource','Know_ref_A.fasta')
    reference_file_b = os.path.join(root_file_path, 'Resource','Know_ref_B.fasta')
    sequence_file_a = os.path.join(working_folder_name, "Sequence_A.fasta")
    sequence_file_b = os.path.join(working_folder_name, "Sequence_B.fasta")
    render_file = os.path.join(root_file_path, 'RenderTree.R')


    # generate subtype A tree
    if len(subtype_a_names) > 0:
        reference_file = os.path.join(root_file_path, 'Resource','Know_ref_A.fasta')
        sequence_file = os.path.join(working_folder_name, "Sequence_A.fasta")
        merge_sequence_file = os.path.join(Temp_folder_name, "Sequence_for_tree_A.fasta")
        merge_alignment_file = os.path.join(Temp_folder_name, "Alignment_for_tree_A.fasta")
        tree_file = os.path.join(Temp_folder_name, "RSV_A.nwk")
        tree_annotation_file = os.path.join(Temp_folder_name, "RSV_A.csv")
        tree_pic_file = os.path.join(Temp_folder_name, "RSV_A.png")
        tree_log = os.path.join(Temp_folder_name, "RSV_A_log.txt")

        with open(tree_annotation_file, 'w') as anno_file:
            text = "VR_1540|ATCC,Reference\nVR_26|ATCC,Reference\nNR_28527|BEI,Reference\nNR_28528|BEI,Reference\nNR_28529|BEI,Reference\nNR_48671|BEI,Reference\nRSVA/Human/USA/78G-104-01/1978|KU316155,Reference\n"
            anno_file.write(text)

        # add additional dataset
        if addtional_results is not None:
            files = glob.glob(os.path.join(addtional_results, '*.fa'))
            #print(files)
            # open a new file for writing
            with open(sequence_file, 'a') as seq_file, open(tree_annotation_file, 'a') as anno_file:
                for filename in files:
                    # open the file for reading
                    with open(filename, 'r') as file:
                        # read all lines from the file
                        lines = file.readlines()
                        for line in lines:
                            if line.startswith('>'):
                                seq_name = re.sub('>','',line)
                                seq_name = re.sub('_R_.+','',seq_name)
                                seq_name = seq_name.rstrip()
                                #print(seq_name)
                                if seq_name in subtype_a_names:
                                    seq_file.write(f'>{seq_name}|NEXT_RSV\n')

                                    text = f"{seq_name},In-house_Results\n{seq_name}|NEXT_RSV,NEXT_RSV_Results\n"
                                    anno_file.write(text)
                                else:
                                    break
                            else:
                                seq_file.write(line)

        # generate tree
        if not os.path.isfile(tree_file):
            # merge sequencing results with reference sequences
            cmd = f"cat {sequence_file} {reference_file} > {merge_sequence_file}"
            subprocess.run(cmd, shell=True, cwd=Temp_folder_name)

            # make multiple sequence alignment
            cmd = f"muscle -in {merge_sequence_file} -out {merge_alignment_file} &> subtypeA_muscle_log.txt"
            subprocess.run(cmd, shell=True, cwd=Temp_folder_name)

            # make phylogenetic tree
            cmd = f"FastTree -nt -gtr {merge_alignment_file} > {tree_file} 2> subtypeA_FastTree_log.txt"
            subprocess.run(cmd, shell=True, cwd=Temp_folder_name)

        # render tree
        if os.path.isfile(tree_file):
            #print("run tree A")
            cmd = f"Rscript {render_file} {tree_file} 'RSVA/Human/USA/78G-104-01/1978|KU316155' {tree_annotation_file} {tree_pic_file} &>{tree_log}"
            subprocess.run(cmd, shell=True, cwd=Temp_folder_name)


    # generate subtype B tree
    if len(subtype_b_names) > 0:
        reference_file = os.path.join(root_file_path, 'Resource','Know_ref_B.fasta')
        sequence_file = os.path.join(working_folder_name, "Sequence_B.fasta")
        merge_sequence_file = os.path.join(Temp_folder_name, "Sequence_for_tree_B.fasta")
        merge_alignment_file = os.path.join(Temp_folder_name, "Alignment_for_tree_B.fasta")
        tree_file = os.path.join(Temp_folder_name, "RSV_B.nwk")
        tree_annotation_file = os.path.join(Temp_folder_name, "RSV_B.csv")
        tree_pic_file = os.path.join(Temp_folder_name, "RSV_B.png")
        tree_log = os.path.join(Temp_folder_name, "RSV_B_log.txt")

        with open(tree_annotation_file, 'w') as anno_file:
            text = "VR_1400|ATCC,Reference\nVR_955|ATCC,Reference\nVR_1580|ATCC,Reference\nRSVB/Human/USA/81G-027-01/1977|KU316116,Reference\n"
            anno_file.write(text)

        # add additional dataset
        if addtional_results is not None:
            files = glob.glob(os.path.join(addtional_results, '*.fa'))
            #print(files)
            # open a new file for writing
            with open(sequence_file_b, 'a') as seq_file, open(tree_annotation_file, 'a') as anno_file:
                for filename in files:
                    # open the file for reading
                    with open(filename, 'r') as file:
                        # read all lines from the file
                        lines = file.readlines()
                        for line in lines:
                            if line.startswith('>'):
                                seq_name = re.sub('>','',line)
                                seq_name = re.sub('_R_.+','',seq_name)
                                seq_name = seq_name.rstrip()
                                #print(seq_name)
                                if seq_name in subtype_b_names:
                                    seq_file.write(f'>{seq_name}|NEXT_RSV\n')

                                    text = f"{seq_name},In-house_Results\n{seq_name}|NEXT_RSV,NEXT_RSV_Results\n"
                                    anno_file.write(text)
                                else:
                                    break
                            else:
                                seq_file.write(line)

        # generate tree
        if not os.path.isfile(tree_file):
            # merge sequencing results with reference sequences
            cmd = f"cat {sequence_file} {reference_file} > {merge_sequence_file}"
            subprocess.run(cmd, shell=True, cwd=Temp_folder_name)

            # make multiple sequence alignment
            cmd = f"muscle -in {merge_sequence_file} -out {merge_alignment_file} &> subtypeB_muscle_log.txt"
            subprocess.run(cmd, shell=True, cwd=Temp_folder_name)

            # make phylogenetic tree
            cmd = f"FastTree -nt -gtr {merge_alignment_file} > {tree_file} 2> subtypeB_FastTree_log.txt"
            subprocess.run(cmd, shell=True, cwd=Temp_folder_name)

        # render tree
        if os.path.isfile(tree_file):
            #print("run tree B")
            cmd = f"Rscript {render_file} {tree_file} 'RSVB/Human/USA/81G-027-01/1977|KU316116' {tree_annotation_file} {tree_pic_file} &>{tree_log}"
            subprocess.run(cmd, shell=True, cwd=Temp_folder_name)

    # job done!
    return True

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

    Temp_folder_name = os.path.join(working_folder_name, 'Temp')
    if not os.path.exists(Temp_folder_name):
        os.mkdir(Temp_folder_name)

    root_file_path = os.path.dirname(os.path.realpath(__file__))

    ###################################################
    ##      start process
    ###################################################

    report = os.path.join(working_folder_name, "Report.csv")
    sequence_file = os.path.join(working_folder_name, "Sequence.fasta")
    sequence_file_a = os.path.join(working_folder_name, "Sequence_A.fasta")
    sequence_file_b = os.path.join(working_folder_name, "Sequence_B.fasta")
    print("Aggregating results ...")
    subtype_a_names, subtype_b_names = generate_csv_fasta(report, sequence_file, reference_folder_name, working_folder_name, sequence_file_a, sequence_file_b)

    ###################################################
    ##      generate phylogenetic tree with reference
    ###################################################

    print("Generating tree ...")
    addtional_results_path = os.path.join(working_folder_name, '..', 'consensus')
    generate_phylogenetic_tree(root_file_path, working_folder_name, subtype_a_names, subtype_b_names, Temp_folder_name, addtional_results_path)

    ###################################################
    ##      generate pdf report
    ###################################################

    print("Generating Pdf report...")
    generate_pdf_report(report, working_folder_name)

    ###################################################
    ##      program finished
    ###################################################

    print("Report and sequences have been generated!")
