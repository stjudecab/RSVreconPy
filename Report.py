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
import pandas as pd
import shutil


###################################################
##      load packages
###################################################

# Assuming RSV_functions is in the same directory as this script
from RSV_functions import get_sub_folders, elements_not_in_array, pct_sum, determine_subtype, processIGV, get_genotype_res, extract_F_protein, detect_F_mutation, extract_gene_covarage
from Report_functions import generate_pdf_report

###################################################
##      functions
###################################################

def generate_csv_fasta(report, sequence_file, reference_folder_name, working_folder_name, sequence_file_a, sequence_file_b, igv_cutoff):
    # create output CSV
    CSV_header = 'Sample name,before_filtering_total_reads, before_filtering_q20_rate, before_filtering_q30_rate, after_filtering_total_reads, after_filtering_q20_rate, after_filtering_q30_rate, QC rate,'
    CSV_header += "Uniquely mapped reads %,MULTI-MAPPING READS %,UNMAPPED READS%,CHIMERIC READS%,"
    CSV_header += "Subtype suggestion,reference_accession,ref_subtype,"
    CSV_header += "F protein mutations,"
    CSV_header += "NS1_cov,NS2_cov,N_cov,P_cov,M_cov,SH_cov,G_cov,F_cov,M2-1_cov,M2-2_cov,L_cov"

    subtype_a_names = []
    subtype_b_names = []

    with open(report, 'w') as out:
        out.write(CSV_header + "\n")

    with open(report, 'a') as out:

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
            if os.path.exists(cur_map_log_file):
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
                genotype_file = os.path.join(working_folder_name, cur_folder, 'Genotype', 'Genotype.txt')
                if os.path.isfile(genotype_file):
                    subtype_str, blast_pct_identity, blast_alignment_length, starin = get_genotype_res(genotype_file)
                    ##### fix later
                    if int(blast_alignment_length) < 1500:
                        subtype_str = 'Not RSV'
                else:
                    subtype_str = 'Not RSV'

                sample = cur_folder.split("/")[-1]

                kma_out_file = os.path.join(working_folder_name, cur_folder, 'KMA', cur_folder + '.res')
                ref_str = determine_subtype(kma_out_file, reference_folder_name, 10)
                ref_str = ref_str.split(",")

                out.write(f"{sample},{QC_str}{map_str}{subtype_str},{ref_str[1]},{ref_str[2]}")

                # Step 4: get assembled genome sequence
                genome_sequence_file = os.path.join(working_folder_name, cur_folder, 'sequence.fasta')
                genome_sequence = ''
                with open(genome_sequence_file, 'r') as read_seq:
                    for line in read_seq:
                        if line.startswith('>'):
                            pass
                        else:
                            genome_sequence += line.strip()   

                with open(sequence_file, 'w') as fasta, open(sequence_file_a, 'w') as fasta_a, open(sequence_file_b, 'w') as fasta_b:
                    if subtype_str == 'Not RSV':
                        sequence_file_error = sequence_file + '.err'
                        with open(sequence_file_error, 'a') as fasta_err:
                            fasta_err.write('>' + sample + "\n" + genome_sequence + "\n")
                    elif subtype_str[0] == 'A':
                        fasta.write('>' + sample + "\n" + genome_sequence + "\n")
                        fasta_a.write('>' + sample + "\n" + genome_sequence + "\n")
                        subtype_a_names.append(sample)
                    else:
                        fasta.write('>' + sample + "\n" + genome_sequence + "\n")
                        fasta_b.write('>' + sample + "\n" + genome_sequence + "\n")
                        subtype_b_names.append(sample)

                # step 5: identify F protein mutations
                reference_sequence_file = os.path.join(working_folder_name, cur_folder, 'reference', f"{ref_str[1]}.fasta")
                gff_file = os.path.join(working_folder_name, cur_folder, 'reference', f"{ref_str[1]}.gff")
                if subtype_str == 'Not RSV':
                    mutations_str = ''
                elif subtype_str[0] == 'A':
                    mutation_file = os.path.join(reference_folder_name, 'Genotype_ref','RSV_A_F_Mutation.csv')
                    assembled_f_protein_sequence = extract_F_protein(sequence_file, reference_sequence_file, gff_file)
                    #print(assembled_f_protein_sequence)
                    F_mutation = detect_F_mutation(assembled_f_protein_sequence, mutation_file)
                    mutations_str = [f"{ele[0]}({ele[1]})" for ele in F_mutation]
                    mutations_str = "|".join(mutations_str)
                else:
                    mutation_file = os.path.join(reference_folder_name, 'Genotype_ref','RSV_B_F_Mutation.csv')
                    assembled_f_protein_sequence = extract_F_protein(sequence_file, reference_sequence_file, gff_file)
                    #print(assembled_f_protein_sequence)
                    F_mutation = detect_F_mutation(assembled_f_protein_sequence, mutation_file)
                    mutations_str = [f"{ele[0]}({ele[1]})" for ele in F_mutation]
                    mutations_str = "|".join(mutations_str)
                
                out.write(","+mutations_str)

                # step 6: calculate coverage for each gene segment
                wig_file = os.path.join(working_folder_name, cur_folder, 'mapping', 'alignments.cov.wig')
                coverage_by_gene = extract_gene_covarage(wig_file, gff_file)
                rsvgenes = ['CDS_1','CDS_2','CDS_3','CDS_4','CDS_5','CDS_6','CDS_7','CDS_8','CDS_9','CDS_10','CDS_11']
                coverages = [coverage_by_gene[gene] for gene in rsvgenes]
                coverages_str = ','.join(str(x) for x in coverages)

                out.write(","+coverages_str + "\n")

            else:
                sample = cur_folder.split("/")[-1]
                QC_str = f"{before['total_reads']},{before['q20_rate']},{before['q30_rate']},{after['total_reads']},{after['q20_rate']},{after['q30_rate']},{rate},"
                map_str = f"0,0,100,0,"
                subtype_str = 'Not RSV'
                out.write(f"{sample},{QC_str}{map_str}{subtype_str},NA,Not RSV,,0,0,0,0,0,0,0,0,0,0,0\n")

    return subtype_a_names, subtype_b_names

# generate tree 
def generate_phylogenetic_tree(root_file_path, reference_folder_name, working_folder_name, subtype_a_names, subtype_b_names, Temp_folder_name, addtional_results = None):
    # set up file path
    render_file = os.path.join(root_file_path, 'RenderTree.R')
    out_group_file = os.path.join(reference_folder_name, 'Genotype_ref','out_group_setting.txt')
    out_group_df = pd.read_csv(out_group_file, sep=',', index_col=0)

    # generate subtype A tree XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    if len(subtype_a_names) > 0:
        reference_file = os.path.join(reference_folder_name, 'Genotype_ref','representative_ref_A.fasta')
        reference_anno_file = os.path.join(reference_folder_name, 'Genotype_ref','representative_ref_A.csv')
        out_group_strain = out_group_df.loc['A','strain']
        color_file = os.path.join(reference_folder_name, 'Genotype_ref','color_A.csv')

        sequence_file = os.path.join(working_folder_name, "Sequence_A.fasta")
        merge_sequence_file = os.path.join(Temp_folder_name, "Sequence_for_tree_A.fasta")
        merge_alignment_file = os.path.join(Temp_folder_name, "Alignment_for_tree_A.fasta")
        tree_file = os.path.join(Temp_folder_name, "RSV_A.nwk")
        tree_annotation_file = os.path.join(Temp_folder_name, "RSV_A.csv")
        tree_pic_file = os.path.join(Temp_folder_name, "RSV_A.png")
        tree_log = os.path.join(Temp_folder_name, "RSV_A_log.txt")

        shutil.copyfile(reference_anno_file, tree_annotation_file)

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
                                    text = f"{seq_name},In-house,Query\n{seq_name}|NEXT_RSV,NEXT_RSV,Query\n"
                                    anno_file.write(text)
                                else:
                                    break
                            else:
                                seq_file.write(line)
        else:
            with open(tree_annotation_file, 'a') as anno_file:
                for seq_name in subtype_a_names:
                    text = f"{seq_name},Query,Query\n"
                    anno_file.write(text)

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
            cmd = f"Rscript {render_file} {tree_file} '{out_group_strain}' {tree_annotation_file} {tree_pic_file} {color_file} &>{tree_log}"
            subprocess.run(cmd, shell=True, cwd=Temp_folder_name)


    # generate subtype B tree XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    if len(subtype_b_names) > 0:
        reference_file = os.path.join(reference_folder_name, 'Genotype_ref','representative_ref_B.fasta')
        reference_anno_file = os.path.join(reference_folder_name, 'Genotype_ref','representative_ref_B.csv')
        out_group_strain = out_group_df.loc['B','strain']
        color_file = os.path.join(reference_folder_name, 'Genotype_ref','color_B.csv')

        sequence_file = os.path.join(working_folder_name, "Sequence_B.fasta")
        merge_sequence_file = os.path.join(Temp_folder_name, "Sequence_for_tree_B.fasta")
        merge_alignment_file = os.path.join(Temp_folder_name, "Alignment_for_tree_B.fasta")
        tree_file = os.path.join(Temp_folder_name, "RSV_B.nwk")
        tree_annotation_file = os.path.join(Temp_folder_name, "RSV_B.csv")
        tree_pic_file = os.path.join(Temp_folder_name, "RSV_B.png")
        tree_log = os.path.join(Temp_folder_name, "RSV_B_log.txt")

        shutil.copyfile(reference_anno_file, tree_annotation_file)

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
                                    text = f"{seq_name},In-house,Query\n{seq_name}|NEXT_RSV,NEXT_RSV,Query\n"
                                    anno_file.write(text)
                                else:
                                    break
                            else:
                                seq_file.write(line)
        else:
            with open(tree_annotation_file, 'a') as anno_file:
                for seq_name in subtype_b_names:
                    text = f"{seq_name},Query,Query\n"
                    anno_file.write(text)

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
            cmd = f"Rscript {render_file} {tree_file} '{out_group_strain}' {tree_annotation_file} {tree_pic_file} {color_file} &>{tree_log}"
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
    igv_cutoff = config.get('IGV_CUTOFF', 50)

    ###################################################
    ##      setup working dir
    ###################################################
    working_folder_name = os.path.normpath(working_folder_name)

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
    mapres_folder = os.path.join(working_folder_name, "Mapping")
    print("Aggregating results ...")
    subtype_a_names, subtype_b_names = generate_csv_fasta(report, sequence_file, reference_folder_name, mapres_folder, sequence_file_a, sequence_file_b, igv_cutoff)

    ###################################################
    ##      generate phylogenetic tree with reference
    ###################################################

    print("Generating tree ...")
    addtional_results_path = os.path.join(working_folder_name, '..', 'consensus')
    if os.path.exists(addtional_results_path):
        pass
    else:
        addtional_results_path = None
    generate_phylogenetic_tree(root_file_path, reference_folder_name, working_folder_name, subtype_a_names, subtype_b_names, Temp_folder_name, addtional_results_path)

    ###################################################
    ##      generate pdf report
    ###################################################

    print("Generating Pdf report...")
    generate_pdf_report(report, working_folder_name, mapres_folder, igv_cutoff)

    ###################################################
    ##      program finished
    ###################################################

    print("Report and sequences have been generated!")
