#!/usr/bin/env python3

import os
import sys
import subprocess
import pandas as pd
import shutil

def genotype_call_whole_genome(query_file_path, ref_db_path, meta_file_path, output_path, reference_folder_name, render_file):
    
    ##############################################################################
    ### step 1    blast against reference database
    ##############################################################################
    blast_out_file = os.path.join(output_path, 'blastn_res.tsv')
    cmd = f"blastn -query {query_file_path} -db {ref_db_path} -out {blast_out_file} -outfmt 6 -num_threads 1 -evalue 1e-5 -perc_identity 90"
    print(cmd)
    subprocess.run(cmd, shell=True, cwd=output_path)

    ##############################################################################
    ### step 2    fetch results 
    ##############################################################################

    blast_df = pd.read_csv(blast_out_file, sep='\t', header=None)
    blast_df.columns = ['Query', 'Subject', 'Pct_identity','Alignment_length','Number_of_mismatches','Number_of_gap','Start_in_query','End_in_query','Start_in_subject','End_in_subject','E_value','Bit_score']
    blast_df = blast_df.sort_values(by=['Bit_score'], ascending = False)

    query_name = blast_df.loc[0, 'Query']

    meta_df = pd.read_csv(meta_file_path, sep='\t', index_col=0)

    clade = ''
    strain = ''
    blast_pct_identity = 0
    blast_alignment_length = 0
    sel_row_index = 0
    while clade == '':
        Selected_ref_name = blast_df.loc[sel_row_index, 'Subject']
        blast_pct_identity = blast_df.loc[sel_row_index, 'Pct_identity']
        blast_alignment_length = blast_df.loc[sel_row_index, 'Alignment_length']

        clade = meta_df.loc[Selected_ref_name, 'clade']
        strain = meta_df.loc[Selected_ref_name, 'strain']
        sel_row_index += 1
    
    genotype_file = os.path.join(output_path, 'Genotype.txt')
    with open(genotype_file, 'w') as geno_file:
        geno_file.write(f"{clade},{blast_pct_identity},{blast_alignment_length},{strain}")


    ##############################################################################
    ### step 3    generate tree
    ##############################################################################

    out_group_file = os.path.join(reference_folder_name, 'Genotype_ref','out_group_setting.txt')
    out_group_df = pd.read_csv(out_group_file, sep=',', index_col=0)

    if clade[0] == 'A':
        reference_file = os.path.join(reference_folder_name, 'Genotype_ref','representative_ref_A.fasta')
        reference_anno_file = os.path.join(reference_folder_name, 'Genotype_ref','representative_ref_A.csv')
        out_group_strain = out_group_df.loc['A','strain']
        color_file = os.path.join(reference_folder_name, 'Genotype_ref','color_A.csv')
    elif clade[0] == 'B':
        reference_file = os.path.join(reference_folder_name, 'Genotype_ref','representative_ref_B.fasta')
        reference_anno_file = os.path.join(reference_folder_name, 'Genotype_ref','representative_ref_B.csv')
        out_group_strain = out_group_df.loc['B','strain']
        color_file = os.path.join(reference_folder_name, 'Genotype_ref','color_B.csv')
    else:
        return 0

    # set file locations
    merge_sequence_file = os.path.join(output_path, "Sequence_for_tree.fasta")
    merge_alignment_file = os.path.join(output_path, "Alignment_for_tree.fasta")
    tree_file = os.path.join(output_path, "genotype.nwk")
    tree_annotation_file = os.path.join(output_path, "genotype.csv")

    tree_pic_file = os.path.join(output_path, "genotype.png")
    tree_log = os.path.join(output_path, "genotype_log.txt")

    shutil.copyfile(reference_anno_file, tree_annotation_file)

    # open a new file for writing
    with open(tree_annotation_file, 'a') as anno_file:
        text = f"{query_name},Query,Query\n"
        anno_file.write(text)

    # generate tree
    if not os.path.isfile(tree_file):
        # merge sequencing results with reference sequences
        cmd = f"cat {query_file_path} {reference_file} > {merge_sequence_file}"
        subprocess.run(cmd, shell=True, cwd=output_path)

        # make multiple sequence alignment
        cmd = f"muscle -in {merge_sequence_file} -out {merge_alignment_file} &> genotype_muscle_log.txt"
        subprocess.run(cmd, shell=True, cwd=output_path)

        # make phylogenetic tree
        cmd = f"FastTree -nt -gtr {merge_alignment_file} > {tree_file} 2> genotype_FastTree_log.txt"
        subprocess.run(cmd, shell=True, cwd=output_path)

    # render tree
    if os.path.isfile(tree_file):
        cmd = f"Rscript {render_file} {tree_file} '{out_group_strain}' {tree_annotation_file} {tree_pic_file} {color_file} &>{tree_log}"
        #print(cmd)
        subprocess.run(cmd, shell=True, cwd=output_path)

if __name__ == "__main__":
    query_file_path = sys.argv[1]
    ref_db_path = sys.argv[2]
    meta_file_path = sys.argv[3]
    output_path = sys.argv[4]
    reference_folder_name = sys.argv[5]
    render_file = sys.argv[6]

    #print(sys.argv)

    genotype_call_whole_genome(query_file_path, ref_db_path, meta_file_path, output_path, reference_folder_name, render_file)