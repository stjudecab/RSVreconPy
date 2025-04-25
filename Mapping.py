#!/usr/bin/env python3

import os
import sys
import subprocess
import pandas as pd
import shutil
from Genotyping import genotype_call_whole_genome, genotype_call_G_protein
from RSV_functions import determine_subtype, processIGV, fetch_record_from_JSON

def sample_mapping(sample_id, working_folder_name, original_read1, original_read2, reference_folder_name, star_ThreadN, igv_cutoff, run_eachstep):
    # skip Undetermined reads
    if "Undetermined" in sample_id:
        return

    cur_folder = os.path.join(working_folder_name, sample_id)
    if not os.path.exists(cur_folder):
        os.mkdir(cur_folder)

    # make log folder
    log_cur_folder = os.path.join(cur_folder, 'log')
    if not os.path.exists(log_cur_folder):
        os.mkdir(log_cur_folder)

    ##############################################################################
    ### step 1    QC using fastp
    ##############################################################################
    print(f"\t######################\tStep 1 Quality check and trimming using fastp ... \t######################")

    read1_trim = "reads_trimed_R1.fastq"
    read2_trim = "reads_trimed_R2.fastq"
    log_file = os.path.join(log_cur_folder, 'fastp_log.txt')

    cmd = f"fastp -i \"{original_read1}\" -I \"{original_read2}\" -o {read1_trim} -O {read2_trim} &> {log_file}"
    print(cmd)
    if run_eachstep['fastp'] is True:
        subprocess.run(cmd, shell=True, cwd=cur_folder)

    ##############################################################################
    ### step 2    KMA and building reference
    ##############################################################################
    print(f"\t######################\tStep 2 Identifying the best reference ...         \t######################")

    kma_db = os.path.join(reference_folder_name, 'RSV_DB')
    kma_folder = os.path.join(cur_folder, 'KMA')
    if not os.path.exists(kma_folder):
        os.mkdir(kma_folder)

    kma_out = os.path.join(kma_folder, sample_id)
    log_file = os.path.join(log_cur_folder, 'kma_log.txt')

    kma_cmd = f"kma -ipe {read1_trim} {read2_trim} -o {kma_out} -t_db {kma_db} -1t1 -t 3 -mrs 0.9 -mq 50 &> {log_file}"
    print(kma_cmd)

    if run_eachstep['kma'] is True:
        subprocess.run(kma_cmd, shell=True, cwd=cur_folder)

    # read info from KMA results
    kma_out_file = kma_out + '.res'
    kma_df = pd.read_csv(kma_out_file, sep='\t')
    kma_df = kma_df.sort_values('Score', ascending=False)
    Selected_ref_name = kma_df.iloc[0, 0]
    #print(Selected_ref_name)

    # fetch GFF and FASTA for reference strain
    ref_gff_file = os.path.join(reference_folder_name, 'RSV_reference_GFF.json')
    ref_fasta_file = os.path.join(reference_folder_name, 'RSV_reference_FASTA.json')
    original_gff = fetch_record_from_JSON(Selected_ref_name, ref_gff_file)
    original_fasta = fetch_record_from_JSON(Selected_ref_name, ref_fasta_file)

    # create output folder
    ref_folder = os.path.join(cur_folder, 'reference')
    if not os.path.exists(ref_folder):
        os.mkdir(ref_folder)

    # write selected reference
    new_fasta_file = os.path.join(ref_folder, Selected_ref_name + '.fasta')
    new_gff_file = os.path.join(ref_folder, Selected_ref_name + '.gff')
    with open(new_fasta_file, 'w') as fasta:
        fasta.write(original_fasta)
    with open(new_gff_file, 'w') as fasta:
        fasta.write(original_gff)

    ref_output_dir = os.path.join(ref_folder, Selected_ref_name)
    log_file = os.path.join(log_cur_folder, 'bwa_build_genome_log.txt')
    star_command = f"bwa index -p {ref_output_dir} {new_fasta_file} &> {log_file}"
    print(star_command)

    if run_eachstep['star_index'] is True:
        subprocess.run(star_command, shell=True, cwd=cur_folder)

    ##############################################################################
    ## step 3    reads mapping
    ##############################################################################
    print(f"\t######################\tStep 3 mapping ...                                \t######################")

    map_folder = os.path.join(cur_folder, 'mapping')
    if not os.path.exists(map_folder):
        os.mkdir(map_folder)

    log_file = os.path.join(log_cur_folder, 'mapping_log.txt')
    sam_file = os.path.join(map_folder, 'sample_aligned.sam')
    cmd_mapping = f"bwa mem -t {star_ThreadN} {ref_output_dir} {read1_trim} {read2_trim} > {sam_file} 2> {log_file}"
    print(cmd_mapping)

    if run_eachstep['star'] is True:
        subprocess.run(cmd_mapping, shell=True, cwd=cur_folder)

    ##############################################################################
    ## step 4    index using Samtools
    ##############################################################################
    print(f"\t######################\tStep 4 Index and Count ...                        \t######################")

    log_file = os.path.join(log_cur_folder, 'samtools_log.txt')
    bam_file = os.path.join(map_folder, 'sample_aligned.bam')
    flagstat_file = os.path.join(map_folder, 'flagstat.txt')
    sort_bam_file = os.path.join(map_folder, 'sample_sorted.bam')
    cmd_2bam = f"samtools view -@ {star_ThreadN} -Sb {sam_file} > {bam_file}"
    print(cmd_2bam)
    cmd_sort = f"samtools sort -@ {star_ThreadN} -o {sort_bam_file} {bam_file}"
    print(cmd_sort)
    cmd_stat = f"samtools flagstat {sort_bam_file} > {flagstat_file}"
    print(cmd_stat)
    cmd_index = f"samtools index {sort_bam_file} &> {log_file}"
    print(cmd_index)

    if run_eachstep['samtools'] is True:
        subprocess.run(cmd_2bam, shell=True, cwd=map_folder)
        subprocess.run(cmd_sort, shell=True, cwd=map_folder)
        subprocess.run(cmd_stat, shell=True, cwd=map_folder)
        subprocess.run(cmd_index, shell=True, cwd=map_folder)

    ##############################################################################
    ## step 4.5    Count using IGVtools 
    ##############################################################################

    log_file = os.path.join(log_cur_folder, 'igvtools_log.txt')
    wig_file = os.path.join(map_folder, 'alignments.cov.wig')
    cmd_count = f"igvtools count -z 5 -w 1 --bases {sort_bam_file} {wig_file} {new_fasta_file} &> {log_file}"
    print(cmd_count)

    if run_eachstep['igvtools'] is True:
        subprocess.run(cmd_count, shell=True, cwd=map_folder)

    ##############################################################################
    ## step 5    Assemble sequence, Genotyping call using whole genome and G sequence
    ##############################################################################
    print(f"\t######################\tStep 5 Calling genotypes ...                      \t######################")

    # assemble genomic sequence
    genome_sequence = processIGV(wig_file, new_fasta_file, int(igv_cutoff))
    query_file_path = os.path.join(cur_folder, 'sequence.fasta')
    with open(query_file_path, 'w') as fasta:
        fasta.write('>' + sample_id + "\n" + genome_sequence + "\n")

    # create folder for genotyping
    genotype_folder = os.path.join(cur_folder, 'Genotype')
    if not os.path.exists(genotype_folder):
        os.mkdir(genotype_folder)

    # get subtype of reference
    subtype_str = determine_subtype(kma_out_file, reference_folder_name, 10)
    print(subtype_str)

    root_file_path = os.path.dirname(os.path.realpath(__file__))
    render_file = os.path.join(root_file_path, 'RenderTree.R')
    
    ## genotype using whole genome (Dumb method)
    if "SubtypeA" in subtype_str:
        ref_db_path = os.path.join(reference_folder_name, 'Genotype_ref', 'BlastDB', 'RSVA')
        meta_file_path = os.path.join(reference_folder_name, 'Genotype_ref', 'metadata_A.tsv')
    elif "SubtypeB" in subtype_str:
        ref_db_path = os.path.join(reference_folder_name, 'Genotype_ref', 'BlastDB', 'RSVB')
        meta_file_path = os.path.join(reference_folder_name, 'Genotype_ref', 'metadata_B.tsv')
    else:
        return

    print(f"genotype_call_whole_genome {query_file_path} {ref_db_path} {genotype_folder} {reference_folder_name} {render_file}")
    if run_eachstep['genotype_whole_genome'] is True:
        genotype_call_whole_genome(query_file_path, ref_db_path, meta_file_path, genotype_folder, reference_folder_name, render_file)

    ## genotype using G protein (Dumb method)
    #ref_db_path = os.path.join(reference_folder_name, 'Genotype_ref', 'G_genotype', 'G_subtype')
    #print(f"genotype_call_G_protein {query_file_path} {ref_db_path} {genotype_folder} {reference_folder_name} {render_file}")
    #if run_eachstep['genotype_G_gene'] is True:
    #    genotype_call_G_protein(query_file_path, ref_db_path, genotype_folder, reference_folder_name, render_file)

    ## blast query against GISAID
    blast_out_file = os.path.join(genotype_folder, 'blastn_res_gisaid.tsv')
    ref_db_path = os.path.join(reference_folder_name, 'Genotype_ref', 'GISAIDDB', 'GISAID')
    cmd = f"blastn -query {query_file_path} -db {ref_db_path} -out {blast_out_file} -outfmt 6 -num_threads 1 -evalue 1e-10 -perc_identity 90"
    print(cmd)
    if run_eachstep['blast_gisaid'] is True:
        subprocess.run(cmd, shell=True, cwd=genotype_folder)

    ## genotype using whole genome using NextClade (Trust this one)
    if "SubtypeA" in subtype_str:
        ref_db_path = os.path.join(reference_folder_name, 'NextClade', 'rsv_a')
    elif "SubtypeB" in subtype_str:
        ref_db_path = os.path.join(reference_folder_name, 'NextClade', 'rsv_b')
    else:
        return
    output_folder = os.path.join(genotype_folder, 'NextClade')
    nextclade_cmd = f"nextclade3 run --input-dataset {ref_db_path} --output-all={output_folder} {query_file_path}"
    print(nextclade_cmd)
    if run_eachstep['nextclade'] is True:
        subprocess.run(nextclade_cmd, shell=True, cwd=genotype_folder)
    

if __name__ == "__main__":
    sample_id = sys.argv[1]
    working_folder_name = sys.argv[2]
    original_read1 = sys.argv[3]
    original_read2 = sys.argv[4]
    reference_folder_name = sys.argv[5]
    star_ThreadN = sys.argv[6]
    igv_cutoff = sys.argv[7]
    
    run_eachstep = {
        "fastp":True, 
        "kma":True,
        "star_index":True,
        "star":True,
        "samtools":True,
        "igvtools":True,
        "genotype_whole_genome":True,
        "genotype_G_gene":True,
        "blast_gisaid":True,
        "nextclade":True
    }
    #print(sys.argv)

    sample_mapping(sample_id, working_folder_name, original_read1, original_read2, reference_folder_name, star_ThreadN, igv_cutoff, run_eachstep)
