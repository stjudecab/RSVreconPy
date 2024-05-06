#!/usr/bin/env python3

import os
import sys
import subprocess
from pathlib import Path
from yaml import safe_load
import argparse
import pandas as pd
import shutil



def sample_mapping(sample_id, working_folder_name, original_read1, original_read2, reference_folder_name)
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
    subprocess.run(cmd, shell=True, cwd=cur_folder)

    #print("######################\nFinish Quality trimming . \n######################")

    ##############################################################################
    ### step 2    KMA and building reference
    ##############################################################################
    print(f"\t######################\tStep 2 Identifying the best reference ... \t######################")

    kma_db = os.path.join(reference_folder_name, 'RSV_DB')
    kma_folder = os.path.join(cur_folder, 'KMA')
    if not os.path.exists(kma_folder):
        os.mkdir(kma_folder)

    kma_out = os.path.join(kma_folder, sample_id)
    log_file = os.path.join(log_cur_folder, 'kma_log.txt')

    '''
    kma_cmd = [
        'kma',
        '-ipe', read1_trim, read2_trim,
        '-o', kma_out,
        '-t_db', kma_db,
        '-1t1', 
        '-t', str(star_ThreadN),
        '-mrs','0.9',
        '-mq','50'
    ]
    '''

    kma_cmd = f"kma -ipe {read1_trim} {read2_trim} -o {kma_out} -t_db {kma_db} -1t1 -t 1 -mrs 0.9 -mq 50 &> {log_file}"
    subprocess.run(kma_cmd, shell=True, cwd=cur_folder)

    # read info from KMA results
    kma_out_file = kma_out + '.res'
    kma_df = pd.read_csv(kma_out_file, sep='\t')
    kma_df = kma_df.sort_values('Score', ascending=False)
    Selected_ref_name = kma_df.iloc[0, 0]
    #print(Selected_ref_name)

    # build reference
    original_fasta_file = os.path.join(reference_folder_name, 'Seq', Selected_ref_name + '.fasta')
    original_gff_file = os.path.join(reference_folder_name, 'Seq', Selected_ref_name + '.gff')

    ref_folder = os.path.join(cur_folder, 'reference')
    if not os.path.exists(ref_folder):
        os.mkdir(ref_folder)
    ref_output_dir = os.path.join(ref_folder, Selected_ref_name)
    new_fasta_file = os.path.join(ref_folder, Selected_ref_name + '.fasta')
    new_gff_file = os.path.join(ref_folder, Selected_ref_name + '.gff')
    shutil.copy(original_fasta_file, new_fasta_file)
    shutil.copy(original_gff_file, new_gff_file)
    '''
    star_command = [
        'STAR',
        '--runMode', 'genomeGenerate',
        '--genomeFastaFiles', new_fasta_file,
        '--genomeDir', ref_output_dir,
        '--runThreadN', str(star_ThreadN)
    ]
    '''
    log_file = os.path.join(log_cur_folder, 'star_build_genome_log.txt')
    star_command = f"STAR --runMode genomeGenerate --genomeFastaFiles {new_fasta_file} --genomeDir {ref_output_dir} --runThreadN 1 --genomeSAindexNbases 6 --sjdbGTFtagExonParentTranscript {new_gff_file}  &> {log_file}"
    subprocess.run(star_command, shell=True, cwd=cur_folder)

    ##############################################################################
    ## step 3    STAR mapping
    ##############################################################################
    print(f"\t######################\tStep 3 mapping ... \t######################")

    map_folder = os.path.join(cur_folder, 'mapping')
    if not os.path.exists(map_folder):
        os.mkdir(map_folder)
    '''
    cmd_mapping = [
        'STAR',
        '--runThreadN', str(star_ThreadN),
        '--genomeDir', ref_output_dir, 
        '--outSAMtype', 'BAM',
        'SortedByCoordinate',
        '--readFilesIn', read1_trim, read2_trim,
        '--outFileNamePrefix', map_folder + '/'
    ]
    '''

    log_file = os.path.join(log_cur_folder, 'star_mapping_log.txt')
    cmd_mapping = f"STAR --runThreadN {star_ThreadN} --genomeDir {ref_output_dir} --outSAMtype BAM SortedByCoordinate --readFilesIn {read1_trim} {read2_trim} --limitBAMsortRAM 12884901888 --outFileNamePrefix {map_folder}/ &> {log_file}"
    subprocess.run(cmd_mapping, shell=True, cwd=cur_folder)

    #print("######################\nFinish STAR mapping . \n######################")
    
    ##############################################################################
    ## step 4    index and Count using Samtools and IGVtools 
    ##############################################################################
    print(f"\t######################\tStep 4 Index and Count ... \t######################")

    '''
    cmd_index = [
        'samtools',
        'index', 'Aligned.sortedByCoord.out.bam'
    ]
    '''
    log_file = os.path.join(log_cur_folder, 'samtools_log.txt')
    cmd_index = f"samtools index Aligned.sortedByCoord.out.bam &> {log_file}"
    subprocess.run(cmd_index, shell=True, cwd=map_folder)

    '''
    cmd_count = [
        'igvtools', 
        'count',
         '-z', '5', 
         '-w', '1',
        '--bases', 
        'Aligned.sortedByCoord.out.bam',
        'alignments.cov.wig', 
        new_fasta_file
    ]
    '''
    log_file = os.path.join(log_cur_folder, 'igvtools_log.txt')
    cmd_count = f"igvtools count -z 5 -w 1 --bases Aligned.sortedByCoord.out.bam alignments.cov.wig {new_fasta_file} &> {log_file}"
    subprocess.run(cmd_count, shell=True, cwd=map_folder)

if __name__ == "__main__":
    sample_id = sys.argv[1]
    working_folder_name = sys.argv[2]
    original_read1 = sys.argv[3]
    original_read2 = sys.argv[4]
    reference_folder_name = sys.argv[5]

    sample_mapping(sample_id, working_folder_name, original_read1, original_read2, reference_folder_name)