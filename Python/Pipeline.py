#!/usr/bin/env python3

import os
import subprocess
from pathlib import Path
from yaml import safe_load

###################################################
##      load packages
###################################################

# Assuming RSV_functions is in the same directory as this script
from RSV_functions import get_sub_folders, elements_not_in_array

###################################################
##      load parameters from YAML file
###################################################
config_file = input("Enter the path to the YAML config file: ")
with open(config_file, 'r') as f:
    config = safe_load(f)

data_folder_name = config['DATA_DIR']
reference_folder_name = config['REFERENCE_DIR']
working_folder_name = config['OUTPUT_DIR']

star_ThreadN = config.get('STAR_THREAD_N', 16)
ref_use = config.get('REF_USE', 'all')

###################################################
##      input data read
###################################################
data_folder_name = os.path.normpath(data_folder_name)

sample_hash = {}

with os.scandir(data_folder_name) as entries:
    for entry in entries:
        if entry.name.endswith('.fastq') and entry.is_file():
            file = entry.name
            sample_id = file.split('_R')[0]

            if '_R1' in file:
                sample_hash[sample_id] = {1: file}
            elif '_R2' in file:
                sample_hash[sample_id][2] = file

###################################################
##      setup working dir
###################################################
working_folder_name = os.path.normpath(working_folder_name)
if not os.path.exists(working_folder_name):
    os.mkdir(working_folder_name)

###################################################
##      start process
###################################################

for sample_id in sorted(sample_hash.keys()):
    print(f"Start process {sample_id} ...")

    # skip Undetermined reads
    if "Undetermined" in sample_id:
        continue

    cur_folder = os.path.join(working_folder_name, sample_id)
    if not os.path.exists(cur_folder):
        os.mkdir(cur_folder)

    original_read1 = os.path.join(data_folder_name, sample_hash[sample_id][1])
    original_read2 = os.path.join(data_folder_name, sample_hash[sample_id][2])

    ##############################################################################
    ### step 1    QC using fastp
    ##############################################################################
    print(f"######################\tStep 1 Quality check and trimming using fastp ... \t######################")

    read1_trim = "reads_trimed_R1.fastq"
    read2_trim = "reads_trimed_R2.fastq"

    cmd = f"cd {cur_folder}; fastp -i {original_read1} -I {original_read2} -o {read1_trim} -O {read2_trim}"
    subprocess.run(cmd, shell=True)

    print("######################\nFinish Quality trimming . \n######################")

    ##############################################################################
    ## step 2    STAR mapping
    ##############################################################################
    print(f"######################\tStep 2 STAR mapping ... \t######################")

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

    for ref_db in all_ref_db:
        refdb_A = os.path.join(reference_folder_name, ref_db)
        map_A_folder = os.path.join(cur_folder, ref_db)

        if not os.path.exists(map_A_folder):
            os.mkdir(map_A_folder)

        cmd = f"cd {cur_folder}; STAR --runThreadN {star_ThreadN} --genomeDir {refdb_A} --outSAMtype BAM SortedByCoordinate --readFilesIn {read1_trim} {read2_trim} --outFileNamePrefix {map_A_folder}"
        subprocess.run(cmd, shell=True)

    print("######################\nFinish STAR mapping . \n######################")

    ##############################################################################
    ## step 3    index and Count using Samtools and IGVtools 
    ##############################################################################
    print(f"######################\tStep 3 Index and Count ... \t######################")

    for ref_db in all_ref_db:
        map_A_folder = os.path.join(cur_folder, ref_db)

        cmd_index = f"cd {map_A_folder}; samtools index Aligned.sortedByCoord.out.bam"
        subprocess.run(cmd_index, shell=True)

        cmd_count = f"cd {map_A_folder}; igvtools count -z 5 -w 1 --bases Aligned.sortedByCoord.out.bam alignments.cov.wig {refdb_A}.fasta"
        subprocess.run(cmd_count, shell=True)

    print("######################\nFinish Index and Count \n######################")
