#!/usr/bin/env python3

import os
import sys
import subprocess
from pathlib import Path
from yaml import safe_load
import argparse
import pandas as pd
import shutil


###################################################
##      load packages
###################################################

# Assuming RSV_functions is in the same directory as this script
from RSV_functions import get_sub_folders, elements_not_in_array, check_tool_availability


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
tool = config.get('TOOL', 'STAR')

if tool not in ['STAR']:
    print(tool + " is not a supported aligner! Please use STAR!\n")
    sys.exit()


###################################################
##      Check tools' availability
###################################################
print(f"######################\tStep 0 Checking tools' availability ... \t######################")
check_tool_availability_res = 0
# Check for IGVtools
check_tool_availability_res += check_tool_availability("igvtools")
# Check for Samtools
check_tool_availability_res += check_tool_availability("samtools")
# Check for fastp
check_tool_availability_res += check_tool_availability("fastp")
# Check for aligners
check_tool_availability_res += check_tool_availability(tool)
# Check for kma
check_tool_availability_res += check_tool_availability('kma')
# Check for blastn
check_tool_availability_res += check_tool_availability('blastn')


if check_tool_availability_res > 0:
    sys.exit()

print(f"######################\tStep 0  PASSED ...                      \t######################")


###################################################
##      input data read
###################################################
data_folder_name = os.path.normpath(data_folder_name)

sample_dict = {}

with os.scandir(data_folder_name) as entries:
    for entry in entries:
        if entry.name.endswith('.fastq.gz') and entry.is_file():
            file = entry.name
            sample_id = file.split('_R')[0]

            if sample_id not in sample_dict:
                sample_dict[sample_id] = ['R1','R2']
                
            if '_R1' in file:
                sample_dict[sample_id][0] = file
            elif '_R2' in file:
                sample_dict[sample_id][1] = file


###################################################
##      setup working dir
###################################################
working_folder_name = os.path.normpath(working_folder_name)
if not os.path.exists(working_folder_name):
    os.mkdir(working_folder_name)

log_folder_name = os.path.join(working_folder_name, 'log')
if not os.path.exists(log_folder_name):
    os.mkdir(log_folder_name)

mapres_folder_name = os.path.join(working_folder_name, 'Mapping')
if not os.path.exists(mapres_folder_name):
    os.mkdir(mapres_folder_name)

###################################################
##      start process
###################################################

for sample_id in sorted(sample_dict.keys()):
    print(f"Start process {sample_id} ...")

    # skip Undetermined reads
    if "Undetermined" in sample_id:
        continue

    sample_id = sample_id
    working_folder_name = working_folder_name
    original_read1 = os.path.join(data_folder_name, sample_dict[sample_id][0])
    original_read2 = os.path.join(data_folder_name, sample_dict[sample_id][1])
    reference_folder_name = reference_folder_name
    #cmd = f"python Mapping.py {sample_id} '{working_folder_name}' '{original_read1}' '{original_read2}' '{reference_folder_name}' {star_ThreadN}"

    root_file_path = os.path.dirname(os.path.realpath(__file__))
    log_file = os.path.join(log_folder_name, sample_id + '.log')
    errlog_file = os.path.join(log_folder_name, sample_id + '.errlog')
    bsub_command = f"bsub -q priority -P CAB -J STAR_Mapping_{sample_id} -M 2G -n {star_ThreadN} -oo {log_file} -eo {errlog_file} \"python {root_file_path}/Mapping.py {sample_id} \'{mapres_folder_name}\' \'{original_read1}\' \'{original_read2}\' \'{reference_folder_name}\' {star_ThreadN} {igv_cutoff}\""
    #print(bsub_command)
    subprocess.run(bsub_command, shell = True)

    print(f"Mapping jobs submitted for {sample_id}!")

###################################################
##      program finished
###################################################


print(f"All mapping jobs are submitted!")
