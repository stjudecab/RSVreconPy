#!/usr/bin/env python3

###################################################
##      load packages and functions
###################################################

import os
import sys
import subprocess
from pathlib import Path
from yaml import safe_load
import argparse
import pandas as pd
import shutil
import concurrent.futures


from RSV_functions import get_sub_folders, elements_not_in_array, check_tool_availability, detect_sequencing_files
from Report_functions import generate_pdf_report, generate_html_report, generate_csv_fasta, generate_phylogenetic_tree


###################################################
##      load parameters from YAML file
###################################################

# Create command-line argument parser
parser = argparse.ArgumentParser(description='Run pipeline to map NGS reads against reference genomes.sortedByCoord')
parser.add_argument('yaml_file', help='Path to the YAML file containing the setting.')

# Parse command-line arguments
args = parser.parse_args()

try:
    with open(args.yaml_file, 'r') as f:
        config = safe_load(f)
except:
    print('Can not open this file!')
    sys.exit()

data_folder_name = config['DATA_DIR']
reference_folder_name = config['REFERENCE_DIR']
working_folder_name = config['OUTPUT_DIR']

star_ThreadN = config.get('THREAD_N', 2)
igv_cutoff = config.get('COV_CUTOFF', 50)
igv_cutoff_low = config.get('COV_CUTOFF_LOW', 10)
MAX_CONCURRENT_JOBS = config.get('MAX_CONCURRENT_JOBS', 2)
additional_results_path = config.get('RSV_NEXT_PIPE_RES', None)

# Only check path existence if it's not None
if additional_results_path is not None and os.path.exists(additional_results_path):
    pass  # Path exists and is valid
else:
    additional_results_path = None

tool = config.get('TOOL', 'BWA')

if tool not in ['BWA']:
    print(tool + " is not a supported aligner! Please use BWA!\n")
    sys.exit()

###################################################
##      Check tools' availability
###################################################
print(f"######################\tStep 1 Checking tools' availability ... \t######################")
check_tool_availability_res = 0
# Check for IGVtools
check_tool_availability_res += check_tool_availability("igvtools")
# Check for Samtools
check_tool_availability_res += check_tool_availability("samtools")
# Check for fastp
check_tool_availability_res += check_tool_availability("fastp")
# Check for aligners
check_tool_availability_res += check_tool_availability("bwa")
# Check for kma
check_tool_availability_res += check_tool_availability('kma')
# Check for blastn
check_tool_availability_res += check_tool_availability('blastn')
# Check for mafft
check_tool_availability_res += check_tool_availability('mafft')

if check_tool_availability_res > 0:
    sys.exit()

#print(f"######################\tStep 0  PASSED ...                      \t######################")

print(f"######################\tStep 2  process each sample ...         \t######################")

###################################################
##      read input data
###################################################
data_folder_name = os.path.normpath(data_folder_name)
sample_dict = detect_sequencing_files(data_folder_name)

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
##      download latest NextClade reference database
###################################################
db_path = os.path.join(reference_folder_name, 'NextClade','rsv_a')
nextclade_cmd = f"nextclade3 dataset get --name rsv_a --output-dir {db_path}"
#print(nextclade_cmd)
subprocess.run(nextclade_cmd, shell=True)

db_path = os.path.join(reference_folder_name, 'NextClade','rsv_b')
nextclade_cmd = f"nextclade3 dataset get --name rsv_b --output-dir {db_path}"
#print(nextclade_cmd)
subprocess.run(nextclade_cmd, shell=True)

print(f"Latest database downloaded from RSV A and RSV B\n")

###################################################
##      start process with concurrent execution
###################################################
def run_mapping(sample_id, mapres_folder_name, original_read1, original_read2, reference_folder_name, star_ThreadN, igv_cutoff, igv_cutoff_low):
    """Run mapping command for a single sample"""
    root_file_path = os.path.dirname(os.path.realpath(__file__))
    log_file = os.path.join(log_folder_name, sample_id + '.log')
    errlog_file = os.path.join(log_folder_name, sample_id + '.err.log')
    job_command = f"python {root_file_path}/Mapping.py {sample_id} '{mapres_folder_name}' '{original_read1}' '{original_read2}' '{reference_folder_name}' {star_ThreadN} {igv_cutoff} {igv_cutoff_low}"
    
    # Redirect output to log files
    with open(log_file, 'w') as log, open(errlog_file, 'w') as errlog:
        process = subprocess.run(job_command, shell=True, stdout=log, stderr=errlog)
    
    return sample_id, process.returncode

# Set maximum concurrent jobs
completed = 0
total_samples = len(sample_dict)

with concurrent.futures.ThreadPoolExecutor(max_workers=MAX_CONCURRENT_JOBS) as executor:
    # Dictionary to track futures
    future_to_sample = {}
    
    # Submit initial batch of jobs
    initial_samples = list(sorted(sample_dict.keys()))[:MAX_CONCURRENT_JOBS]
    remaining_samples = list(sorted(sample_dict.keys()))[MAX_CONCURRENT_JOBS:]
    
    for sample_id in initial_samples:
        if "Undetermined" in sample_id:
            continue
            
        original_read1 = os.path.join(data_folder_name, sample_dict[sample_id][0])
        original_read2 = os.path.join(data_folder_name, sample_dict[sample_id][1])
        
        future = executor.submit(
            run_mapping,
            sample_id,
            mapres_folder_name,
            original_read1,
            original_read2,
            reference_folder_name,
            star_ThreadN,
            igv_cutoff
        )
        future_to_sample[future] = sample_id
        print(f"Submitted job for sample: {sample_id}")
    
    # Process remaining samples as jobs complete
    while future_to_sample:
        # Wait for the next job to complete
        done, _ = concurrent.futures.wait(
            future_to_sample.keys(),
            return_when=concurrent.futures.FIRST_COMPLETED
        )
        
        # Process completed jobs
        for future in done:
            sample_id = future_to_sample.pop(future)
            try:
                sample_id, returncode = future.result()
                completed += 1
                if returncode == 0:
                    print(f"Completed ({completed}/{total_samples}): {sample_id} - Success")
                else:
                    print(f"Completed ({completed}/{total_samples}): {sample_id} - Failed (return code: {returncode})")
            except Exception as e:
                completed += 1
                print(f"Error processing {sample_id}: {str(e)}")
            
            # Submit a new job if there are remaining samples
            if remaining_samples:
                next_sample = remaining_samples.pop(0)
                if "Undetermined" in next_sample:
                    continue
                    
                original_read1 = os.path.join(data_folder_name, sample_dict[next_sample][0])
                original_read2 = os.path.join(data_folder_name, sample_dict[next_sample][1])
                
                future = executor.submit(
                    run_mapping,
                    next_sample,
                    mapres_folder_name,
                    original_read1,
                    original_read2,
                    reference_folder_name,
                    star_ThreadN,
                    igv_cutoff
                )
                future_to_sample[future] = next_sample
                print(f"Submitted job for sample: {next_sample}")

print(f"All mapping jobs are completed! Total: {completed}/{total_samples}")


###################################################
##      Aggregate results and generate reports
###################################################

print(f"######################\tStep 3  Aggregate results ...           \t######################")

root_file_path = os.path.dirname(os.path.realpath(__file__))
Temp_folder_name = os.path.join(working_folder_name, 'Temp')
if not os.path.exists(Temp_folder_name):
    os.mkdir(Temp_folder_name)

report = os.path.join(working_folder_name, "Report.csv")
sequence_file = os.path.join(working_folder_name, "Sequence.fasta")
sequence_file_a = os.path.join(working_folder_name, "Sequence_A.fasta")
sequence_file_b = os.path.join(working_folder_name, "Sequence_B.fasta")
mapres_folder = os.path.join(working_folder_name, "Mapping")

print("Aggregating results ...")
subtype_a_names, subtype_b_names = generate_csv_fasta(report, sequence_file, reference_folder_name, mapres_folder, sequence_file_a, sequence_file_b, igv_cutoff)

# generate phylogenetic tree with reference

print("Generating tree ...")
generate_phylogenetic_tree(root_file_path, reference_folder_name, working_folder_name, subtype_a_names, subtype_b_names, Temp_folder_name, additional_results_path)

# generate pdf report

print("Generating Pdf report...")
generate_pdf_report(root_file_path, report, working_folder_name, mapres_folder, igv_cutoff)

# generate html report

print("Generating HTML report...")
generate_html_report(root_file_path, report, working_folder_name, mapres_folder, igv_cutoff)

print("Report and sequences have been generated!")

print(f"######################\tSuccessfully completed!                 \t######################")






