#!/usr/bin/env python3

import os
import json
import re
#from Bio.Seq import Seq
#from Bio.SeqRecord import SeqRecord
#from Bio import SeqIO
from yaml import safe_load

###################################################
##      load packages
###################################################

# Assuming RSV_functions is in the same directory as this script
from RSV_functions import get_sub_folders, elements_not_in_array, pct_sum, determine_subtype, processIGV

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
##      setup working dir
###################################################
working_folder_name = os.path.normpath(working_folder_name)
if not os.path.exists(working_folder_name):
    os.mkdir(working_folder_name)

###################################################
##      start process
###################################################

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
CSV_header_lv_1 = ',before_filtering, before_filtering, before_filtering, after_filtering, after_filtering, after_filtering, QC rate,'
CSV_header_lv_2 = 'Sample name,total_reads, q20_rate, q30_rate, total_reads, q20_rate, q30_rate, QC rate,'
for ref_db in all_ref_db:
    CSV_header_lv_1 += f"Mapping {ref_db},Mapping {ref_db},Mapping {ref_db},Mapping {ref_db},"
    CSV_header_lv_2 += 'Uniquely mapped reads %, MULTI-MAPPING READS %, UNMAPPED READS%, CHIMERIC READS%,'
CSV_header_lv_1 += "Subtype suggestion"
CSV_header_lv_2 += "Subtype suggestion"

report = os.path.join(working_folder_name, "Report.csv")
with open(report, 'w') as out:
    out.write(CSV_header_lv_1 + "\n")
    out.write(CSV_header_lv_2 + "\n")

sequence_file = os.path.join(working_folder_name, "Sequence.fasta")
with open(sequence_file, 'w') as fasta:

    # collect information from each sample
    Sample_folders = [f for f in os.listdir(working_folder_name) if os.path.isdir(os.path.join(working_folder_name, f))]
    for cur_folder in Sample_folders:
        if os.path.isfile(cur_folder):
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
            uniquely_mapped_reads_hash[all_ref_db[i]] = result['Uniquelymappedreads']

        subtype_str = determine_subtype(uniquely_mapped_reads_hash)

        sample = cur_folder
        sample = sample.split("/")[-1]

        with open(report, 'a') as out:
            out.write(f"{sample},{QC_str}{map_str}{subtype_str}\n")

        # Step 4: assemble genome sequence
        wig_file = os.path.join(working_folder_name, cur_folder, subtype_str, 'alignments.cov.wig')
        genome_sequence = processIGV(wig_file)

        fasta.write('>' + sample + "\n" + genome_sequence + "\n")

###################################################
##      program finished
###################################################

print("Report and sequences have been generated!")
