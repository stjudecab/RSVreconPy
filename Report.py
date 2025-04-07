#!/usr/bin/env python3

import os
import sys
import argparse
from yaml import safe_load


###################################################
##      load packages and functions
###################################################

from Report_functions import generate_pdf_report, generate_html_report, generate_csv_fasta, generate_phylogenetic_tree

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

    star_ThreadN = config.get('THREAD_N', 2)
    igv_cutoff = config.get('COV_CUTOFF', 50)

    addtional_results_path = os.path.join(working_folder_name, '..', 'consensus')
    if os.path.exists(addtional_results_path):
        pass
    else:
        addtional_results_path = None

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
    generate_phylogenetic_tree(root_file_path, reference_folder_name, working_folder_name, subtype_a_names, subtype_b_names, Temp_folder_name, addtional_results_path)

    ###################################################
    ##      generate pdf report
    ###################################################

    print("Generating Pdf report...")
    generate_pdf_report(root_file_path, report, working_folder_name, mapres_folder, igv_cutoff)

    ###################################################
    ##      generate html report
    ###################################################

    print("Generating HTML report...")
    generate_html_report(root_file_path, report, working_folder_name, mapres_folder, igv_cutoff)

    ###################################################
    ##      program finished
    ###################################################

    print("Report and sequences have been generated!")
