import os
import glob
import json
import re
import sys
import shutil
import subprocess
import fileinput
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.colors as mcolors
from PIL import Image as PILImage
from reportlab.lib.styles import getSampleStyleSheet, TA_LEFT, ParagraphStyle
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Flowable, Table, Image, TableStyle, PageBreak, ActionFlowable
from reportlab.lib.colors import black
from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm
from reportlab.lib import colors
from reportlab.lib.units import inch
from RSV_functions import parse_gff, find_gene_at_position, find_gff_files_in_path, get_genotype_res, get_version, get_sub_folders, elements_not_in_array, pct_sum, determine_subtype, processIGV, extract_gene_seq, detect_F_mutation, extract_key_residue_fgene, extract_gene_covarage, translate_nt_to_aa, get_version
from SNP import SNP_calling
import seaborn as sns
import base64

def image_to_base64(image_path):
    with open(image_path, "rb") as image_file:
        encoded_string = base64.b64encode(image_file.read()).decode('utf-8')
    return encoded_string

class Line(Flowable):
    """Line flowable --- draws a line in a flowable"""

    #----------------------------------------------------------------------
    def __init__(self, width, height=0):
        Flowable.__init__(self)
        self.width = width
        self.height = height

    #----------------------------------------------------------------------
    def draw(self):
        """draw the line"""
        self.canv.line(0, self.height, self.width, self.height)

def rgba_to_hex(rgba):
    return '#' + ''.join(f'{int(x*255):02x}' for x in rgba[:3])

# Create a color gradient function
def color_gradient(value):
    """Return a color from green to red based on the input value (0 to 1)."""
    return colors.Color(1 - value, value, 0)

def color_gradient_matplotlib(value):
    colors = ["#ff4d4d", "#4dff4d"]
    cmap = mcolors.LinearSegmentedColormap.from_list("", colors)
    return cmap(value / 100)

def color_gradient_blue(value):
    """Return a color from green to red based on the input value (0 to 1)."""
    return colors.Color(1 - value, 1 - value, 1)

def int_to_gradient_color(value):
    """
    Convert an integer between 0 and 100 to a gradient color between gray and green.

    Parameters:
    value (int): The input integer between 0 and 100.

    Returns:
    str: The hexadecimal color code representing the gradient color.
    """
    if not (0 <= value <= 100):
        raise ValueError("Input value must be between 0 and 100")

    # Define the RGB values for gray and green
    gray_rgb = (211, 211, 211)  # Light Gray color
    green_rgb = (144, 238, 144)     # Light Green color

    # Calculate the interpolation factor
    factor = value / 100.0

    # Interpolate the RGB values
    r = int(green_rgb[0] + factor * (gray_rgb[0] - green_rgb[0]))
    g = int(green_rgb[1] + factor * (gray_rgb[1] - green_rgb[1]))
    b = int(green_rgb[2] + factor * (gray_rgb[2] - green_rgb[2]))

    # Convert the RGB values to a hexadecimal color code
    hex_color = f'#{r:02x}{g:02x}{b:02x}'

    return hex_color

def percentage_to_number(percentage):
    return float(percentage.strip('%')) / 100

def float_to_percentage(value):
    if isinstance(value, np.float64):
        return "{:.2%}".format(value)
    else:
        return "{:.2%}".format(value)

def count_greater_than_n(array, n):
    return len([x for x in array if x > n])

def get_best_blast_hit(blast_res):
    blast_df = pd.read_csv(blast_res, sep='\t', header=None)
    blast_df.columns = ['Query', 'Subject', 'Pct_identity','Alignment_length','Number_of_mismatches','Number_of_gap','Start_in_query','End_in_query','Start_in_subject','End_in_subject','E_value','Bit_score']
    blast_df = blast_df.sort_values(by=['Bit_score'], ascending = False)

    Selected_ref_name = blast_df.loc[0, 'Subject']
    return Selected_ref_name, blast_df.loc[0, 'Alignment_length'], blast_df.loc[0, 'Pct_identity']

def make_coverage_heatmap(df, png_file):
    number_of_rows = df.shape[0]
    height = 2 + number_of_rows * 0.15
    df.columns = df.columns.str.replace('_cov', '')

    # make a cmap
    colors = ['#D3D3D3', '#d1fbd8', '#a3ee5b', '#40ba0f']
    n_bins = [0, 50, 500, 1000, 50000]  # Include np.inf for values > 5000

    # Create a custom colormap
    cmap_name = 'custom_cmap'
    cmap = LinearSegmentedColormap.from_list(cmap_name, colors, N=len(n_bins) - 1)
    norm = BoundaryNorm(boundaries=n_bins, ncolors=len(n_bins) - 1)

    plt.figure(figsize=(9, height))  # Optional: Set the size of the figure
    plt.rcParams['font.size'] = 7
    ax = sns.heatmap(df, annot=False, cmap=cmap, norm=norm, cbar_kws={"aspect": 20, "shrink": 0.5}, linewidths=0.2, linecolor='black')

    # Rotate the x-tick labels
    #ax.set_ylabel('')
    plt.yticks(rotation=0)
    plt.xticks(rotation=0)

    xticks = ax.get_xticklabels()
    xticks[6].set_fontweight('bold')
    xticks[7].set_fontweight('bold')
    ax.set_xticklabels(xticks)


    # Draw custom grid lines for a specific column
    for y in range(len(df)):
        plt.plot([6, 6], [y - 1, y + 1], color='black', linewidth=1)  # Adjust column index and line width
        plt.plot([8, 8], [y - 1, y + 1], color='black', linewidth=1)  # Adjust column index and line width

    # save figure
    plt.savefig(png_file, dpi = 300, bbox_inches='tight')
    plt.close()


def array_to_html_table(array, header, color=None, table_id=None, table_class=None):
    """
    Convert a 2D array into an HTML table with optional ID and class.

    Parameters:
    array (list of lists): The 2D array to convert.
    table_id (str): Optional ID for the table.
    table_class (str): Optional class for the table.

    Returns:
    str: The HTML string representing the table.
    """
    # Start the HTML table with optional ID and class
    id_attr = f' id="{table_id}"' if table_id else ''
    class_attr = f' class="{table_class}"' if table_class else ''
    html = f'<table border="1"{id_attr}{class_attr}>\n'

    
    # make header
    html += '  <thead>\n'
    html += '  <tr>\n'
    # Iterate over the cells in the row
    for cell in header:
        html += f'    <td>{cell}</td>\n'
    html += '  </tr>\n'
    html += '  </thead>\n'

    # Iterate over the rows in the array
    html += '  <tbody>\n'
    for i, row in enumerate(array):
        html += '  <tr>\n'
        # Iterate over the cells in the row
        for j, cell in enumerate(row):
            # Determine the background color for the cell
            bg_color = color[i][j] if color and i < len(color) and j < len(color[i]) else ''
            style_attr = f' style="background-color: {bg_color};"' if bg_color else ''
            html += f'    <td{style_attr}>{cell}</td>\n'
        html += '  </tr>\n'
    html += '  </tbody>\n'
    # End the HTML table
    html += '</table>'

    return html

def generate_empty_2d_array(array):
    """
    Generate a 2D array with the same dimensions as the input array, filled with empty strings.

    Parameters:
    array (list of lists): The input 2D array.

    Returns:
    list of lists: A new 2D array with the same dimensions, filled with empty strings.
    """
    # Determine the number of rows and columns in the input array
    num_rows = len(array)
    num_cols = len(array[0]) if num_rows > 0 else 0

    # Create a new 2D array with the same dimensions, filled with empty strings
    empty_array = [["" for _ in range(num_cols)] for _ in range(num_rows)]

    return empty_array

##########################################
# generate CSV report
##########################################

def generate_csv_fasta(report, sequence_file, reference_folder_name, working_folder_name, sequence_file_a, sequence_file_b, igv_cutoff):
    # get version info
    current_script_path = os.path.dirname(os.path.abspath(__file__))
    version_file_path = os.path.join(current_script_path, 'version.txt')
    current_version = get_version(version_file_path)

    # create output CSV
    CSV_header = f"Pipeline verions: {current_version}\n"
    CSV_header += "Sample name,before_filtering_total_reads, before_filtering_q20_rate, before_filtering_q30_rate, after_filtering_total_reads, after_filtering_q20_rate, after_filtering_q30_rate, QC rate,"
    CSV_header += "Uniquely mapped reads %,MULTI-MAPPING READS %,UNMAPPED READS%,CHIMERIC READS%,"
    CSV_header += "Subtype,reference_accession,ref_subtype,"
    CSV_header += "F protein mutations,"
    CSV_header += "NS1_cov,NS2_cov,N_cov,P_cov,M_cov,SH_cov,G_cov,F_cov,M2-1_cov,M2-2_cov,L_cov,"
    CSV_header += "G-protein Clade(NextClade),Whole Genome Clade(NextClade),G-protein Clade(Blast),Whole Genome Clade(Blast)"

    subtype_a_names = []
    subtype_b_names = []

    F_protein_dict_A = {}
    F_protein_dict_B = {}

    F_sequence_A = sequence_file.replace('Sequence.fasta', 'F_protein_seq_A.fasta')
    F_sequence_B = sequence_file.replace('Sequence.fasta', 'F_protein_seq_B.fasta')

    with open(report, 'w') as out, open(F_sequence_A, 'w') as f_protein_a, open(F_sequence_B, 'w') as f_protein_b:
        out.write(CSV_header + "\n")

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
                QC_str = f"{before['total_reads']},{before['q20_rate']},{before['q30_rate']},{after['total_reads']},{after['q20_rate']},{after['q30_rate']},{rate}"

                map_str = f"{cur_map_hash['Uniquelymappedreads']},{pct_sum(cur_map_hash['ofreadsmappedtomultipleloci'], cur_map_hash['ofreadsmappedtotoomanyloci'])},{pct_sum(cur_map_hash['ofreadsunmappedtoomanymismatches'], cur_map_hash['ofreadsunmappedtooshort'], cur_map_hash['ofreadsunmappedother'])},{cur_map_hash['ofchimericreads']}"

                # reference subtype
                genotype_file = os.path.join(working_folder_name, cur_folder, 'Genotype', 'Genotype.txt')
                if os.path.isfile(genotype_file):
                    subtype_str, blast_pct_identity, blast_alignment_length, starin = get_genotype_res(genotype_file)
                    ##### fix later
                    if int(blast_alignment_length) < 1500:
                        subtype_str = 'Not RSV'
                    else:
                        subtype_str = subtype_str[0]
                else:
                    subtype_str = 'Not RSV'

                sample = cur_folder.split("/")[-1]

                kma_out_file = os.path.join(working_folder_name, cur_folder, 'KMA', cur_folder + '.res')
                ref_str = determine_subtype(kma_out_file, reference_folder_name, 10)
                ref_str = ref_str.split(",")

                out.write(f"{sample},{QC_str},{map_str},{subtype_str},{ref_str[1]},{ref_str[2]}")

                # Step 4: get assembled genome sequence
                genome_sequence_file = os.path.join(working_folder_name, cur_folder, 'sequence.fasta')
                genome_sequence = ''
                with open(genome_sequence_file, 'r') as read_seq:
                    for line in read_seq:
                        if line.startswith('>'):
                            pass
                        else:
                            genome_sequence += line.strip()   

                with open(sequence_file, 'a') as fasta, open(sequence_file_a, 'a') as fasta_a, open(sequence_file_b, 'a') as fasta_b:
                    if subtype_str == 'Not RSV':
                        sequence_file_error = sequence_file + '.err'
                        with open(sequence_file_error, 'a') as fasta_err:
                            fasta_err.write('>' + sample + "\n" + genome_sequence + "\n")
                    elif subtype_str == 'A':
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
                    seq_name, assembled_f_gene_sequence = extract_gene_seq(genome_sequence_file, reference_sequence_file, gff_file, 'CDS_8')
                    assembled_f_protein_sequence = translate_nt_to_aa(assembled_f_gene_sequence)
                    #print(assembled_f_protein_sequence)
                    f_protein_a.write(f">{cur_folder}\n{assembled_f_protein_sequence}\n")
                    F_mutation = detect_F_mutation(assembled_f_protein_sequence, mutation_file)
                    mutations_str = [f"{ele[0]}({ele[1]})" for ele in F_mutation]
                    F_protein_dict_A[cur_folder] = mutations_str
                    mutations_str = "|".join(mutations_str)
                else:
                    mutation_file = os.path.join(reference_folder_name, 'Genotype_ref','RSV_B_F_Mutation.csv')
                    seq_name, assembled_f_gene_sequence = extract_gene_seq(genome_sequence_file, reference_sequence_file, gff_file, 'CDS_8')
                    assembled_f_protein_sequence = translate_nt_to_aa(assembled_f_gene_sequence)
                    #print(assembled_f_protein_sequence)
                    f_protein_b.write(f">{cur_folder}\n{assembled_f_protein_sequence}\n")
                    F_mutation = detect_F_mutation(assembled_f_protein_sequence, mutation_file)
                    mutations_str = [f"{ele[0]}({ele[1]})" for ele in F_mutation]
                    F_protein_dict_B[cur_folder] = mutations_str
                    mutations_str = "|".join(mutations_str)

                # step 6: calculate coverage for each gene segment
                wig_file = os.path.join(working_folder_name, cur_folder, 'mapping', 'alignments.cov.wig')
                coverage_by_gene = extract_gene_covarage(wig_file, gff_file)
                rsvgenes = ['CDS_1','CDS_2','CDS_3','CDS_4','CDS_5','CDS_6','CDS_7','CDS_8','CDS_9','CDS_10','CDS_11']
                coverages = [coverage_by_gene[gene] for gene in rsvgenes]
                coverages_str = ','.join(str(x) for x in coverages)

                if mutations_str == '' and coverage_by_gene['CDS_8'] < 20:
                    mutations_str = 'Low coverage of F gene'

                out.write(f",{mutations_str},{coverages_str}")

                # Step 7: Fetch genotype information from NextClade
                if subtype_str != 'Not RSV':
                    nextclade_file = os.path.join(working_folder_name, cur_folder, 'Genotype', 'NextClade', 'nextclade.tsv')
                    df = pd.read_csv(nextclade_file, index_col=0, header=0, sep='\t')
                    G_subtype_str = df.iloc[0,2]
                    W_subtype_str = df.iloc[0,1]
                else:
                    G_subtype_str = 'Not RSV'
                    W_subtype_str = 'Not RSV'
                
                out.write(f",{G_subtype_str},{W_subtype_str}")

                # step 8: Fetch G protein genotype from blast
                if subtype_str != 'Not RSV':
                    genotype_g_file = os.path.join(working_folder_name, cur_folder, 'Genotype', 'Genotype_G.txt')
                    if os.path.isfile(genotype_g_file):
                        G_subtype_str, blast_pct_identity, blast_alignment_length, starin = get_genotype_res(genotype_g_file)
                    else:
                        G_subtype_str = 'Not RSV'
                else:
                    G_subtype_str = 'Not RSV'
                out.write(f",{G_subtype_str}")

                # step 9: Fetch whole genome genotype from blast
                if subtype_str != 'Not RSV':
                    genotype_file = os.path.join(working_folder_name, cur_folder, 'Genotype', 'Genotype.txt')
                    if os.path.isfile(genotype_file):
                        W_subtype_str, blast_pct_identity, blast_alignment_length, starin = get_genotype_res(genotype_file)
                        ##### fix later
                        if int(blast_alignment_length) < 1500:
                            W_subtype_str = 'Not RSV'
                    else:
                        W_subtype_str = 'Not RSV'
                else:
                    W_subtype_str = 'Not RSV'
                out.write(f",{W_subtype_str}\n")

            else:
                sample = cur_folder.split("/")[-1]
                QC_str = f"{before['total_reads']},{before['q20_rate']},{before['q30_rate']},{after['total_reads']},{after['q20_rate']},{after['q30_rate']},{rate},"
                map_str = f"0,0,100,0,"
                subtype_str = 'Not RSV'
                out.write(f"{sample},{QC_str}{map_str}{subtype_str},NA,Not RSV,,0,0,0,0,0,0,0,0,0,0,0,{subtype_str},{subtype_str},{subtype_str},{subtype_str}\n")

    # make F protein CSV report
    F_report_A = report.replace('Report.csv', 'F_mutation_A_report.csv')
    F_report_B = report.replace('Report.csv', 'F_mutation_B_report.csv')

    F_sequence_key_position_A = sequence_file.replace('Sequence.fasta', 'F_protein_key_position_A.csv')
    F_sequence_key_position_B = sequence_file.replace('Sequence.fasta', 'F_protein_key_position_B.csv')
    extract_key_residue_fgene(F_sequence_A, os.path.join(reference_folder_name, 'Genotype_ref','RSV_A_F_Mutation.csv'), F_sequence_key_position_A)
    extract_key_residue_fgene(F_sequence_B, os.path.join(reference_folder_name, 'Genotype_ref','RSV_B_F_Mutation.csv'), F_sequence_key_position_B)

    with open(F_report_A, 'w') as f_a_out:
        # write version info
        version_info = f"Pipeline verions: {current_version}\n"
        f_a_out.write(version_info)

        # make header
        mutation_file = os.path.join(reference_folder_name, 'Genotype_ref','RSV_A_F_Mutation.csv')
        csv_header = ''
        csv_pos = []
        with open(mutation_file, 'r') as file:
            for line in file:
                csv_header += ',' + re.sub(',','',line.strip())
                parts = line.strip().split(',')
                csv_pos.append(parts[1])

        f_a_out.write(csv_header + '\n')

        # each line for a sample
        for key, value in F_protein_dict_A.items():
            csv_line = key
            result_list = [''] * len(csv_pos)
            for mutation in value:
                match = re.search(r'\d+', mutation)  # Find the first sequence of digits
                if match:
                    number = match.group()
                    index = csv_pos.index(number)
                    mutation_str = mutation.replace('(Reported)','')
                    result_list[index] = mutation_str

            csv_line = csv_line + ',' + ','.join(result_list)
            f_a_out.write(csv_line + '\n')


    with open(F_report_B, 'w') as f_b_out:
        # write version info
        version_info = f"Pipeline verions: {current_version}\n"
        f_b_out.write(version_info)

        # make header
        mutation_file = os.path.join(reference_folder_name, 'Genotype_ref','RSV_B_F_Mutation.csv')
        csv_header = ''
        csv_pos = []
        with open(mutation_file, 'r') as file:
            for line in file:
                csv_header += ',' + re.sub(',','',line.strip())
                parts = line.strip().split(',')
                csv_pos.append(parts[1])

        f_b_out.write(csv_header + '\n')

        # each line for a sample
        for key, value in F_protein_dict_B.items():
            csv_line = key
            result_list = [''] * len(csv_pos)
            for mutation in value:
                match = re.search(r'\d+', mutation)  # Find the first sequence of digits
                if match:
                    number = match.group()
                    index = csv_pos.index(number)
                    mutation_str = mutation.replace('(Reported)','')
                    result_list[index] = mutation_str

            csv_line = csv_line + ',' + ','.join(result_list)
            f_b_out.write(csv_line + '\n')

    return subtype_a_names, subtype_b_names

##########################################
# generate trees
##########################################

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

##########################################
# generate pdf report
##########################################

def generate_pdf_report(csv_file, working_folder, mapres_folder, igv_cutoff):
    # get version info
    current_script_path = os.path.dirname(os.path.abspath(__file__))
    version_file_path = os.path.join(current_script_path, 'version.txt')
    current_version = get_version(version_file_path)

    line_width = 440

    # make a temp folder under working folder
    temp_folder = os.path.join(working_folder, 'Temp')
    if not os.path.exists(temp_folder):
        # Create the folder
        os.makedirs(temp_folder)
        print(f"Folder '{temp_folder}' created successfully for temporary files.")
    else:
        print(f"Using existing folder '{temp_folder}' for temporary files.")

    # resource path and file path
    file_path = os.path.dirname(os.path.realpath(__file__))
    Logo = os.path.join(file_path, 'Resource','CAB.png')
    pdf_report = os.path.join(working_folder, "Report.pdf")

    # create document
    doc = SimpleDocTemplate(pdf_report)
    styles = getSampleStyleSheet()
    elements = []

    # title
    title_style = styles['Title']
    title_style.fontSize = 20
    title_style.alignment = TA_LEFT
    title = Paragraph("Detection of RSV from clinical samples", title_style)
    elements.append(title)

    # logo
    pil_img = PILImage.open(Logo)
    original_width, original_height = pil_img.size
    new_height = 60
    new_width = original_width * (new_height / original_height)
    img_logo = Image(Logo, width=new_width, height=new_height)
    img_logo.hAlign = 'LEFT'
    elements.append(img_logo)

    # divide line
    line = Line(line_width)  # width of line
    elements.append(line)

    ######################################################################
    # section 1, Summary
    ######################################################################
    subtitle_summary_section = Paragraph("Summary", styles['Heading2'])
    elements.append(subtitle_summary_section)

    #spacer = Spacer(1, 12)  # width and height in points
    #elements.append(spacer)

    table_info_text = f"Pipeline Verions: {current_version}<br/> Subtypes of each reference are highlighted in different colors: <font color='red'>Subtype A</font> and <font color='blue'>Subtype B</font>"
    table_info_text += "<br/> Genotype calling is based on <a href='https://nextstrain.org/rsv/a/genome'><b>Nextstrain</b> and <b>Nextclade3</b>, Data updated 2024-05-21</a><br/> A clade lable with * indicates the clade of best blast hit strain due to negative results from NextClade3 <br/>Click the sample name to jump to the detail section<br/> "
    paragraph = Paragraph(table_info_text, styles['BodyText'])
    elements.append(paragraph)
    spacer = Spacer(1, 6)  # width and height in points
    elements.append(spacer)

    # load table for summary section
    df = pd.read_csv(csv_file, skiprows = 1)
    # Format the column as percentages
    df.iloc[:,7] = df.iloc[:,7].apply(lambda x: '{:.2%}'.format(x/100)) # QC rate
    #df.iloc[:,8] = df.iloc[:,8].apply(lambda x: '{:.2%}'.format(x/100))
    #df.iloc[:,12] = df.iloc[:,12].apply(lambda x: '{:.2%}'.format(x/100))

    df = df.iloc[:, [0, 7, 8, 12,13,14, 27, 28, 29, 30, 22]] # name, QC rate, mapping rate, subtype, reference_accession, ref_subtype, Gtype, Wtype, Gtype_blast, Wtype_blast, G_cov
    data = df.values.tolist()

    custom_style = ParagraphStyle(
        name='CustomStyle',
        parent=styles['BodyText'],
        fontSize=9  # Set the font size here
    )
    
    # make figure for mapping rate for each sample
    for row in range(0, len(data)):
        cur_sample_name = data[row][0]
        cur_sample_map = data[row][2]
        cur_sample_subtype = data[row][3]
        ref_id = data[row][4]
        ref_type = data[row][5]

        # data
        color_dict = {'SubtypeA': 'red', 'SubtypeB': 'blue', 'Not RSV': 'gray'}
        x = [f"{ref_id}"]
        y = [cur_sample_map]
        my_color = color_gradient_matplotlib(cur_sample_map)
        mapping_colors = [my_color]

        # plot
        plt.figure(figsize=(7, 0.3))
        plt.barh(x, y, color=mapping_colors)
        plt.xlim(0, 100)

        # hide tick, label and borders
        plt.xlabel('')
        plt.ylabel('')
        plt.title('')
        plt.xticks([])
        plt.yticks(color = color_dict[ref_type])
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['bottom'].set_visible(False)
        plt.gca().spines['left'].set_visible(False)

        # set text label
        for i, v in enumerate(y):
            plt.text(v + 1, i, str(v) + '%', color='black', va='center')
        
        # save fig
        png_file = os.path.join(temp_folder, cur_sample_name + '_mapping_figure.png')
        plt.savefig(png_file, dpi = 300)
        plt.close()

        cur_sample_png = Image(png_file, width=4*inch, height=0.2*inch)

        if cur_sample_subtype == 'Not RSV':
            png_file = os.path.join(file_path, 'Resource','error.png')
        else:
            if int(cur_sample_map) > 80:
                png_file = os.path.join(file_path, 'Resource','correct.png')
            else:
                png_file = os.path.join(file_path, 'Resource','warning.png')
        sign_png = Image(png_file, width=0.2*inch, height=0.2*inch)
        
        cur_sample_link = Paragraph(f"<a href='#{cur_sample_name}'>{cur_sample_name}</a>", custom_style)
        
        # genotype for whole genome
        if data[row][7] in ['A','B']:
            genotype_text = data[row][9] + '*'
        else:
            genotype_text = data[row][7]
        
        # genotype for G gene
        if data[row][6] == 'unassigned':
            g_genotype_text = data[row][8] + '*'
        else:
            g_genotype_text = data[row][6]
        if data[row][10] < 20:
            g_genotype_text = data[row][8]
        if g_genotype_text == 'Low G coverage':
            g_genotype_text = 'Low cov'

        data[row] = [cur_sample_link, data[row][1], cur_sample_png, genotype_text, g_genotype_text, sign_png]

    df_columns = ['Sample name', 'Pass QC', 'Mapping rate', 'Clade', 'G-Clade', 'Sign']
    data.insert(0, df_columns)
    col_widths = [1.5*inch, 0.7*inch, 4*inch, 0.8*inch, 0.8*inch, 0.4*inch]

    table = Table(data, colWidths = col_widths)
    #table._argW = [100,60,150,80]

    # Create a TableStyle and add it to the table
    style = TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), colors.grey),  # set background color for the first row to grey
        ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),  # set text color for the first row to white

        ('LINEBEFORE', (0, 0), (-1, -1), 0, (1, 1, 1)),  # Hide left border
        ('LINEAFTER', (0, 0), (-1, -1), 0, (1, 1, 1)),   # Hide right border

        ('ALIGN', (0, 0), (-1, -1), 'CENTER'),  # align all cells to the center
        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),  # set font for the first row to Helvetica-Bold
        ('FONTSIZE', (0, 0), (-1, -1), 9),  # set font size entire table

        ('BOTTOMPADDING', (0, 0), (-1, 0), 12),  # add more space below the text of the first row
        ('BACKGROUND', (0, 1), (-1, -1), colors.white),  # set background color for the other rows to beige
        ('GRID', (0,0), (-1,-1), 1, colors.black)  # set grid color to black and line width to 1
    ])

    # Set the background color of each cell based on its value
    for row in range(1, len(data)):
        for col in [1]:
            value = percentage_to_number(data[row][col])
            style.add('BACKGROUND', (col,row), (col,row), color_gradient(value))

    color_dict = {'A': '#ffcccc', 'B': '#ccccff', 'N': '#d3d3d3'}
    for row in range(1, len(data)):
        for col in [3]:
            value = data[row][col]
            value = value[0]
            style.add('BACKGROUND', (col,row), (col,row), color_dict[value])
            style.add('BACKGROUND', (col + 1,row), (col + 1,row), color_dict[value])

    table.setStyle(style)

    elements.append(table)

    spacer = Spacer(1, 12)  # width and height in points
    elements.append(spacer)

    # divide line
    #line = Line(line_width)  # width of line
    #elements.append(line)

    ######################################################################
    # section 1.5, coverage heatmap
    ######################################################################
    elements.append(PageBreak())
    subtitle_cov_section = Paragraph("Coverage by genes", styles['Heading2'])
    elements.append(subtitle_cov_section)

    df = pd.read_csv(csv_file, skiprows = 1)
    plot_df = df.iloc[:, [16,17,18,19,20,21,22,23,24,25,26]]
    plot_df.index = df.iloc[:,0]
    png_file = os.path.join(temp_folder, 'coverage_heatmap.png')
    make_coverage_heatmap(plot_df, png_file)

    pil_img = PILImage.open(png_file)
    original_width, original_height = pil_img.size
    new_width = 480
    new_height = original_height * (new_width / original_width)
    if new_height > 655:
        new_height = 650
        new_width = original_width * (new_height / original_height)

    img_heatmap = Image(png_file, width=new_width, height=new_height)
    img_heatmap.hAlign = 'LEFT'
    elements.append(img_heatmap)

    ######################################################################
    # section 2, phylogenetic tree
    ######################################################################
    tree_A_file = os.path.join(temp_folder, "RSV_A.png")
    tree_B_file = os.path.join(temp_folder, "RSV_B.png")
    if os.path.isfile(tree_A_file) or os.path.isfile(tree_B_file):
        elements.append(PageBreak())
        subtitle_tree_section = Paragraph("Phylogenetic Analysis", styles['Heading2'])
        elements.append(subtitle_tree_section)

        phylogenetic_text = "Phylogenetic trees are genertated using whole genome sequences. "
        paragraph = Paragraph(phylogenetic_text, styles['BodyText'])
        elements.append(paragraph)

        # subtype A tree
        if os.path.isfile(tree_A_file):
            subtitle_tree_section = Paragraph("Subtype A", styles['Heading3'])
            elements.append(subtitle_tree_section)

            pil_img = PILImage.open(tree_A_file)
            original_width, original_height = pil_img.size
            new_width = 480
            new_height = original_height * (new_width / original_width)
            img_tree = Image(tree_A_file, width=new_width, height=new_height)
            img_tree.hAlign = 'LEFT'
            elements.append(img_tree)

        # subtype B tree
        if os.path.isfile(tree_B_file):
            if os.path.isfile(tree_A_file):
                elements.append(PageBreak())

            subtitle_tree_section = Paragraph("Subtype B", styles['Heading3'])
            elements.append(subtitle_tree_section)
            
            pil_img = PILImage.open(tree_B_file)
            original_width, original_height = pil_img.size
            new_width = 480
            new_height = original_height * (new_width / original_width)
            img_tree = Image(tree_B_file, width=new_width, height=new_height)
            img_tree.hAlign = 'LEFT'
            elements.append(img_tree)

    ######################################################################
    # section 3, Details for each sample
    ######################################################################
    elements.append(PageBreak())
    subtitle_summary_section = Paragraph("Details", styles['Heading2'])
    elements.append(subtitle_summary_section)

    df = pd.read_csv(csv_file, skiprows = 1, index_col=0)
    Sample_folders = [f for f in os.listdir(mapres_folder) if os.path.isdir(os.path.join(mapres_folder, f))]
    sample_count = 0
    for cur_folder in sorted(Sample_folders):
        if os.path.isfile(cur_folder):
            continue
        if cur_folder == "Temp":
            continue
        #print(cur_folder + "\n")
        
        if sample_count != 0:
            elements.append(PageBreak())
        sample_count += 1
         
        subtitle_cur_sample = Paragraph(f"<a name='{cur_folder}'/>Sample: {cur_folder}", styles['Heading3'])
        #elements.append(ActionFlowable('SetPageDestination', name = cur_folder))
        elements.append(subtitle_cur_sample)

        # ######################################## genotype calling
        subtitle = Paragraph('Genotype calls', styles['Heading4'])
        elements.append(subtitle)

        # genotype for whole genome
        if df.loc[cur_folder][27] in ['A','B']:
            genotype_text = df.loc[cur_folder][29] + '*'
        else:
            genotype_text = df.loc[cur_folder][27]
        
        # genotype for G gene
        if df.loc[cur_folder][26] == 'unassigned':
            g_genotype_text = df.loc[cur_folder][28] + '*'
        else:
            g_genotype_text = df.loc[cur_folder][26]

        if df.loc[cur_folder][22] < 20:
            g_genotype_text = df.loc[cur_folder][28]

        F_protein_mutation_text = df.loc[cur_folder][14]
        if isinstance(F_protein_mutation_text, str):
            print(cur_folder)
            print(F_protein_mutation_text)
            F_protein_mutation_text = F_protein_mutation_text.replace("|", ", ")
            F_protein_mutation_text = F_protein_mutation_text.replace('(Reported)','')
        else:
            F_protein_mutation_text = ''
            
        fig_size = '20'
        if genotype_text == "Not RSV":
            genotype_para  = '<img src="' + os.path.join(file_path, 'Resource','error.png') + '" valign="middle" width="' + fig_size + '" height="' + fig_size + '"/>  ' + genotype_text
        else:
            cur_sample_map = int(df.loc[cur_folder][7])
            if cur_sample_map > 80:
                genotype_para  = '<img src="' + os.path.join(file_path, 'Resource','correct.png') + '" valign="middle" width="' + fig_size + '" height="' + fig_size + '"/>  <b>' + genotype_text + '</b> (based on whole genome)'
            else:
                genotype_para  = '<img src="' + os.path.join(file_path, 'Resource','warning.png') + '" valign="middle" width="' + fig_size + '" height="' + fig_size + '"/>  <b>' + genotype_text + '</b> (based on whole genome)'
            
            genotype_para += f";  <b>{g_genotype_text}</b> (based on G-ectodomain)<br/><br/>"
            #genotype_para += f"Genotype Resource:   <b><a href='https://nextstrain.org/rsv/a/genome'>Nextstrain (click for details), Data updated 2024-05-21</a></b>"
            genotype_para += f"F protein mutations:  <b>{F_protein_mutation_text}</b> <br/><br/>"

        paragraph = Paragraph(genotype_para)
        
        spacer = Spacer(1, 4)  # width and height in points
        elements.append(spacer)
        elements.append(paragraph)

        style = TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.grey),  # set background color for the first row to grey
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),  # set text color for the first row to white

            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),  # align all cells to the center
            ('LEFTMARGIN', (0,0), (-1,-1), 0),     # align table to the left of page
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),  # set font for the first row to Helvetica-Bold
            ('FONTSIZE', (0, 0), (-1, 0), 9),  # set font size for the first row to 10

            ('BOTTOMPADDING', (0, 0), (-1, 0), 12),  # add more space below the text of the first row
            ('BACKGROUND', (0, 1), (-1, -1), colors.beige),  # set background color for the other rows to beige
            ('GRID', (0,0), (-1,-1), 1, colors.black)  # set grid color to black and line width to 1
        ])

        # most identical strains on GISAID and Next strain
        if genotype_text == "Not RSV":
            pass
        else:
            # blast hit of GISAID
            blast_out_file = os.path.join(mapres_folder, cur_folder, 'Genotype', 'blastn_res_gisaid.tsv')
            gisaid_name, blast_alignment_length_gisaid, blast_pct_identity_gisaid = get_best_blast_hit(blast_out_file)
            gisaid_name = gisaid_name.split('|')
            #text = Paragraph(f"GISAID best hit: <b>{gisaid_name[0]}</b>, Identity: {blast_pct_identity_gisaid}%, Alignment Length: {blast_alignment_length_gisaid}")
            #elements.append(text)
            
            # blast hit of NextStrain
            genotype_file = os.path.join(mapres_folder, cur_folder, 'Genotype', 'Genotype.txt')
            subtype_str, blast_pct_identity_nextstrain, blast_alignment_length_nextstrain, nextstrain_name = get_genotype_res(genotype_file)
            #text = Paragraph(f"NextStrain best hit: <b>{nextstrain_name}</b>, Identity: {blast_pct_identity_nextstrain}%, Alignment Length: {blast_alignment_length_nextstrain}")
            #elements.append(text)

            data = [
                ['','Best hit','Identity(%)','Alignment Length'],
                ['GISAID', gisaid_name[0], blast_pct_identity_gisaid, blast_alignment_length_gisaid],
                ['NextStrain', nextstrain_name, blast_pct_identity_nextstrain, blast_alignment_length_nextstrain],
            ]

            table = Table(data)
            table.setStyle(style)
            elements.append(table)

        # tree figure for current strain
        if genotype_text == "Not RSV":
            pass
        else:
            # trees
            tree_png_file = os.path.join(mapres_folder, cur_folder, 'Genotype','genotype.png')
            if os.path.exists(tree_png_file):
                subtitle = Paragraph('Phylogenetic analysis: Whole genome', styles['Heading4'])
                elements.append(subtitle)
                
                pil_img = PILImage.open(tree_png_file)
                original_width, original_height = pil_img.size
                new_width = 430
                new_height = original_height * (new_width / original_width)
                img_tree = Image(tree_png_file, width=new_width, height=new_height)
                img_tree.hAlign = 'CENTER'
                elements.append(img_tree)
                elements.append(PageBreak())

            tree_png_file = os.path.join(mapres_folder, cur_folder, 'Genotype','G_gene_genotype.png')
            if os.path.exists(tree_png_file):
                subtitle = Paragraph('Phylogenetic analysis: G-ectodomain', styles['Heading4'])
                elements.append(subtitle)

                pil_img = PILImage.open(tree_png_file)
                original_width, original_height = pil_img.size
                new_width = 480
                new_height = original_height * (new_width / original_width)
                img_tree = Image(tree_png_file, width=new_width, height=new_height)
                img_tree.hAlign = 'CENTER'
                elements.append(img_tree)

        # ######################################## QC details
        subtitle = Paragraph('QC Details', styles['Heading4'])
        elements.append(subtitle)

        data = [
                ['','Total reads', 'Q20 pct', 'Q30 pct'],
                ['Raw data', df.loc[cur_folder][0], float_to_percentage(df.loc[cur_folder][1]), float_to_percentage(df.loc[cur_folder][2])], 
                ['Filtered data', df.loc[cur_folder][3], float_to_percentage(df.loc[cur_folder][4]), float_to_percentage(df.loc[cur_folder][5])]
            ]

        table = Table(data, hAlign='LEFT')
        table.setStyle(style)
        elements.append(table)

        # ######################################## coverage summary
        if genotype_text == "Not RSV":
            pass
        else:
            elements.append(PageBreak())

        ref_genotype_text = df.loc[cur_folder][13]
        if ref_genotype_text == "Not RSV":
            pass
        else:
            subtitle = Paragraph('Coverage summary', styles['Heading4'])
            elements.append(subtitle)

            reference_genome_path = os.path.join(mapres_folder, cur_folder, 'reference')
            subfolders = [ f.name for f in os.scandir(reference_genome_path) if f.is_dir() ]
            reference_genome_accession = subfolders[0]

            ref_text = f"Reference genome:   <b><a href='https://www.ncbi.nlm.nih.gov/nuccore/{reference_genome_accession}'>{reference_genome_accession}, click for details</a></b>"
            wig_file = os.path.join(mapres_folder, cur_folder, 'mapping', 'alignments.cov.wig')

            if genotype_text == 'SubtypeA':
                genome_intervals = [(70,420),(599,375),(1111,1176),(2318,726),(3226,771),(4266,195),(4652,966),(5697,1725),(7640,585),(8158,273),(8532,6498)] # subtype A
                genome_xticks = [70,800,1600,2700,3550,4366,5000,6300,7650,8400,11500]
                gene_names = ['NS1','NS2','N','P','M','SH','G','F','M2-1','M2-2','L']
            elif genotype_text == 'SubtypeB':
                genome_intervals = [(57,420),(584,375),(1097,1176),(2305,726),(3220,771),(4259,198),(4646,933),(5676,1725),(7627,588),(8180,273),(8518,6501)] # subtype B
                genome_xticks = [70,800,1600,2700,3550,4366,5000,6300,7650,8400,11500]
                gene_names = ['NS1','NS2','N','P','M','SH','G','F','M2-1','M2-2','L']
            else:
                genome_intervals = [(70,420),(599,375),(1111,1176),(2318,726),(3226,771),(4266,195),(4652,966),(5697,1725),(7640,585),(8158,273),(8532,6498)] # subtype A
                genome_xticks = [70,800,1600,2700,3550,4366,5000,6300,7650,8400,11500]
                gene_names = ['NS1','NS2','N','P','M','SH','G','F','M2-1','M2-2','L']
            
            paragraph = Paragraph(ref_text, styles['BodyText'])
            elements.append(paragraph)

            # read wig file
            cov_df = pd.read_csv(wig_file, skiprows=3, delimiter='\t', index_col=0)
            cov = cov_df.sum(axis=1)
            cov_pos = cov_df.index.tolist()

            # coverage details
            line = Line(300)  # width of line
            elements.append(line)
            cov_text = '<b>Genome size (BP):</b>   ' + str(max(cov_pos))
            paragraph = Paragraph(cov_text, styles['BodyText'])
            elements.append(paragraph)
            cov_text = '<b>20x coverage (BP):</b>   ' + str(count_greater_than_n(cov, 19)) + ' ('+ str(float_to_percentage(count_greater_than_n(cov, 19)/ max(cov_pos))) + ')'
            paragraph = Paragraph(cov_text, styles['BodyText'])
            elements.append(paragraph)
            cov_text = '<b>50x coverage (BP):</b>   ' + str(count_greater_than_n(cov, 49)) + ' ('+ str(float_to_percentage(count_greater_than_n(cov, 49)/ max(cov_pos))) + ')'
            paragraph = Paragraph(cov_text, styles['BodyText'])
            elements.append(paragraph)
            cov_text = '<b>100x coverage (BP):</b>   ' + str(count_greater_than_n(cov, 99)) + ' ('+ str(float_to_percentage(count_greater_than_n(cov, 99)/ max(cov_pos))) + ')'
            paragraph = Paragraph(cov_text, styles['BodyText'])
            elements.append(paragraph)
            cov_text = '<b>500x coverage (BP):</b>   ' + str(count_greater_than_n(cov, 499)) + ' ('+ str(float_to_percentage(count_greater_than_n(cov, 499)/ max(cov_pos))) + ')'
            paragraph = Paragraph(cov_text, styles['BodyText'])
            elements.append(paragraph)
            cov_text = '<b>Average sequencing depth:</b>   ' +  str(int(sum(cov)/len(cov)))
            paragraph = Paragraph(cov_text, styles['BodyText'])
            elements.append(paragraph)
            line = Line(300)  # width of line
            spacer = Spacer(1, 2)  # width and height in points
            elements.append(spacer)
            elements.append(line)
            spacer = Spacer(1, 4)  # width and height in points
            elements.append(spacer)

            # coverage figure
            genome_colors = ['red', 'yellow', 'green', 'cyan', 'blue', 'magenta', '#CCCC00', '#800080', '#D2B48C', '#8B4513', '#ADD8E6'] 
            fig = plt.figure(figsize=(12, 4))
            gs = gridspec.GridSpec(2, 1, height_ratios=[4, 0.5], hspace=0.05)
            font_axis = {'family': 'sans-serif', 'color': 'black', 'weight': 'bold', 'size': 12}

            gs1 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs[0])
            ax_cov = plt.subplot(gs1[0])
            ax_cov.set_ylabel('Coverage', fontdict=font_axis)

            #ax_cov.bar(cov_pos, cov, color = 'gray')
            ax_cov.fill_between(cov_pos, cov, color = 'gray')
            ax_cov.set_xticks([])  # Set the positions for the ticks
            ax_cov.set_yscale('log')    # Set Y scale to log

            ax_cov.axhline(y=50, color='red', linestyle='--')
            ax_cov.axhline(y=500, color='blue', linestyle='--')

            gs2 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[1], hspace=0.05)
            ax_genome = plt.subplot(gs2[0], sharex=ax_cov)
            for interval, color in zip(genome_intervals, genome_colors):
                ax_genome.broken_barh([interval], (10, 9), facecolors=color, edgecolors=color)
            ax_genome.set_xticks(genome_xticks)  # Set the positions for the ticks
            ax_genome.set_xticklabels(gene_names)  # Set the labels for the ticks
            # Remove the y-ticks
            ax_genome.set_yticks([])
            ax_genome.set_ylabel('Genes', rotation=0, fontdict=font_axis, labelpad=30)
            ax_genome.yaxis.set_label_coords(0,-0.1)

            # Remove the border
            ax_genome.spines['top'].set_visible(False)
            ax_genome.spines['right'].set_visible(False)
            ax_genome.spines['bottom'].set_visible(False)
            ax_genome.spines['left'].set_visible(False)

            png_file = os.path.join(temp_folder, cur_folder + '_coverage_figure.png')
            plt.savefig(png_file, dpi = 300)
            plt.close()

            img_cov = Image(png_file, width=600, height=200)
            cov_title_text = f"The coverage plot is under <b>Log scale</b>, the red and blue lines indicate coverage = 50 and 500 <br/>"
            paragraph = Paragraph(cov_title_text, styles['BodyText'])
            elements.append(paragraph)
            elements.append(img_cov)

        # ######################################## SNP section
        if "Not RSV" in (genotype_text, ref_genotype_text):
            pass
        else:
            cutoff = 0.2
            gff_file = find_gff_files_in_path(os.path.join(mapres_folder, cur_folder, 'reference'))
            gff_path = os.path.join(mapres_folder, cur_folder, 'reference',gff_file[0])
            table_data = SNP_calling(wig_file, cutoff, genotype_text, gff_path, temp_folder, cur_folder, igv_cutoff)
            # heading
            subtitle = Paragraph('SNP details', styles['Heading4'])
            elements.append(subtitle)
            # figure
            png_file = os.path.join(temp_folder, cur_folder + '_snp_figure.png')
            img_cov = Image(png_file, width=600, height=200)
            elements.append(img_cov)
            # table
            style = TableStyle([
                ('BACKGROUND', (0, 0), (-1, 0), colors.lightblue),  # Header background color
                ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),  # Header text color
                ('ALIGN', (0, 0), (-1, -1), 'CENTER'),  # Center align all cells
                ('LINEABOVE', (0, 0), (-1, 0), 2, colors.black),  # Top border
                ('LINEBELOW', (0, 0), (-1, 0), 1, colors.black),  # Border under header
                ('LINEBELOW', (0, -1), (-1, -1), 2, colors.black),  # Bottom border
            ])
            #print(table_data)
            table = Table(table_data)
            table.setStyle(style)
            elements.append(table)


    # build the doc
    doc.build(elements)

##########################################
# generate html report
##########################################

def generate_html_report(file_path, csv_file, working_folder, mapres_folder, igv_cutoff):
    # get version info
    current_script_path = os.path.dirname(os.path.abspath(__file__))
    version_file_path = os.path.join(current_script_path, 'version.txt')
    current_version = get_version(version_file_path)

    sidebar_div = '<div class="sidebar">\n'
    main_content_div = '<div class="main-content">\n'

    # load template
    template_file = os.path.join(file_path, 'template.html')
    with open(template_file, 'r') as tem_handle:
        html_report_str = tem_handle.read()

    # set the temp folder under working folder
    temp_folder = os.path.join(working_folder, 'Temp')
    if not os.path.exists(temp_folder):
        # Create the folder
        exit(f"Folder '{temp_folder}' does not exists!")

    # resource path and file path
    html_report = os.path.join(working_folder, "Report.html")

    ######################################################################
    # section 1, title and Summary
    ######################################################################
    # make sidebar_div
    sidebar_div += f"<ul><li><h2><a onclick=\"showSection('summary_section')\">Summary</a></h2></li></ul>\n"

    title = "Detection of RSV from clinical samples"
    
    Logo = os.path.join(file_path, 'Resource','CAB.png')

    table_info_text = f"Pipeline Verions: {current_version}<br/> Subtypes of each reference are highlighted in different colors: <font color='red'>Subtype A</font> and <font color='blue'>Subtype B</font>"
    table_info_text += "<br/> Genotype calling is based on <a href='https://nextstrain.org/rsv/a/genome'  target=\"_blank\"><b>Nextstrain</b> and <b>Nextclade3</b>, Data updated 2024-05-21</a><br/> A clade lable with * indicates the clade of best blast hit strain due to negative results from NextClade3 <br/>Click the sample name to jump to the detail section<br/> "

    main_content_div += '<div id="summary_section" class="content-section">'
    main_content_div += f"<h1>{title}</h1>\n"
    base64_string = image_to_base64(Logo)
    main_content_div += f"<img src='data:image/png;base64,{base64_string}' alt=\"Logo\" width=1000px>\n"
    main_content_div += f"<h2>Summary</h2>\n"
    main_content_div += f"<p>{table_info_text}</p>\n"

    # load table for summary section
    df = pd.read_csv(csv_file, skiprows=1)
    df.iloc[:,7] = df.iloc[:,7].apply(lambda x: '{:.2%}'.format(x/100)) # QC rate

    df = df.iloc[:, [0, 7, 8, 12,13,14, 27, 28, 29, 30, 22]] # name, QC rate, mapping rate, subtype, reference_accession, ref_subtype, Gtype, Wtype, Gtype_blast, Wtype_blast, G_cov
    data = df.values.tolist()

    for row in range(0, len(data)):
        cur_sample_name = data[row][0]
        cur_sample_map = data[row][2]
        cur_sample_subtype = data[row][3]
        ref_id = data[row][4]
        ref_type = data[row][5]

        cur_sample_png = os.path.join(temp_folder, cur_sample_name + '_mapping_figure.png')

        if cur_sample_subtype == 'Not RSV':
            sign_png = os.path.join(file_path, 'Resource','error.png')
        else:
            if int(cur_sample_map) > 80:
                sign_png = os.path.join(file_path, 'Resource','correct.png')
            else:
                sign_png = os.path.join(file_path, 'Resource','warning.png')
        
        cur_sample_link = f"<a onclick=\"showSection('{cur_sample_name}')\">{cur_sample_name}</a>"
        
        # genotype for whole genome
        if data[row][7] in ['A','B']:
            genotype_text = data[row][9] + '*'
        else:
            genotype_text = data[row][7]
        
        # genotype for G gene
        if data[row][6] == 'unassigned':
            g_genotype_text = data[row][8] + '*'
        else:
            g_genotype_text = data[row][6]
        if data[row][10] < 20:
            g_genotype_text = data[row][8]
        if g_genotype_text == 'Low G coverage':
            g_genotype_text = 'Low cov'

        mapping_fig_base64_string = image_to_base64(cur_sample_png)
        mapping_fig = f"<img src='data:image/png;base64,{mapping_fig_base64_string}' alt=\"mapping_fig\" style='margin-top:0px;height:30px'>\n"
        sign_fig_base64_string = image_to_base64(sign_png)
        sign_fig = f"<img src='data:image/png;base64,{sign_fig_base64_string}' alt=\"sign_fig\" style='margin-top:0px;height:30px'>\n"

        data[row] = [cur_sample_link, data[row][1], mapping_fig, genotype_text, g_genotype_text, sign_fig]

    df_columns = ['Sample name', 'Pass QC', 'Mapping rate', 'Clade', 'G-Clade', 'Sign']
    
    # Set the background color of each cell based on its value
    bg_color = generate_empty_2d_array(data)
    color_dict = {'A': '#ffcccc', 'B': '#ccccff', 'N': '#d3d3d3'}
    for row in range(0, len(data)):
        col = 1
        value = int(percentage_to_number(data[row][col]))
        bg_color[row][col] = int_to_gradient_color(value)

        col = 3
        value = data[row][col]
        value = value[0]
        bg_color[row][col] = color_dict[value]
        bg_color[row][col+1] = color_dict[value]


    html_table = array_to_html_table(data, df_columns, color = bg_color, table_id = 'summary_table', table_class = 'summary_table_class')
    main_content_div += f"{html_table}\n"

    # coverage heatmap
    main_content_div += f"<h2>Coverage by genes</h2>\n"

    png_file = os.path.join(temp_folder, 'coverage_heatmap.png')
    base64_string = image_to_base64(png_file)
    main_content_div += f"<img src='data:image/png;base64,{base64_string}' alt=\"coverage_heatmap\" style='margin-top:0px;width:1500px'>\n"

    # phylogenetic tree
    tree_A_file = os.path.join(temp_folder, "RSV_A.png")
    tree_B_file = os.path.join(temp_folder, "RSV_B.png")
    if os.path.isfile(tree_A_file) or os.path.isfile(tree_B_file):
        main_content_div += f"<h2>Phylogenetic Analysis</h2>\n"
        main_content_div += f"<p>Phylogenetic trees are genertated using whole genome sequences.</p>\n"

        # subtype A tree
        if os.path.isfile(tree_A_file):
            main_content_div += f"<h3>Subtype A</h3>\n"
            base64_string = image_to_base64(tree_A_file)
            main_content_div += f"<img src='data:image/png;base64,{base64_string}' alt=\"tree_A_file\" style='margin-top:0px;width:1500px'>\n"

        # subtype B tree
        if os.path.isfile(tree_B_file):
            main_content_div += f"<h3>Subtype B</h3>\n"
            base64_string = image_to_base64(tree_B_file)
            main_content_div += f"<img src='data:image/png;base64,{base64_string}' alt=\"tree_B_file\" style='margin-top:0px;width:1500px'>\n"

    main_content_div += '</div>\n'

    ######################################################################
    # section 3, Details for each sample
    ######################################################################
    sidebar_div += '<h2>Individual samples:</h2>\n'
    sidebar_div += '<ul>\n'

    df = pd.read_csv(csv_file, skiprows = 1, index_col=0)
    Sample_folders = [f for f in os.listdir(mapres_folder) if os.path.isdir(os.path.join(mapres_folder, f))]
    for cur_folder in sorted(Sample_folders):
        if os.path.isfile(cur_folder):
            continue
        if cur_folder == "Temp":
            continue
        
        section_id = cur_folder

        sidebar_div += f"<li><a onclick=\"showSection('{section_id}')\">{section_id}</a></li>"
        main_content_div += f"<div id=\"{section_id}\" class=\"content-section\">"

        main_content_div += f"<h2>Sample: {section_id}</h2>\n"

        # ######################################## genotype calling
        main_content_div += f"<h3>Genotype calls</h3>\n"

        # genotype for whole genome
        if df.loc[cur_folder][27] in ['A','B']:
            genotype_text = df.loc[cur_folder][29] + '*'
        else:
            genotype_text = df.loc[cur_folder][27]
        
        # genotype for G gene
        if df.loc[cur_folder][26] == 'unassigned':
            g_genotype_text = df.loc[cur_folder][28] + '*'
        else:
            g_genotype_text = df.loc[cur_folder][26]

        if df.loc[cur_folder][22] < 20:
            g_genotype_text = df.loc[cur_folder][28]

        F_protein_mutation_text = df.loc[cur_folder][14]
        if isinstance(F_protein_mutation_text, str):
            print(cur_folder)
            print(F_protein_mutation_text)
            F_protein_mutation_text = F_protein_mutation_text.replace("|", ", ")
            F_protein_mutation_text = F_protein_mutation_text.replace('(Reported)','')
        else:
            F_protein_mutation_text = ''
            
        if genotype_text == "Not RSV":
            base64_string = image_to_base64(os.path.join(file_path, 'Resource','error.png'))
            genotype_para  = f"<img src='data:image/png;base64,{base64_string}' style='margin-top:0px;width:30px'><b>{genotype_text}</b>"
        else:
            cur_sample_map = int(df.loc[cur_folder][7])
            if cur_sample_map > 80:
                base64_string = image_to_base64(os.path.join(file_path, 'Resource','correct.png'))
            else:
                base64_string = image_to_base64(os.path.join(file_path, 'Resource','warning.png'))
            genotype_para  = f"<img src='data:image/png;base64,{base64_string}' style='margin-top:0px;width:30px'><b>{genotype_text}</b> (based on whole genome)"
            

            genotype_para += f";  <b>{g_genotype_text}</b> (based on G-ectodomain)<br/><br/>"
            #genotype_para += f"Genotype Resource:   <b><a href='https://nextstrain.org/rsv/a/genome'>Nextstrain (click for details), Data updated 2024-05-21</a></b>"
            genotype_para += f"F protein mutations:  <b>{F_protein_mutation_text}</b> <br/><br/>"

        main_content_div += genotype_para


        # most identical strains on GISAID and Next strain
        if genotype_text == "Not RSV":
            pass
        else:
            # blast hit of GISAID
            blast_out_file = os.path.join(mapres_folder, cur_folder, 'Genotype', 'blastn_res_gisaid.tsv')
            gisaid_name, blast_alignment_length_gisaid, blast_pct_identity_gisaid = get_best_blast_hit(blast_out_file)
            gisaid_name = gisaid_name.split('|')
            
            # blast hit of NextStrain
            genotype_file = os.path.join(mapres_folder, cur_folder, 'Genotype', 'Genotype.txt')
            subtype_str, blast_pct_identity_nextstrain, blast_alignment_length_nextstrain, nextstrain_name = get_genotype_res(genotype_file)

            data = [
                ['GISAID', gisaid_name[0], blast_pct_identity_gisaid, blast_alignment_length_gisaid],
                ['NextStrain', nextstrain_name, blast_pct_identity_nextstrain, blast_alignment_length_nextstrain],
            ]

            table_id = f"blast_table_{section_id}"
            html_table = array_to_html_table(data, ['Database','Best hit','Identity(%)','Alignment Length'], color = None, table_id = table_id, table_class = 'blast_table_class')
            main_content_div += f"{html_table}\n"

        # tree figure for current strain
        if genotype_text == "Not RSV":
            pass
        else:
            # trees
            tree_png_file = os.path.join(mapres_folder, cur_folder, 'Genotype','genotype.png')
            if os.path.exists(tree_png_file):
                main_content_div += f"<h3>Phylogenetic analysis: Whole genome</h3>\n"
                base64_string = image_to_base64(tree_png_file)
                main_content_div += f"<img src='data:image/png;base64,{base64_string}' alt=\"{section_id}\" style='margin-top:0px;width:1200px'>\n"

            tree_png_file = os.path.join(mapres_folder, cur_folder, 'Genotype','G_gene_genotype.png')
            if os.path.exists(tree_png_file):
                main_content_div += f"<h3>Phylogenetic analysis: G-ectodomain</h3>\n"
                base64_string = image_to_base64(tree_png_file)
                main_content_div += f"<img src='data:image/png;base64,{base64_string}' alt=\"{section_id}\" style='margin-top:0px;width:1200px'>\n"

        # ######################################## QC details
        main_content_div += f"<h2>QC Details</h2>\n"

        data = [
                ['Raw data', df.loc[cur_folder][0], float_to_percentage(df.loc[cur_folder][1]), float_to_percentage(df.loc[cur_folder][2])], 
                ['Filtered data', df.loc[cur_folder][3], float_to_percentage(df.loc[cur_folder][4]), float_to_percentage(df.loc[cur_folder][5])]
            ]

        table_id = f"qc_table_{section_id}"
        html_table = array_to_html_table(data, ['Data type','Total reads', 'Q20 pct(%)', 'Q30 pct(%)'], color = None, table_id = table_id, table_class = 'qc_table_class')
        main_content_div += f"{html_table}\n"

        # ######################################## coverage summary
        ref_genotype_text = df.loc[cur_folder][13]
        if ref_genotype_text == "Not RSV":
            pass
        else:
            main_content_div += f"<h2>Coverage summary</h2>\n"

            reference_genome_path = os.path.join(mapres_folder, cur_folder, 'reference')
            subfolders = [ f.name for f in os.scandir(reference_genome_path) if f.is_dir() ]
            reference_genome_accession = subfolders[0]

            ref_text = f"Reference genome:   <b><a href='https://www.ncbi.nlm.nih.gov/nuccore/{reference_genome_accession}'>{reference_genome_accession}, click for details</a></b>"
            main_content_div += f"<p>{ref_text}</p>\n"
            
            # read wig file
            wig_file = os.path.join(mapres_folder, cur_folder, 'mapping', 'alignments.cov.wig')
            cov_df = pd.read_csv(wig_file, skiprows=3, delimiter='\t', index_col=0)
            cov = cov_df.sum(axis=1)
            cov_pos = cov_df.index.tolist()

            # coverage details
            cov_text = '<b>Genome size (BP):</b>   ' + str(max(cov_pos))
            main_content_div += f"<p>{cov_text}</p>\n"
            cov_text = '<b>20x coverage (BP):</b>   ' + str(count_greater_than_n(cov, 19)) + ' ('+ str(float_to_percentage(count_greater_than_n(cov, 19)/ max(cov_pos))) + ')'
            main_content_div += f"<p>{cov_text}</p>\n"
            cov_text = '<b>50x coverage (BP):</b>   ' + str(count_greater_than_n(cov, 49)) + ' ('+ str(float_to_percentage(count_greater_than_n(cov, 49)/ max(cov_pos))) + ')'
            main_content_div += f"<p>{cov_text}</p>\n"
            cov_text = '<b>100x coverage (BP):</b>   ' + str(count_greater_than_n(cov, 99)) + ' ('+ str(float_to_percentage(count_greater_than_n(cov, 99)/ max(cov_pos))) + ')'
            main_content_div += f"<p>{cov_text}</p>\n"
            cov_text = '<b>500x coverage (BP):</b>   ' + str(count_greater_than_n(cov, 499)) + ' ('+ str(float_to_percentage(count_greater_than_n(cov, 499)/ max(cov_pos))) + ')'
            main_content_div += f"<p>{cov_text}</p>\n"
            cov_text = '<b>Average sequencing depth:</b>   ' +  str(int(sum(cov)/len(cov)))
            main_content_div += f"<p>{cov_text}</p>\n"
           
            cov_title_text = f"The coverage plot is under <b>Log scale</b>, the red and blue lines indicate coverage = 50 and 500 <br/>"
            main_content_div += f"<p>{cov_title_text}</p>\n"
            png_file = os.path.join(temp_folder, cur_folder + '_coverage_figure.png')
            base64_string = image_to_base64(png_file)
            main_content_div += f"<img src='data:image/png;base64,{base64_string}' alt=\"{section_id}\" style='margin-top:0px;width:1500px'>\n"

        # ######################################## SNP section
        if "Not RSV" in (genotype_text, ref_genotype_text):
            pass
        else:
            main_content_div += f"<h2>SNP details</h2>\n"
            # figure
            png_file = os.path.join(temp_folder, cur_folder + '_snp_figure.png')
            base64_string = image_to_base64(png_file)
            main_content_div += f"<img src='data:image/png;base64,{base64_string}' alt=\"{section_id}\" style='margin-top:0px;width:1500px'>\n"

            # table
            cutoff = 0.2
            gff_file = find_gff_files_in_path(os.path.join(mapres_folder, cur_folder, 'reference'))
            gff_path = os.path.join(mapres_folder, cur_folder, 'reference',gff_file[0])
            table_data = SNP_calling(wig_file, cutoff, genotype_text, gff_path, temp_folder, cur_folder, igv_cutoff)
            
            header = table_data.pop(0)
            table_id = f"snp_table_{section_id}"
            html_table = array_to_html_table(table_data, header, color = None, table_id = table_id, table_class = 'snp_table_class')
            main_content_div += f"{html_table}\n"


        main_content_div += '</div>\n'  # close div for current sample

    # close div, finish the html file

    sidebar_div += '</ul></div>'    # close ul and div for sidebar
    main_content_div += '</div>'    # close div for main content
    main_content_div += '</body></html>'    # end of body and html

    # generate html report
    html_report = os.path.join(working_folder, "Report.html")
    with open(html_report, 'w') as report:
        report.write(html_report_str)
        report.write(sidebar_div)
        report.write(main_content_div)

#csv_file = '../RSV_run/3083501/Report/Report.csv'
#working_folder = '../RSV_run/3083501/Report'
#generate_pdf_report(csv_file, working_folder)