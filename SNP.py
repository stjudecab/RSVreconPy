import os
import re
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
import matplotlib.colors as mcolors
from RSV_functions import get_sub_folders, elements_not_in_array, pct_sum, determine_subtype, processIGV, find_gene_at_position, parse_gff, last_pos_gff

def float_to_percentage(value):
    return "{:.2%}".format(value)

def SNP_calling(wig_file, cutoff, genotype_text, gff_path, out_path, prefix_name, cov_cutoff):
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

    cov_df = pd.read_csv(wig_file, skiprows=3, delimiter='\t', index_col=0, header=None, usecols=range(6))
    column_names = ['A', 'C', 'G','T','N']
    cov_df.columns = column_names
    cov = cov_df.sum(axis=1)
    percentage_df = cov_df.divide(cov, axis=0)
    percentage_df['cov'] = cov

    # remove columns that are lowly covered
    cov_df = cov_df[percentage_df['cov'] > cov_cutoff]
    percentage_df = percentage_df[percentage_df['cov'] > cov_cutoff]
    del percentage_df['cov']

    num_columns_to_check = 2  # Number of columns larger than the threshold required
    # Filter rows based on the condition
    filtered_pct_df = percentage_df[(percentage_df > cutoff).sum(axis=1) >= num_columns_to_check]
    filtered_df = cov_df[(percentage_df > cutoff).sum(axis=1) >= num_columns_to_check]
    
    # get a cutoff from gtf file
    end_pos_cds = last_pos_gff(gff_path)
    filtered_pct_df = filtered_pct_df[filtered_pct_df.index < end_pos_cds]
    filtered_df = filtered_df[filtered_df.index < end_pos_cds]
    cov = filtered_df.sum(axis=1)
    filtered_df['cov'] = cov

    genome_colors = ['red', 'yellow', 'green', 'cyan', 'blue', 'magenta', '#CCCC00', '#800080', '#D2B48C', '#8B4513', '#ADD8E6'] 
    fig = plt.figure(figsize=(12, 4))
    gs = gridspec.GridSpec(2, 1, height_ratios=[4, 0.5], hspace=0.05)
    font_axis = {'family': 'sans-serif', 'color': 'black', 'weight': 'bold', 'size': 12}

    gs1 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs[0])
    ax_cov = plt.subplot(gs1[0])
    ax_cov.set_ylabel('Coverage', fontdict=font_axis)

    # plot 1
    ax_cov.bar(filtered_pct_df.index, filtered_pct_df['A'], label='A', width = 100)
    ax_cov.bar(filtered_pct_df.index, filtered_pct_df['C'], bottom=filtered_pct_df['A'], label='C', width = 100)
    ax_cov.bar(filtered_pct_df.index, filtered_pct_df['G'], bottom=filtered_pct_df['A'] + filtered_pct_df['C'], label='G', width = 100)
    ax_cov.bar(filtered_pct_df.index, filtered_pct_df['T'], bottom=filtered_pct_df['A'] + filtered_pct_df['C'] + filtered_pct_df['G'], label='T', width = 100)
    ax_cov.legend()

    ax_cov.set_xticks([])  # Set the positions for the ticks

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

    png_file = os.path.join(out_path, prefix_name + '_snp_figure.png')
    plt.savefig(png_file, dpi = 300)
    plt.close()


    # table
    gene_positions = parse_gff(gff_path)
    data = [['Position','Coverage','A', 'A%', 'C', 'C%', 'G', 'G%', 'T', '%T', 'Gene']]

    for index_label in filtered_df.index:
        sub_list = [int(index_label), int(filtered_df.loc[index_label, 'cov']),
            int(filtered_df.loc[index_label, 'A']),float_to_percentage(filtered_pct_df.loc[index_label, 'A']),
            int(filtered_df.loc[index_label, 'C']),float_to_percentage(filtered_pct_df.loc[index_label, 'C']),
            int(filtered_df.loc[index_label, 'G']),float_to_percentage(filtered_pct_df.loc[index_label, 'G']),
            int(filtered_df.loc[index_label, 'T']),float_to_percentage(filtered_pct_df.loc[index_label, 'T']),
            find_gene_at_position(gene_positions, index_label)
        ]

        data.append(sub_list)

    return(data)
    

#wig_file = '/research/groups/cab/projects/automapper/common/lli75/RSV_run/3083501/Report/1_S1/mapping/alignments.cov.wig'
#cutoff = 0.2
#gff_path = '/research/groups/cab/projects/automapper/common/lli75/RSV_run/3083501/Report/1_S1/reference/MG813995.gff'
#out_path = './'
#prefix_name = 'test'
#SNP_calling(wig_file, cutoff, 'SubtypeA', gff_path, out_path, prefix_name)