import os
import glob
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
from RSV_functions import parse_gff, find_gene_at_position, find_gff_files_in_path, get_genotype_res
from SNP import SNP_calling
import seaborn as sns

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

# generate pdf report
def generate_pdf_report(csv_file, working_folder, mapres_folder, igv_cutoff):
    
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

    table_info_text = "Subtypes of each reference are highlighted in different colors: <font color='red'>Subtype A</font> and <font color='blue'>Subtype B</font>"
    table_info_text += "<br/> Genotype calling is based on <a href='https://nextstrain.org/rsv/a/genome'>Nextstrain (click for details), Data updated 2024-05-21</a><br/> Click the sample name to jump to the detail section<br/> "
    paragraph = Paragraph(table_info_text, styles['BodyText'])
    elements.append(paragraph)
    spacer = Spacer(1, 6)  # width and height in points
    elements.append(spacer)

    # load table for summary section
    df = pd.read_csv(csv_file, header=0)
    # Format the column as percentages
    df.iloc[:,7] = df.iloc[:,7].apply(lambda x: '{:.2%}'.format(x/100)) # QC rate
    #df.iloc[:,8] = df.iloc[:,8].apply(lambda x: '{:.2%}'.format(x/100))
    #df.iloc[:,12] = df.iloc[:,12].apply(lambda x: '{:.2%}'.format(x/100))

    df = df.iloc[:, [0, 7, 8, 12,13,14]] # name, QC rate, mapping rate, subtype
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
        plt.figure(figsize=(7, 0.4))
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

        cur_sample_png = Image(png_file, width=4.5*inch, height=0.25*inch)

        if cur_sample_subtype == 'Not RSV':
            png_file = os.path.join(file_path, 'Resource','error.png')
        else:
            if int(cur_sample_map) > 80:
                png_file = os.path.join(file_path, 'Resource','correct.png')
            else:
                png_file = os.path.join(file_path, 'Resource','warning.png')
        sing_png = Image(png_file, width=0.2*inch, height=0.2*inch)
        
        cur_sample_link = Paragraph(f"<a href='#{cur_sample_name}'>{cur_sample_name}</a>", custom_style)
        data[row] = [cur_sample_link, data[row][1], cur_sample_png, data[row][3], sing_png]

    df_columns = ['Sample name', 'Pass QC', 'Mapping rate', 'Subtype', 'Sign']
    data.insert(0, df_columns)
    col_widths = [1.5*inch, 0.7*inch, 4.7*inch, 0.9*inch, 0.4*inch]

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

    df = pd.read_csv(csv_file, header=0)
    plot_df = df.iloc[:, [16,17,18,19,20,21,22,23,24,25,26]]
    plot_df.index = df.iloc[:,0]
    png_file = os.path.join(temp_folder, 'coverage_heatmap.png')
    make_coverage_heatmap(plot_df, png_file)

    pil_img = PILImage.open(png_file)
    original_width, original_height = pil_img.size
    new_width = 480
    new_height = original_height * (new_width / original_width)
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

    df = pd.read_csv(csv_file, skiprows = 0, header=0, index_col=0)
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

        #print(cur_folder)
        genotype_text = df.loc[cur_folder][11]
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
            genotype_para  = '<img src="' + os.path.join(file_path, 'Resource','correct.png') + '" valign="middle" width="' + fig_size + '" height="' + fig_size + '"/>  ' + genotype_text + '(based on whole genome) <br/> <br/>'
            #genotype_para += f"Genotype Resource:   <b><a href='https://nextstrain.org/rsv/a/genome'>Nextstrain (click for details), Data updated 2024-05-21</a></b>"
            genotype_para += f"<b>F protein mutations:  {F_protein_mutation_text}</b> <br/>"

        paragraph = Paragraph(genotype_para)
        
        spacer = Spacer(1, 4)  # width and height in points
        elements.append(spacer)
        elements.append(paragraph)

        # most identical strains on GISAID and Next strain
        if genotype_text == "Not RSV":
            pass
        else:
            # blast hit of GISAID
            blast_out_file = os.path.join(mapres_folder, cur_folder, 'Genotype', 'blastn_res_gisaid.tsv')
            gisaid_name, blast_alignment_length, blast_pct_identity = get_best_blast_hit(blast_out_file)
            gisaid_name = gisaid_name.split('|')
            text = Paragraph(f"GISAID best hit: <b>{gisaid_name[0]}</b>, Identity: {blast_pct_identity}%, Alignment Length: {blast_alignment_length}")
            elements.append(text)
            
            # blast hit of NextStrain
            genotype_file = os.path.join(mapres_folder, cur_folder, 'Genotype', 'Genotype.txt')
            subtype_str, blast_pct_identity, blast_alignment_length, nextstrain_name = get_genotype_res(genotype_file)
            text = Paragraph(f"NextStrain best hit: <b>{nextstrain_name}</b>, Identity: {blast_pct_identity}%, Alignment Length: {blast_alignment_length}")
            elements.append(text)

        # tree figure for current strain
        if genotype_text == "Not RSV":
            pass
        else:
            # trees
            tree_png_file = os.path.join(mapres_folder, cur_folder, 'Genotype','genotype.png')
            pil_img = PILImage.open(tree_png_file)
            original_width, original_height = pil_img.size
            new_width = 400
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

#csv_file = '../RSV_run/3083501/Report/Report.csv'
#working_folder = '../RSV_run/3083501/Report'
#generate_pdf_report(csv_file, working_folder)