import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from PIL import Image as PILImage
from reportlab.lib.styles import getSampleStyleSheet, TA_LEFT
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Flowable, Table, Image, TableStyle
from reportlab.lib.colors import black
from reportlab.lib import colors
from reportlab.lib.units import inch

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

# Create a color gradient function
def color_gradient(value):
    """Return a color from green to red based on the input value (0 to 1)."""
    return colors.Color(1 - value, value, 0)

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

# generate pdf report
def generate_pdf_report(csv_file, working_folder):
    
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
    Logo = os.path.join(file_path, 'Resource','CAB_Hz.png')
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

    # section 1, summary
    subtitle_summary_section = Paragraph("Summary", styles['Heading2'])
    elements.append(subtitle_summary_section)


    spacer = Spacer(1, 12)  # width and height in points
    elements.append(spacer)

    # load table for summary section
    df = pd.read_csv(csv_file, header=0)
    # Format the column as percentages
    df.iloc[:,7] = df.iloc[:,7].apply(lambda x: '{:.2%}'.format(x/100))
    #df.iloc[:,8] = df.iloc[:,8].apply(lambda x: '{:.2%}'.format(x/100))
    #df.iloc[:,12] = df.iloc[:,12].apply(lambda x: '{:.2%}'.format(x/100))

    df = df.iloc[:, [0, 7, 8, 12,16]]
    data = df.values.tolist()
    
    # make figure for mapping rate for each sample
    for row in range(0, len(data)):
        cur_sample_name = data[row][0]
        cur_sample_mapA = data[row][2]
        cur_sample_mapB = data[row][3]

        # data
        x = ["Subtype A","Subtype B"]
        y = [cur_sample_mapA, cur_sample_mapB]
        mapping_colors = ['blue', 'green']

        # plot
        plt.figure(figsize=(6.5, 0.5))
        plt.barh(x, y, color=mapping_colors)
        plt.xlim(0, 100)

        # hide tick, label and borders
        plt.xlabel('')
        plt.ylabel('')
        plt.title('')
        plt.xticks([])
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

        cur_sample_png = Image(png_file, width=4.5*inch, height=0.4*inch)
        data[row] = [cur_sample_name, data[row][1], cur_sample_png, data[row][4]]

    df_columns = ['Sample name', 'Pass QC', 'Mapping rate', 'Subtype']
    data.insert(0, df_columns)
    table = Table(data)
    #table._argW = [100,60,150,80]

    # Create a TableStyle and add it to the table
    style = TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), colors.grey),  # set background color for the first row to grey
        ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),  # set text color for the first row to white

        ('LINEBEFORE', (0, 0), (-1, -1), 0, (1, 1, 1)),  # Hide left border
        ('LINEAFTER', (0, 0), (-1, -1), 0, (1, 1, 1)),   # Hide right border

        ('ALIGN', (0, 0), (-1, -1), 'CENTER'),  # align all cells to the center
        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),  # set font for the first row to Helvetica-Bold
        ('FONTSIZE', (0, 0), (-1, 0), 9),  # set font size for the first row to 8

        ('BOTTOMPADDING', (0, 0), (-1, 0), 12),  # add more space below the text of the first row
        ('BACKGROUND', (0, 1), (-1, -1), colors.white),  # set background color for the other rows to beige
        ('GRID', (0,0), (-1,-1), 1, colors.black)  # set grid color to black and line width to 1
    ])

    # Set the background color of each cell based on its value
    for row in range(1, len(data)):
        for col in [1]:
            value = percentage_to_number(data[row][col])
            style.add('BACKGROUND', (col,row), (col,row), color_gradient(value))

    #for row in range(1, len(data)):
    #    for col in [2,3]:
    #        value = percentage_to_number(data[row][col])
    #        style.add('BACKGROUND', (col,row), (col,row), color_gradient_blue(value))

    table.setStyle(style)

    elements.append(table)

    spacer = Spacer(1, 12)  # width and height in points
    elements.append(spacer)

    # divide line
    line = Line(line_width)  # width of line
    elements.append(line)

    # section 2, Details
    subtitle_summary_section = Paragraph("Details", styles['Heading2'])
    elements.append(subtitle_summary_section)

    df = pd.read_csv(csv_file, skiprows = 0, header=0, index_col=0)
    Sample_folders = [f for f in os.listdir(working_folder) if os.path.isdir(os.path.join(working_folder, f))]
    for cur_folder in Sample_folders:
        if os.path.isfile(cur_folder):
            continue
        if cur_folder == "Temp":
            continue
        #print(cur_folder + "\n")

        subtitle_cur_sample = Paragraph("Sample: " + cur_folder, styles['Heading3'])
        elements.append(subtitle_cur_sample)

        # genotype calling
        subtitle = Paragraph('Genotype calls', styles['Heading4'])
        elements.append(subtitle)

        genotype_text = df.loc[cur_folder][15]
        paragraph = Paragraph(genotype_text, styles['BodyText'])
        elements.append(paragraph)

        # QC details
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


        # coverage summary
        subtitle = Paragraph('Coverage summary', styles['Heading4'])
        elements.append(subtitle)

        if genotype_text == 'SubTypeA':
            ref_text = 'Reference genome:   <b>EPI_ISL_412866_hRSV/A/England/397/2017</b>'
            wig_file = os.path.join(working_folder, cur_folder, 'SubTypeA', 'alignments.cov.wig')
            genome_intervals = [(70,420),(599,375),(1111,1176),(2318,726),(3226,771),(4266,195),(4652,966),(5697,1725),(7640,585),(8158,273),(8532,6498)] # subtype A
            genome_xticks = [70,800,1600,2700,3550,4366,5000,6300,7650,8400,11500]
            gene_names = ['NS1','NS2','N','P','M','SH','G','F','M2-1','M2-2','L']
        elif genotype_text == 'SubTypeB':
            ref_text = 'Reference genome:   <b>EPI_ISL_1653999_hRSV/B/Australia/VIC-RCH056/2019</b>'
            wig_file = os.path.join(working_folder, cur_folder, 'SubTypeB', 'alignments.cov.wig')
            genome_intervals = [(57,420),(584,375),(1097,1176),(2305,726),(3220,771),(4259,198),(4646,933),(5676,1725),(7627,588),(8180,273),(8518,6501)] # subtype B
            genome_xticks = [70,800,1600,2700,3550,4366,5000,6300,7650,8400,11500]
            gene_names = ['NS1','NS2','N','P','M','SH','G','F','M2-1','M2-2','L']
        else:
            ref_text = 'Reference genome:   <b>EPI_ISL_412866_hRSV/A/England/397/2017</b>'
            wig_file = os.path.join(working_folder, cur_folder, 'SubTypeA', 'alignments.cov.wig')
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
        cov_text = '<b>Genome size (BP)p:</b>   ' + str(max(cov_pos))
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
        elements.append(img_cov)

    # build the doc
    doc.build(elements)

#csv_file = '../RSV_run/3083501/Report/Report.csv'
#working_folder = '../RSV_run/3083501/Report'
#generate_pdf_report(csv_file, working_folder)
