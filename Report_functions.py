import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image as PILImage
from reportlab.lib.styles import getSampleStyleSheet, TA_LEFT
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Flowable, Table, Image, TableStyle
from reportlab.lib.colors import black
from reportlab.lib import colors


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
        return "Input is not a numpy.float64 value"

# generate pdf report
def generate_pdf_report(csv_file, working_folder):
    
    line_width = 440

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
    new_height = 80
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

    #text = "qqq"
    #paragraph = Paragraph(text, styles['BodyText'])
    #elements.append(paragraph)

    spacer = Spacer(1, 12)  # width and height in points
    elements.append(spacer)

    # load table for summary section
    df = pd.read_csv(csv_file, skiprows=2)
    # Format the column as percentages
    df.iloc[:,7] = df.iloc[:,7].apply(lambda x: '{:.2%}'.format(x/100))
    df.iloc[:,8] = df.iloc[:,8].apply(lambda x: '{:.2%}'.format(x/100))
    df.iloc[:,12] = df.iloc[:,12].apply(lambda x: '{:.2%}'.format(x/100))

    df = df.iloc[:, [0, 7, 8, 12,16]]
    df.columns = ['Sample name', 'Pass QC rate', 'Mapping SubTypeA', 'Mapping SubTypeB', 'Subtype']
    data = df.values.tolist()
    data.insert(0, df.columns.tolist())
    table = Table(data)
    table._argW = [100,100,100,100,100]

    # Create a TableStyle and add it to the table
    style = TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), colors.grey),  # set background color for the first row to grey
        ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),  # set text color for the first row to white

        ('ALIGN', (0, 0), (-1, -1), 'CENTER'),  # align all cells to the center
        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),  # set font for the first row to Helvetica-Bold
        ('FONTSIZE', (0, 0), (-1, 0), 9),  # set font size for the first row to 8

        ('BOTTOMPADDING', (0, 0), (-1, 0), 12),  # add more space below the text of the first row
        ('BACKGROUND', (0, 1), (-1, -1), colors.beige),  # set background color for the other rows to beige
        ('GRID', (0,0), (-1,-1), 1, colors.black)  # set grid color to black and line width to 1
    ])

    # Set the background color of each cell based on its value
    for row in range(1, len(data)):
        for col in [1]:
            value = percentage_to_number(data[row][col])
            style.add('BACKGROUND', (col,row), (col,row), color_gradient(value))

    for row in range(1, len(data)):
        for col in [2,3]:
            value = percentage_to_number(data[row][col])
            style.add('BACKGROUND', (col,row), (col,row), color_gradient_blue(value))

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

    df = pd.read_csv(csv_file, skiprows = 1, header=0, index_col=0)
    Sample_folders = [f for f in os.listdir(working_folder) if os.path.isdir(os.path.join(working_folder, f))]
    for cur_folder in Sample_folders:
        if os.path.isfile(cur_folder):
            continue

        subtitle_cur_sample = Paragraph(cur_folder, styles['Heading3'])
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
        elif genotype_text == 'SubTypeB':
            ref_text = 'Reference genome:   <b>EPI_ISL_1653999_hRSV/B/Australia/VIC-RCH056/2019</b>'
            wig_file = os.path.join(working_folder, cur_folder, 'SubTypeB', 'alignments.cov.wig')
        else:
            ref_text = 'Reference genome:   <b>EPI_ISL_412866_hRSV/A/England/397/2017</b> and <b>EPI_ISL_1653999_hRSV/B/Australia/VIC-RCH056/2019</b>'
            wig_file = os.path.join(working_folder, cur_folder, 'SubTypeA', 'alignments.cov.wig')
        paragraph = Paragraph(ref_text, styles['BodyText'])
        elements.append(paragraph)

        # read wig file
        cov_df = pd.read_csv(wig_file, skiprows=3, delimiter='\t', index_col=0)
        cov = cov_df.sum(axis=1)

        cov.plot(kind='bar')
        # Add labels and title
        plt.xlabel('Base pair')
        plt.ylabel('Coverage')
        plt.title('Sum of rows')
        plt.savefig('figure.png')
        elements.append(Image('figure.png'))


    # build the doc
    doc.build(elements)

csv_file = '../RSV_run/3083501/Report/Report.csv'
working_folder = '../RSV_run/3083501/Report'
generate_pdf_report(csv_file, working_folder)