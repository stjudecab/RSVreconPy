import os
import re
import subprocess
import pandas as pd

# get sub folders for a dir
def get_sub_folders(folder_path):
    sub_folders = []

    with os.scandir(folder_path) as entries:
        for entry in entries:
            if entry.is_dir() and not entry.name.startswith(('.', '..')):
                sub_folders.append(entry.name)

    return sub_folders

# test if elements in array A in array B, if not, return elements that are not in B
def elements_not_in_array(array_a, array_b):
    not_in_b = [element for element in array_a if element not in array_b]
    return not_in_b

# sum the pct values
def pct_sum(*args):
    total = sum(float(item.replace('%', '')) for item in args)
    return total

# assemble sequences from IGV counts
def processIGV(file_name, cutoff):
    sequence = ""

    with open(file_name, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('track') or line.startswith('#') or line.startswith('variableStep chrom='):
                continue
            elif line[0].isdigit():
                match = re.match(r'^(\d+)\t(\d+)\.0\t(\d+)\.0\t(\d+)\.0\t(\d+)\.0.+', line)
                if match:
                    pos, a, c, g, t = map(int, match.groups())
                    
                    count_dict = {'A':a,'C':c,'G':g,'T':t}
                    cov = a + c + g + t

                    if cov > cutoff:
                        max_key = max(count_dict, key=count_dict.get)
                        sequence += max_key
                    else:
                        sequence += "N"
                else:
                    print("error")
            else:
                print("error1")

    return sequence

# determine the subtype
def determine_subtype(kma_out_file, reference_folder_name, cutoff=90):
    kma_df = pd.read_csv(kma_out_file, sep='\t')
    kma_df = kma_df.sort_values('Score', ascending=False)
    Selected_ref_name = kma_df.iloc[0, 0]
    template_identity = kma_df.iloc[0, 4]

    info_table = os.path.join(reference_folder_name, 'RSV.csv')
    info_df = pd.read_csv(info_table, delimiter=',', index_col=0)
    subtype = info_df.loc[Selected_ref_name,'Subtype']

    if float(template_identity) > cutoff:
        return f"{subtype},{Selected_ref_name},{subtype}"
    else:
        return f"Not RSV,{Selected_ref_name},{subtype}"

# check tool availability
def check_tool_availability(tool_name):
    try:
        # Use subprocess to run the command with the "--version" flag
        subprocess.run([tool_name, '--version'], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print(f"{tool_name} is available.")
        return 0
    except:
        try:
            # Use subprocess to run the command with the "-h" flag
            subprocess.run([tool_name, '-h'], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            print(f"{tool_name} is available.")
            return 0
        except:
            print(f"{tool_name} is not available. Please install it and try again.")
            return 1

# translate RNA sequence to DNA
def rna_to_dna(rna_seq):
    dna_seq = rna_seq.replace('U', 'T')
    return dna_seq

# parse gff file (not tested)
def parse_gff(gff_file):
    gene_positions = {}
    with open(gff_file, 'r') as file:
        for line in file:
            if not line.startswith('#'):
                fields = line.strip().split('\t')
                if len(fields) >= 9 and fields[2] == 'gene':
                    seq_id = fields[0]
                    start = int(fields[3])
                    end = int(fields[4])
                    gene_name = fields[8].split(';')[0].split('=')[1]
                    gene_positions[gene_name] = (seq_id, start, end)
    return gene_positions

# find_gene_at_position (not tested)
def find_gene_at_position(gene_positions, position):
    for gene_name, (seq_id, start, end) in gene_positions.items():
        if start <= position <= end:
            return gene_name
    return 'Not in CDS'

# find gff file under a path
def find_gff_files_in_path(path):
    # Get a list of all files and directories in the specified path
    entries = os.listdir(path)

    # Filter the entries to only include files with a ".gff" extension
    gff_files = [entry for entry in entries if os.path.isfile(os.path.join(path, entry)) and entry.endswith(".gff")]

    # Return the list of GFF files
    return gff_files

