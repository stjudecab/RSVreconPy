import os
import re
import json
import subprocess
import pandas as pd
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq

# fetch record from JSON
def fetch_record_from_JSON(file_id, json_file):
    """Fetches the record by unique ID from the JSON file."""
    with open(json_file, "r", encoding="utf-8") as f:
        data = json.load(f)
    
    return data.get(file_id, "Record not found")

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
def processIGV(file_name, ref_file, cutoff):
    # Load the ref genome sequence
    ref_genome_record = next(SeqIO.parse(ref_file, "fasta"))
    ref_genome_sequence = str(ref_genome_record.seq)
    ref_genome_name = str(ref_genome_record.id)
    ref_length = len(ref_genome_sequence)

    sequence = ""
    sequence_array = ["N"] * ref_length

    with open(file_name, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('track') or line.startswith('#') or line.startswith('variableStep chrom='):
                continue
            elif line[0].isdigit():
                match = re.match(r'^(\d+)\t(\d+)\.0\t(\d+)\.0\t(\d+)\.0\t(\d+)\.0.+', line)
                if match:
                    pos, a, c, g, t = map(int, match.groups())
                    pos = pos - 1
                    
                    count_dict = {'A':a,'C':c,'G':g,'T':t}
                    cov = a + c + g + t

                    if cov >= cutoff:   # for high confidance position
                        max_key = max(count_dict, key=count_dict.get)
                        sequence_array[pos] = max_key
                    elif cov > 10 and any(x / cov > 0.9 for x in (a, c, t, g)):  # for position between 10 and cutoff, only trust the result if it's highly consisstant 
                        max_key = max(count_dict, key=count_dict.get)
                        sequence_array[pos] = max_key

                else:
                    print("error")
            else:
                print("error1")

    sequence = ''.join(sequence_array)

    return sequence

# determine the subtype
def determine_subtype(kma_out_file, reference_folder_name, cutoff=90):
    kma_df = pd.read_csv(kma_out_file, sep='\t')
    kma_df = kma_df.sort_values('Score', ascending=False)
    Selected_ref_name = kma_df.iloc[0, 0]
    template_identity = kma_df.iloc[0, 4]
    Selected_ref_name = re.sub(' ','',Selected_ref_name)
    #template_identity = re.sub(' ','',template_identity)

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

# parse gff file 
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

# find_gene_at_position 
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

# get genotype results
def get_genotype_res(path):
    with open(path, 'r') as file:
        lines = file.readlines()
        subtype_str = lines[0].split(',')
        return subtype_str

def last_pos_gff(gff_file):
    try:
        gff_df = pd.read_csv(gff_file, sep='\t', comment='#', header=None)

        # Assign column names
        gff_df.columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']

        # Filter for CDS entries
        cds_df = gff_df[gff_df['type'] == 'CDS']

        # Get the end position of the last CDS
        last_cds_end = cds_df['end'].iloc[-1]
    except:
        last_cds_end = 15000

    return int(last_cds_end)

# parse gff file and return ID
def parse_gff_return_id(gff_file):
    gene_positions = {}
    with open(gff_file, 'r') as file:
        for line in file:
            if not line.startswith('#'):
                fields = line.strip().split('\t')
                if len(fields) >= 9 and fields[2] == 'gene':
                    seq_id = fields[0]
                    start = int(fields[3])
                    end = int(fields[4])
                    ################################################
                    gene_name = fields[8].split(';')[-1].split('=')[1]
                    gene_positions[gene_name] = (seq_id, start, end)
    return gene_positions

# parse gff file and return ID
def parse_gff_return_CDSid(gff_file):
    gene_positions = {}
    with open(gff_file, 'r') as file:
        for line in file:
            if not line.startswith('#'):
                fields = line.strip().split('\t')
                if len(fields) >= 9 and fields[2] == 'CDS':
                    seq_id = fields[0]
                    start = int(fields[3])
                    end = int(fields[4])
                    ################################################
                    match = re.search(r'ID=([^;]+)', fields[8])
                    if match:
                        gene_name = match.group(1)
                    gene_name = fields[8].split(';')[-1].split('=')[1]
                    gene_positions[gene_name] = (seq_id, start, end)
    return gene_positions

def extract_sequence_from_fasta(fasta_file, start, end):
    for record in SeqIO.parse(fasta_file, "fasta"):
        return str(record.seq[start-1:end])

def align_sequences(ref_sequence, assembled_sequence):
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.open_gap_score = -1
    aligner.extend_gap_score = 0
    alignments = aligner.align(ref_sequence, assembled_sequence)
    best_alignment = alignments[0]  # Get the best alignment
    
    aligned_query_seq = best_alignment.aligned[1]   # get aligned region of query
    aligned_query_parts = ''.join([assembled_sequence[start:end] for start, end in aligned_query_seq])  # make aligned parts

    return aligned_query_parts

def translate_nt_to_aa(nucleotide_sequence):
    # Create a Seq object from the nucleotide sequence
    seq = Seq(nucleotide_sequence)
    
    # Translate the sequence to amino acids
    amino_acid_sequence = seq.translate()

    return str(amino_acid_sequence)

def extract_gene_seq(assembled_sequence_file, reference_sequence_file, gff_file, cds_name):
    
    # Load the assembled genome sequence
    assembled_genome_record = next(SeqIO.parse(assembled_sequence_file, "fasta"))
    assembled_genome_sequence = str(assembled_genome_record.seq)
    assembled_genome_name = str(assembled_genome_record.id)

    # Load the reference
    reference = next(SeqIO.parse(assembled_sequence_file, "fasta"))
    reference_sequence = str(reference.seq)
    reference_name = str(reference.id)

    if len(reference_sequence) == len(assembled_genome_sequence):
        # Parse GFF to get F gene coordinates
        gene_coordinates = parse_gff_return_CDSid(gff_file)

        start = gene_coordinates[cds_name][1]
        end = gene_coordinates[cds_name][2]

        extracted_sequence = assembled_genome_sequence[start-1:end]

    return assembled_genome_name, extracted_sequence


def extract_gene_seq_old(assembled_sequence_file, reference_sequence_file, gff_file, cds_name):
    # Parse GFF to get F gene coordinates
    gene_coordinates = parse_gff_return_CDSid(gff_file)

    start = gene_coordinates[cds_name][1]
    end = gene_coordinates[cds_name][2]

    # Extract F gene sequence from reference genome
    ref_gene_sequence = extract_sequence_from_fasta(reference_sequence_file, start, end)

    #print(ref_f_gene_sequence)

    # Load the assembled genome sequence
    assembled_genome_record = next(SeqIO.parse(assembled_sequence_file, "fasta"))
    assembled_genome_sequence = str(assembled_genome_record.seq)
    assembled_genome_name = str(assembled_genome_record.id)

    # Align the F gene sequence with the assembled genome
    alignment_query = align_sequences(ref_gene_sequence, assembled_genome_sequence)

    return assembled_genome_name, alignment_query

def detect_F_mutation(assembled_f_protein_sequence, mutation_file):
    # load mutation table
    df = pd.read_csv(mutation_file, index_col=1, header=None)
    df.columns = ['original', 'alternative']

    # strart screening
    screen_res = []
    for position in list(df.index):
        aa_cur_seq = assembled_f_protein_sequence[position - 1]
        if aa_cur_seq in df.loc[position]['original']:
            next
        else:
            #print(f"{df.loc[position]['original']}{position}{aa_cur_seq}")
            if aa_cur_seq == 'X':
                pass
            elif aa_cur_seq in df.loc[position]['alternative']:
                screen_res.append([f"{df.loc[position]['original']}{position}{aa_cur_seq}", 'Reported'])
            else:
                screen_res.append([f"{df.loc[position]['original']}{position}{aa_cur_seq}", 'Novel'])
        
    return screen_res

def extract_key_residue_fgene(fasta, key_mutation_file, output_csv):
    # Step 1: Read the FASTA file into a DataFrame
    sequence_names = []
    sequences = []
    for record in SeqIO.parse(fasta, "fasta"):
        sequence_names.append(record.id)  # Store the sequence name (header)
        sequences.append(list(str(record.seq)))  # Convert sequence to a list of amino acids

    seq_df = pd.DataFrame(sequences, index=sequence_names)      # Create a DataFrame
    seq_df.columns = range(1, seq_df.shape[1] + 1)  # Rename columns to start from 1

    # Step 2: Read the CSV file and extract positions
    data = pd.read_csv(key_mutation_file, header=None, names=["Residue", "Position", "Substitutions"])
    positions = data["Position"].astype(int).tolist()   # Extract positions as a list of integers

    # Step 3: Subset the DataFrame to keep only specified positions
    seq_df = seq_df[positions]

    # Step 4: Write the subset DataFrame to a CSV file
    seq_df.to_csv(output_csv, index=True, header=True)

def read_wig(wig_file):
    cov_df = pd.read_csv(wig_file, skiprows=3, delimiter='\t', index_col=0, header=None)
    cov_df = cov_df.iloc[:, :5]
    column_names = ['A', 'C', 'G','T','N']
    cov_df.columns = column_names
    cov = cov_df.sum(axis=1)

    return cov

def extract_gene_covarage(wig_file, gff_file):
    gene_coordinates = parse_gff_return_CDSid(gff_file)
    coverage = read_wig(wig_file)
    index_pool = list(coverage.index)

    gene_cov = {}
    for gene in gene_coordinates:
        start = gene_coordinates[gene][1]
        end = gene_coordinates[gene][2]
        gene_len = end - start + 1

        index_cur_gene = [x for x in index_pool if start <= x <= end]
        gene_average_coverage = sum(coverage[index_cur_gene].values) / gene_len

        gene_cov[gene] = gene_average_coverage

    return gene_cov

def fetch_file_by_type(folder, pattern):
    found_files = []
    # List only the files in the given directory (no recursion)
    for filename in os.listdir(folder):
        if filename.endswith(pattern):
            found_files.append(os.path.join(folder, filename))
    return found_files

def get_version(version_file):
    with open(version_file, "r") as file:
    # Read the first line
    first_line = file.readline().strip()  # strip() removes leading/trailing whitespace
    return first_line