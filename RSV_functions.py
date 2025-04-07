import os
import re
import json
import shutil
from io import StringIO
import subprocess
import pandas as pd
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from collections import defaultdict

# detect sequencing files
def detect_sequencing_files(data_folder):
    """
    Detect sequencing files in a directory, compatible with various common naming conventions
    
    Args:
        data_folder: Path to the directory containing sequencing files
        
    Returns:
        A dictionary where keys are sample IDs and values are lists containing R1 and R2 file paths
    """
    data_folder = os.path.normpath(data_folder)
    sample_dict = defaultdict(list)
    
    # Supported sequencing file extensions
    extensions = ('.fastq.gz', '.fq.gz', '.fastq', '.fq')
    
    # Common paired-end sequencing file naming patterns
    patterns = [
        # Format: sample_R1.fastq.gz, sample_R2.fastq.gz
        re.compile(r'(.+)_R?([12])[_.].*'),
        # Format: sample.1.fastq.gz, sample.2.fastq.gz
        re.compile(r'(.+)\.([12])\..*'),
        # Format: sample_1.fastq.gz, sample_2.fastq.gz
        re.compile(r'(.+)_([12])\..*'),
        # Format: sample-read1.fastq.gz, sample-read2.fastq.gz
        re.compile(r'(.+)[_-]read?([12])[_.].*'),
        # Format: sample_S1_L001_R1_001.fastq.gz (Illumina format)
        re.compile(r'(.+)_S\d+_L\d+_R?([12])_\d+.*')
    ]
    
    with os.scandir(data_folder) as entries:
        for entry in entries:
            if entry.name.endswith(extensions) and entry.is_file():
                file = entry.name
                file_path = os.path.join(data_folder, file)
                
                # Try to match against various naming patterns
                matched = False
                for pattern in patterns:
                    match = pattern.fullmatch(file)
                    if match:
                        sample_id = match.group(1)
                        read_num = match.group(2)
                        
                        # Ensure read_num is either 1 or 2
                        if read_num not in ['1', '2']:
                            continue
                            
                        # Initialize sample entry if not exists
                        if sample_id not in sample_dict:
                            sample_dict[sample_id] = [None, None]
                        
                        # Store file path
                        sample_dict[sample_id][int(read_num)-1] = file_path
                        matched = True
                        break
                
                # Fallback to original method if no patterns matched
                if not matched:
                    if '_R1' in file or '_1.' in file:
                        sample_id = re.split(r'_R1|_1\.', file)[0]
                        if sample_id not in sample_dict:
                            sample_dict[sample_id] = [None, None]
                        sample_dict[sample_id][0] = file_path
                    elif '_R2' in file or '_2.' in file:
                        sample_id = re.split(r'_R2|_2\.', file)[0]
                        if sample_id not in sample_dict:
                            sample_dict[sample_id] = [None, None]
                        sample_dict[sample_id][1] = file_path
    
    # Remove incomplete samples (only R1 or only R2 present)
    complete_samples = {k: v for k, v in sample_dict.items() if all(v)}
    
    return complete_samples

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

    info_table = os.path.join(reference_folder_name, 'RSV.csv')
    info_df = pd.read_csv(info_table, delimiter=',', index_col=0)
    subtype = info_df.loc[Selected_ref_name,'Subtype']

    if float(template_identity) > cutoff:
        return f"{subtype},{Selected_ref_name},{subtype}"
    else:
        return f"Not RSV,{Selected_ref_name},{subtype}"

# check tool availability
def check_tool_availability(tool_name):
    # First check if the tool exists in PATH
    if shutil.which(tool_name) is None:
        print(f"{tool_name} is not installed or not in system PATH.")
        return 1
    
    # Special handling for BWA
    if tool_name.lower() == "bwa":
        try:
            # BWA shows version info when run without arguments (to stderr)
            result = subprocess.run([tool_name], check=False,
                                  stdout=subprocess.PIPE, 
                                  stderr=subprocess.PIPE,
                                  text=True)
            if "version" in result.stderr.lower():
                print(f"{tool_name} is available.")
                return 0
            raise Exception("BWA didn't show version info")
        except:
            print(f"{tool_name} is present but couldn't verify version info.")
            return 1
    
    # For all other tools
    try:
        # First try with no arguments
        subprocess.run([tool_name], check=True, 
                      stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print(f"{tool_name} is available.")
        return 0
    except:
        pass
    
    # If that fails, try various version/help flags
    for flag in ['--version', '-v', '-h', '--help', '-V', 'version']:
        try:
            subprocess.run([tool_name, flag], check=True,
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            print(f"{tool_name} is available (responded to {flag}).")
            return 0
        except:
            continue
    
    print(f"{tool_name} is present but couldn't verify it works properly.")
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
    single_mutation_index = []
    single_mutation_array = []
    co_mutation_array = []
    with open(mutation_file, 'r') as file:
        for line in file:
            cur_line = line.strip()
            if "+" in cur_line:
                # co-mutations
                ele = cur_line.split('+')
                co_mutation_array.append(ele)
            else:
                # single mutations
                ele = cur_line.split(',')
                if len(ele) == 3:
                    single_mutation_array.append([ele[0], ele[2]])
                    single_mutation_index.append(ele[1])

    screen_res = []
    # strart screening single mutations
    single_mutation_df = pd.DataFrame(single_mutation_array, index = single_mutation_index, columns = ['original', 'alternative'])
    for position in list(single_mutation_df.index):
        aa_cur_seq = assembled_f_protein_sequence[int(position) - 1]
        if aa_cur_seq in single_mutation_df.loc[position]['original']:
            next
        else:
            #print(f"{df.loc[position]['original']}{position}{aa_cur_seq}")
            if aa_cur_seq == 'X':
                pass
            elif aa_cur_seq in single_mutation_df.loc[position]['alternative']:
                screen_res.append([f"{single_mutation_df.loc[position]['original']}{position}{aa_cur_seq}", 'Reported'])
            else:
                screen_res.append([f"{single_mutation_df.loc[position]['original']}{position}{aa_cur_seq}", 'Novel'])
    
    # strart screening co-mutations
    for ele in co_mutation_array:
        co_mutation_names = []
        detect_sign = 1
        for mutations in ele:
            cur_mutation_name = mutations.replace(",", "")
            co_mutation_names.append(cur_mutation_name)
            sub_ele = mutations.split(',')
            aa_cur_seq = assembled_f_protein_sequence[int(sub_ele[1]) - 1]
            if aa_cur_seq not in sub_ele[2]:
                detect_sign = 0
                break
        if detect_sign == 1:
            co_mutation_name = '+'.join(co_mutation_names)
            screen_res.append([co_mutation_name, 'Reported'])

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
    # Read the file content and replace multiple line separators
    with open(key_mutation_file, 'r') as file:
        content = file.read().replace('+', '\n')
    cleaned_file = StringIO(content)

    data = pd.read_csv(cleaned_file, header=None, names=["Residue", "Position", "Substitutions"])
    positions = data["Position"].astype(int).tolist()   # Extract positions as a list of integers

    # Step 3: Subset the DataFrame to keep only specified positions
    unique_sorted_positions = sorted(set(positions))
    seq_df = seq_df[unique_sorted_positions]

    # Step 4: Write the subset DataFrame to a CSV file
    seq_df.to_csv(output_csv, index=True, header=True)


def read_wig(wig_file):
    cov_df = pd.read_csv(wig_file, skiprows=3, delimiter='\t', index_col=0, header=None, usecols=range(6))
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

