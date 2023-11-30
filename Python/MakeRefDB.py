#!/usr/bin/env python3

import subprocess
import argparse

def build_star_reference(fasta_file, output_dir, num_threads=4):
    # Command to build STAR reference index
    star_command = [
        'STAR',
        '--runMode', 'genomeGenerate',
        '--genomeFastaFiles', fasta_file,
        '--genomeDir', output_dir,
        '--runThreadN', str(num_threads)
    ]

    # Run the STAR command
    subprocess.run(star_command)


if __name__ == "__main__":
    # Create command-line argument parser
    parser = argparse.ArgumentParser(description='Build STAR reference database from a Fasta file.')
    parser.add_argument('fasta_file', help='Path to the Fasta file containing the reference genome.')
    parser.add_argument('output_dir', help='Path to the output directory for STAR reference.')
    parser.add_argument('--num_threads', type=int, default=4, help='Number of threads to use (default: 4)')

    # Parse command-line arguments
    args = parser.parse_args()

    # Call the function to build STAR reference
    build_star_reference(args.fasta_file, args.output_dir, args.num_threads)
    
    # finished
    print("Reference building finished!\n")
