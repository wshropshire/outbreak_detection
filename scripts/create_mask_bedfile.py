#!/usr/bin/env python3

import sys
from Bio import SeqIO

def create_mask_bed(input_file, output_file, mask_char='N'):
    # Detect file format based on extension
    file_format = "genbank" if input_file.endswith(".gbk") else "fasta"

    with open(output_file, "w") as bed:
        for record in SeqIO.parse(input_file, file_format):
            seq = str(record.seq)
            header = record.id
            
            in_mask_block = False
            start_pos = None
            
            for i, base in enumerate(seq):
                if base == mask_char:
                    if not in_mask_block:
                        start_pos = i
                        in_mask_block = True
                else:
                    if in_mask_block:
                        end_pos = i
                        bed.write(f"{header}\t{start_pos}\t{end_pos}\n")
                        in_mask_block = False
            
            # If the sequence ends with the mask character, make sure to close the last block
            if in_mask_block:
                bed.write(f"{header}\t{start_pos}\t{len(seq)}\n")

    print(f"BED file created: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: create_mask_bed.py <input_file> <output_file> [<mask_char>]")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    mask_char = sys.argv[3] if len(sys.argv) > 3 else 'N'

    create_mask_bed(input_file, output_file, mask_char)

