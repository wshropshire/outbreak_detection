import os
import subprocess
import argparse
from pathlib import Path

# Function to parse command-line arguments
def parse_args():
    parser = argparse.ArgumentParser(description="Pipeline for identifying repeat regions, running Snippy, and Gubbins.")
    parser.add_argument('--reference', required=True, help='Reference genome fasta file')
    parser.add_argument('--prefix', required=True, help='Output prefix for generated files')
    parser.add_argument('--output_dir', required=True, help='Directory to store output files')
    parser.add_argument('--threads', required=False, default=16, type=int, help='Number of threads to use (default: 16)')
    return parser.parse_args()

# Step 1: Run NUCmer for All-vs-All comparison to identify repeat regions
def run_nucmer(reference, output_dir, prefix):
    nucmer_cmd = f"nucmer --maxmatch --nosimplify {reference} {reference} -p {output_dir}/{prefix}"
    subprocess.run(nucmer_cmd, shell=True)
    print("NUCmer All-vs-All comparison completed.")

# Step 2: Generate BED file of repeat regions
def generate_bed(delta_file, bed_file):
    show_coords_cmd = f"show-coords -r -T {delta_file} | awk '{{if ($1 != $3 && $2 != $4) print $0}}' | awk '{{print $8\"\\t\"$1\"\\t\"$2}}' > {bed_file}"
    subprocess.run(show_coords_cmd, shell=True)
    print(f"BED file {bed_file} generated.")

# Step 3: Mask repeat regions in the reference genome using BED file
def mask_reference(reference, bed_file, masked_reference):
    mask_cmd = f"bedtools maskfasta -fi {reference} -bed {bed_file} -fo {masked_reference}"
    subprocess.run(mask_cmd, shell=True)
    print(f"Masked reference genome saved as {masked_reference}.")

# Step 4: Prepare input tab file for Snippy
def create_snippy_input_tab(output_dir):
    names = []
    read1_files = []
    read2_files = []

    # Get file names for trimmed reads
    for read1 in sorted(Path().glob("*trim_1*")):
        name = read1.stem.split('_')[0]
        read2 = read1.with_name(read1.name.replace("trim_1", "trim_2"))

        names.append(name)
        read1_files.append(str(read1.resolve()))
        read2_files.append(str(read2.resolve()))

    # Write to input.tab
    input_tab_path = f"{output_dir}/input.tab"
    with open(input_tab_path, "w") as f:
        for name, read1, read2 in zip(names, read1_files, read2_files):
            f.write(f"{name}\t{read1}\t{read2}\n")

    print("Input tab file for Snippy created.")
    return input_tab_path

# Step 5: Run Snippy using the generated input.tab
def run_snippy(input_tab, masked_reference, threads):
    snippy_cmd = f"snippy-multi {input_tab} --cpus {threads} --ref {masked_reference} > runme.sh"
    subprocess.run(snippy_cmd, shell=True)
    subprocess.run("sh ./runme.sh", shell=True)
    print("Snippy pipeline completed.")

# Step 6: Clean Snippy alignment
def clean_snippy_alignment():
    clean_cmd = "snippy-clean_full_aln core.full.aln > clean.full.aln"
    subprocess.run(clean_cmd, shell=True)
    print("Snippy alignment cleaned.")

# Step 7: Run Gubbins for recombination analysis
def run_gubbins():
    gubbins_cmd = "run_gubbins.py -v -p VRE_masked_gubbins clean.full.aln"
    subprocess.run(gubbins_cmd, shell=True)
    print("Gubbins recombination analysis completed.")

# Step 8: Mask recombinant regions using a Python script
def mask_recombinant_regions():
    mask_gubbins_cmd = "python mask_gubbins_aln.py --aln clean.full.aln --gff VRE_masked_gubbins.recombination_predictions.gff --out NC_017022.1_gubbins_mask_MSA.fasta --out-fmt fasta"
    subprocess.run(mask_gubbins_cmd, shell=True)
    print("Recombinant regions masked.")

# Step 9: Create a fully masked file with 'N' for gaps
def create_total_masked_fasta():
    mask_total_cmd = f"awk '/^>Reference/ {{if (flag == 0) print $0; flag=1; next}} flag {{gsub(/-/, \"N\"); if (/^>/) {{flag=0; next}} else {{print}}}}' NC_017022.1_gubbins_mask_MSA.fasta > NC_017022.1_total_mask.fasta"
    subprocess.run(mask_total_cmd, shell=True)
    print("Total masked fasta file created.")

# Step 10: Check the number of Ns in the masked file
def check_masked_N_count():
    count_Ns_cmd = "awk '/^>/ {next} {n += gsub(/N/, \"N\") } END {print n}' NC_017022.1_total_mask.fasta"
    subprocess.run(count_Ns_cmd, shell=True)
    print("Checked the number of Ns in the total masked fasta file.")

# Step 11: Create mask BED file from masked FASTA
def create_mask_bedfile():
    bedfile_cmd = "python3 create_mask_bedfile.py NC_017022.1_total_mask.fasta final_mask_file.bed"
    subprocess.run(bedfile_cmd, shell=True)
    print("Mask BED file created.")

# Main execution workflow
def main():
    args = parse_args()

    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)

    delta_file = f"{args.output_dir}/{args.prefix}.delta"
    bed_file = f"{args.output_dir}/{args.prefix}_repeats.bed"
    masked_reference = f"{args.output_dir}/{args.prefix}_masked.fasta"

    run_nucmer(args.reference, args.output_dir, args.prefix)
    generate_bed(delta_file, bed_file)
    mask_reference(args.reference, bed_file, masked_reference)

    input_tab = create_snippy_input_tab(args.output_dir)
    run_snippy(input_tab, masked_reference, args.threads)

    clean_snippy_alignment()
    run_gubbins()
    mask_recombinant_regions()
    create_total_masked_fasta()
    check_masked_N_count()
    create_mask_bedfile()

if __name__ == "__main__":
    main()
