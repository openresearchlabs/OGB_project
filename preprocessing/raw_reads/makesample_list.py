import os
import csv
import sys

def create_samplesheet(input_dir, output_file):
    with open(output_file, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["sample", "run_accession", "instrument_platform", "fastq_1", "fastq_2", "fasta"])

        for root, dirs, files in os.walk(input_dir):
            for dir_name in dirs:
                sample = dir_name
                dir_path = os.path.join(root, dir_name)
                fastq_files = [f for f in os.listdir(dir_path) if f.endswith('_1.fq.gz')]
                for fq1 in fastq_files:
                    fq2 = fq1.replace('_1.fq.gz', '_2.fq.gz')
                    if fq2 in os.listdir(dir_path):
                        run_accession = '_'.join(fq1.split('_')[1:4])
                        fastq_1 = os.path.join(input_dir,sample, fq1)
                        fastq_2 = os.path.join(input_dir,sample, fq2)
                        instrument_platform = "ILLUMINA"
                        writer.writerow([sample, run_accession, instrument_platform, fastq_1, fastq_2, ""])

    print(f"Samplesheet created: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python create_samplesheet.py <input_dir> <output_file>")
        sys.exit(1)

    input_dir = sys.argv[1]
    output_file = sys.argv[2]

    create_samplesheet(input_dir, output_file)
