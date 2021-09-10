printf "\nStarting VCF files analysis...\n\n"

python3 run_mirnaome_analysis.py input_files/input_vcf/examplary_vcf_set \
output_folder/examplary_output_from_vcf \
reference_files/coordinates_hg38.bed \
reference_files/localizations_hg38.csv \
reference_files/coordinates_with_seq_hg38.bed

printf "\nStarting CSV file analysis...\n\n"

python3 run_mirnaome_analysis.py input_files/input_vcf/examplary_vcf_set \
output_folder/examplary_output_from_csv \
reference_files/coordinates_hg38.bed \
reference_files/localizations_hg38.csv \
reference_files/coordinates_with_seq_hg38.bed \
-c input_files/input_csv/examplary_csv_file/examplary_csv.csv

printf "\nAnalyses ended.\n\n"
