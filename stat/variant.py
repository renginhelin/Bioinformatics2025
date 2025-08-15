import vcfpy
import os
import matplotlib.pyplot as plt
import pandas as pd

vcf_files = ["filtered_mutect_bwa_NoBR.vcf", "filtered_mutect_bowtie_NoBR.vcf", 
             "filtered_strelka_bwa_NoBR.vcf", "filtered_strelka_bowtie_NoBR.vcf", "filtered_somatic_bwa_NoBR.vcf", 
             "filtered_somatic_bowtie_NoBR.vcf", "filtered_mutect_bwa.vcf", "filtered_mutect_bowtie.vcf", 
             "filtered_strelka_bwa.vcf", "filtered_strelka_bowtie.vcf", "filtered_somaticsniper_bwa.vcf", 
             "filtered_somaticsniper_bowtie.vcf", "galaxy_filtered_bowtie.vcf"]

def count_variants_in_vcf(vcf_file):
    reader = vcfpy.Reader.from_path(vcf_file)
    return sum(1 for _ in reader)

variant_counts = []

for vcf_file in vcf_files:
    count = count_variants_in_vcf(vcf_file)
    print(f"File: {vcf_file}, Variant Count: {count}")
    variant_counts.append(count)

df = pd.DataFrame({
    'VCF File': vcf_files,
    'Variant Count': variant_counts
})

plt.figure(figsize=(10, 6))
plt.bar(df['VCF File'], df['Variant Count'], color='skyblue')
plt.xlabel('VCF File')
plt.ylabel('Variant Count')
plt.title('Variant Count in Each VCF File')
plt.xticks(rotation=90)  
plt.tight_layout()
plt.show()
