import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import pysam

def extract_variants(vcf_file):
    variants = set()
    vcf = pysam.VariantFile(vcf_file)
    for record in vcf.fetch():
        variant = (record.chrom, record.pos, record.ref, str(record.alts[0]))
        variants.add(variant)
    return variants

vcf_mapping = {
    "filtered_mutect_bwa_NoBR.vcf":"Mutect BWA No BR", 
    "filtered_mutect_bowtie_NoBR.vcf":"Mutect Bowtie No BR", 
    "filtered_strelka_bwa_NoBR.vcf":"Strelka BWA No BR", 
    "filtered_strelka_bowtie_NoBR.vcf":"Strelka Bowtie No BR", 
    "filtered_somatic_bwa_NoBR.vcf":"SomaticSniper BWA No BR", 
    "filtered_somatic_bowtie_NoBR.vcf":"SomaticSniper Bowtie No BR", 
    "galaxy_filtered_bowtie.vcf":"Galaxy Bowtie",
    "filtered_mutect_bwa.vcf":"Mutect BWA", 
    "filtered_mutect_bowtie.vcf":"Mutect Bowtie", 
    "filtered_strelka_bwa.vcf":"Strelka BWA", 
    "filtered_strelka_bowtie.vcf":"Strelka Bowtie", 
    "filtered_somaticsniper_bwa.vcf":"SomaticSniper BWA", 
    "filtered_somaticsniper_bowtie.vcf":"SomaticSniper Bowtie"
}
high_confidence_vcf = "hc_bed_filtered.recode.vcf"
high_confidence_name = "High Confidence"

all_variants = set()
vcf_variants = {}

for file, name in vcf_mapping.items():
    variants = extract_variants(file)
    vcf_variants[name] = variants
    all_variants.update(variants)

high_confidence_variants = extract_variants(high_confidence_vcf)
all_variants.update(high_confidence_variants)

all_variants = sorted(list(all_variants))  
binary_matrix = []

for name, variants in vcf_variants.items():
    row = [1 if variant in variants else 0 for variant in all_variants]
    binary_matrix.append(row)

high_confidence_row = [1 if variant in high_confidence_variants else 0 for variant in all_variants]
binary_matrix.append(high_confidence_row)

pipeline_names = list(vcf_mapping.values()) + [high_confidence_name]
df = pd.DataFrame(binary_matrix, columns=[f"Variant {i+1}" for i in range(len(all_variants))], index=pipeline_names)

pca = PCA(n_components=2)
pca_results = pca.fit_transform(df)

df_pca = pd.DataFrame(pca_results, columns=['PC1', 'PC2'], index=df.index)

print("PCA Coordinates:")
for name, coords in zip(df_pca.index, df_pca.values):
    print(f"{name}: PC1 = {coords[0]:.3f}, PC2 = {coords[1]:.3f}")

colors = [
    "red", "blue", "green", "purple", "orange", "brown", "pink", "gray", 
    "olive", "cyan", "teal", "gold", "lime", "black"
]
color_map = dict(zip(pipeline_names, colors))

plt.figure(figsize=(10, 7))
for i, label in enumerate(df.index):
    color = color_map[label]
    if label == high_confidence_name:
        plt.scatter(df_pca.iloc[i, 0], df_pca.iloc[i, 1], label=label, marker='x', s=200, c=color)
    else:
        plt.scatter(df_pca.iloc[i, 0], df_pca.iloc[i, 1], label=label, s=100, c=color)
              
plt.title('PCA of VCF Files and High-Confidence Variants')
plt.xlabel(f'Principal Component 1 ({pca.explained_variance_ratio_[0]:.2%})')
plt.ylabel(f'Principal Component 2 ({pca.explained_variance_ratio_[1]:.2%})')
plt.legend()
plt.grid()
plt.show()
