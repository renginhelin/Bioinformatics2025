import vcfpy
import os
import pandas as pd

def load_vcf_to_set(vcf_file):
    vcf_set = set()
    reader = vcfpy.Reader.from_path(vcf_file)
    for record in reader:
        vcf_set.add((record.CHROM, record.POS, record.REF, tuple(record.ALT)))
    return vcf_set

def calculate_metrics(ground_truth_set, predicted_set):
    TP = len(ground_truth_set & predicted_set)  # True positives: intersection
    FP = len(predicted_set - ground_truth_set)  # False positives: predicted but not in ground truth
    FN = len(ground_truth_set - predicted_set)  # False negatives: in ground truth but not predicted

    precision = TP / (TP + FP) if (TP + FP) > 0 else 0
    recall = TP / (TP + FN) if (TP + FN) > 0 else 0
    f1_score = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
    accuracy = TP / (TP + FP + FN) if (TP + FP + FN) > 0 else 0

    return TP, TN, FP, FN, precision, recall, f1_score, accuracy

ground_truth_vcf = "hc_bed_filtered.recode.vcf"
ground_truth_set = load_vcf_to_set(ground_truth_vcf)

vcf_directory = "/home/rengin/last"
vcf_files = [
    "filtered_mutect_bwa_NoBR.vcf", "filtered_mutect_bowtie_NoBR.vcf", 
    "filtered_strelka_bwa_NoBR.vcf", "filtered_strelka_bowtie_NoBR.vcf", "filtered_somatic_bwa_NoBR.vcf", 
    "filtered_somatic_bowtie_NoBR.vcf", "galaxy_filtered_bowtie.vcf",
    "filtered_mutect_bwa.vcf", "filtered_mutect_bowtie.vcf", "filtered_strelka_bwa.vcf", 
    "filtered_strelka_bowtie.vcf", "filtered_somaticsniper_bwa.vcf", "filtered_somaticsniper_bowtie.vcf"
]

results = []

for vcf_file in vcf_files:
    predicted_vcf_path = os.path.join(vcf_directory, vcf_file)
    
    if os.path.isfile(predicted_vcf_path):
        predicted_set = load_vcf_to_set(predicted_vcf_path)

        TP, TN, FP, FN, precision, recall, f1_score, accuracy = calculate_metrics(ground_truth_set, predicted_set)

        results.append({
            'VCF File': vcf_file,
            'TP': TP,
            'TN': TN,
            'FP': FP,
            'FN': FN,
            'Precision': f"{precision:.4f}",
            'Recall': f"{recall:.4f}",
            'F1 Score': f"{f1_score:.4f}",
            'Accuracy': f"{accuracy:.4f}"
        })
    else:
        print(f"File {vcf_file} not found in the directory {vcf_directory}")

df = pd.DataFrame(results)
print(df)
