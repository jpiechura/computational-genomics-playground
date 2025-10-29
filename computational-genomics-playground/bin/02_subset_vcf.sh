#!/usr/bin/env bash

set -euo pipefail

# Check for correct number of arguments
if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <input_vcf.gz> <sample_list> <out_vcf.gz>"
  echo "Example: $0 data/1kg/vcf/ALL.chr22....vcf.gz data/1kg/meta/EUR_500.samples data/1kg/vcf/chr22.EUR500.vcf.gz"
  exit 1
fi

# Input arguments
# Input VCF file (bgzipped)
# Sample list file (one sample ID per line)
# Output VCF file (bgzipped)
INVCF="$1"
SAMPLES="$2"
OUTVCF="$3"

tmp="${OUTVCF%.vcf.gz}.tmp.vcf.gz"

# Create output directory if it doesn't exist
mkdir -p "$(dirname "$OUTVCF")"

# Subset VCF using bcftools
# -S specifies the sample list. sample variant status is represented in columns in the vcf.
# sample columns are all the columns after the FORMAT column
# -Oz outputs a bgzipped VCF
# -o specifies the output file
# last argument is the input VCF, which is bgzipped. bcftools will look for the .tbi index file automatically
bcftools view -S "$SAMPLES" -Oz -o "${tmp}" "$INVCF"

#add unique ids to variants with missing ids
bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -Oz -o "$OUTVCF" "${tmp}"

# Index the output VCF

bcftools index --tbi "$OUTVCF"

rm -f "${tmp}"

echo "Subset VCF → $OUTVCF (and .tbi)"