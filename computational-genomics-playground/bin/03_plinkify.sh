#!/usr/bin/env bash
set -euo pipefail

# Check for correct number of arguments
if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <input_vcf.gz> <out_prefix>"
  echo "Example: $0 data/1kg/vcf/chr22.EUR500.vcf.gz data/1kg/plink/chr22.EUR500"
  exit 1
fi


VCF="$1"
OUTPREFIX="$2"

# Create output directory if it doesn't exist
mkdir -p "$(dirname "$OUTPREFIX")"

# Convert VCF to PLINK format. this generates three files:
# .bed : binary genotype table
# .bim : variant information
# .fam : sample information
plink --vcf "$VCF" --make-bed --out "$OUTPREFIX"

echo "PLINK files → ${OUTPREFIX}.bed/.bim/.fam"