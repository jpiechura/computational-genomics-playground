#!/usr/bin/env bash
set -euo pipefail # Exit on error (-e), treat unset vars as errors (-u), and fail if any command in a pipeline fails (pipefail)

# Check for correct number of arguments
if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <chrom> 
  echo "Example: $0 22 
  exit 1
fi

CHR="$1"

# Define base FTP path and file names
BASE_FTP="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"
VCF="ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
TBI="${VCF}.tbi"
PANEL="integrated_call_samples_v3.20130502.ALL.panel"

# Download VCF, TBI, and panel files
# vcf file stores the variant calls
# tbi is the index file for the vcf
# panel file contains sample metadata
curl -L "${BASE_FTP}/${VCF}" -o "${VCF}"
curl -L "${BASE_FTP}/${TBI}"  -o "${TBI}"
curl -L "${BASE_FTP}/${PANEL}" -o "${PANEL}"

echo "Downloaded:"
echo " - ${VCF}"
echo " - ${TBI}"
echo " - ${PANEL}"