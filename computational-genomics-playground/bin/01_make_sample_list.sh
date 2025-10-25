#!/usr/bin/env bash

set -euo pipefail # Exit on error (-e), treat unset vars as errors (-u), and fail if any command in a pipeline fails (pipefail)

# Check for correct number of arguments
if [[ $# -lt 4 ]]; then
  echo "Usage: $0 <panel_file> <super_pop: AFR|AMR|EAS|EUR|SAS> <n_samples> <out_samples_file>"
  echo "Example: $0 data/1kg/meta/integrated_call_samples_v3.20130502.ALL.panel EUR 500 data/1kg/meta/EUR_500.samples"
  exit 1
fi

PANEL="$1"
SUPER="$2"
NSAMP="$3"
OUT="$4"

mkdir -p "$(dirname "$OUT")" # Create output directory (and parents) if it doesn't exist

awk -v sp="$SUPER" '$3==sp {print $1}' "$PANEL" | shuf | head -n "$NSAMP" > "$OUT" 
# From the panel file, extract IDs in the chosen super-population, randomize, take N samples, and write to output
# set awk variable 'sp' to SUPER for comparison
# $3 is the super-pop column, $1 is the sample ID column
# ask awk to print the first column where the third column matches the super-population
# use shuf to randomize the order
# use head to take the first N samples

echo "Wrote sample list ($NSAMP from $SUPER) → $OUT"