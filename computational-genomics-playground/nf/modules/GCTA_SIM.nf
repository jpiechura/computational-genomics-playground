

nextflow.enable.dsl=2

process GCTA_SIM {
  tag "gcta_sim_chr${params.chr}_${params.pop}${params.n_samples}"

  input:
    path bed
    path bim
    path fam

  publishDir "${params.phendir}", mode: 'copy', pattern: "sim.*"

  output:
    path "sim.phen", emit: pheno_file

  script:
'''
set -euo pipefail

bed_file="''' + bed + '''"
base="$(basename "$bed_file")"
base="${base%.bed}"

echo "Running GCTA GRM..."
gcta64 --bfile "$base" --make-grm --out sim --thread-num ''' + task.cpus + '''

# Count unique SNP IDs
nuniq=$(awk '!seen[$2]++{c++} END{print c+0}' "${base}.bim")
if [ "$nuniq" -lt ''' + params.n_causal + ''' ]; then
  echo "Not enough unique SNP IDs in ${base}.bim: have $nuniq, need ''' + params.n_causal + '''" >&2
  exit 1
fi

# Build causal list
set +o pipefail
awk '!seen[$2]++{print $2}' "${base}.bim" \\
  | awk -v seed=''' + params.seed + ''' 'BEGIN{srand(seed)} {printf "%0.12f\\t%s\\n", rand(), $0}' \\
  | LC_ALL=C sort -k1,1n -k2,2 \\
  | cut -f2 \\
  | head -n ''' + params.n_causal + ''' > causal.snplist
set -o pipefail
echo "Selected $(wc -l < causal.snplist) unique causal SNPs."

echo "Simulating phenotype..."
gcta64 --bfile "$base" \
       --simu-qt \
       --grm sim \
       --simu-hsq ''' + params.h2 + ''' \
       --simu-causal-loci causal.snplist \
       --simu-seed ''' + params.seed + ''' \
       --out sim --thread-num 2
'''
}