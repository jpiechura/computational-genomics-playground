nextflow.enable.dsl=2

process BUILD_COVAR {
  tag "build_covar.${params.chr}_${params.pop}${params.n_samples}"

  input:
    path eigenvec
    path panel
    path append_sex_script

  output:
    path "*.covar",        emit: covar_file

  publishDir "covar", mode: 'copy'

  script:
"""
  set -euo pipefail

  OUT='${eigenvec.baseName}.with_sex.covar'
  AWK_SCRIPT='${append_sex_script}'

  # Defensive: normalize line endings in case the awk file has CRLF
  tr -d '\\r' < "\$AWK_SCRIPT" > __append_sex.awk && mv __append_sex.awk "\$AWK_SCRIPT"

  # Sanity checks
  echo "Using awk: \$(awk --version 2>/dev/null | head -1 || true)" >&2
  echo "AWK script path: \$AWK_SCRIPT" >&2
  head -5 "\$AWK_SCRIPT" >&2

  awk -v sc='1' \
      -v xc='4' \
      -f "\$AWK_SCRIPT" \
      '${panel}' '${eigenvec}' > "\$OUT"

  awk '(\$NF==1){m++} (\$NF==2){f++} (\$NF==-9){u++} \
       END{printf("SEX counts  male:%d  female:%d  unknown:%d\\n", m+0,f+0,u+0)}' \
       "\$OUT" >&2
  """
}

