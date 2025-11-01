nextflow.enable.dsl=2

process MAKE_PY_ENV {   
  publishDir "${params.workdir}/pyenv", mode: 'copy'
  container 'python:3.11'

  output:
    path 'pyenv'  // directory

  script:
  """
  python -m venv pyenv
  ./pyenv/bin/python -m pip install --upgrade pip >/dev/null
  ./pyenv/bin/python -m pip install --no-cache-dir -q numpy pandas
  # quick sanity check
  ./pyenv/bin/python - <<'PY'
  import numpy, pandas
  PY
  """
}
