# GWAS Demo Pipeline (Nextflow + Docker)

This repo demonstrates a minimal, reproducible setup for GWAS-style workflows using **Nextflow** with **containerized tools**.  In short, the pipeline downloads a subset of the **1000 genomes** Human Genome dataset based on inputs and filters the data for various QC criteria. Currently, **GCTA** is used to simulate a quantitative signal for demonstrative purposes. **PLINK** is used to test for association of individual SNPs with the simulated quantitative phenotype and produce a list of linked regions of statistically significant SNPs that are associated with the phenotype. Finally, **SuSiE-RSS** is used to fine-map SNPs within each linked region of SNPs to determine whether multiple independent causal variants are present and to compute posterior inclusion probabilities and credible sets. I plan to end some further plotting functionality to improve exploration of results.

The pipeline assumes **system prerequisites** (Java, Nextflow, a container runtime) and then runs everything inside containers for reproducibility.

## Why this repo?

- **Reproducible by default:** all heavy lifting runs inside containers.
- **Production-minded design:** clear inputs/outputs, automation-ready structure, and CI-ready layout.
- **Pedagogical:** meant to be easy to clone, run, and extend.

## At a glance

- **Orchestrator:** Nextflow (≥ 24.x)
- **Runtime:** Docker Desktop / Colima / Podman
- **Language mix:** Nextflow + R + shell + a touch of Python for utilities
- **Typical stages:** ingest → QC → covariates/PCA → association → fine-mapping (e.g., SuSiE-RSS)

## Pipeline diagrame
```mermaid
flowchart LR
A --> B;
B --> C;


## Prerequisites

- **Java** 11–21 (recommend Temurin 17)
- **Nextflow** ≥ 24.x
- **Container runtime**: Docker Desktop (recommended) or Colima/Podman

### Install Java
- Download and install: https://adoptium.net/temurin/releases/?version=17  
- If you run on a mac, it's easy to install manually than via homebrew
- Verify:
  ```bash
  java -version

### Install Nextflow
```bash
curl -s https://get.nextflow.io | bash
mkdir -p ~/bin && mv nextflow ~/bin/
echo 'export PATH="$HOME/bin:$PATH"' >> ~/.zshrc
exec $SHELL -l
nextflow -version
```

### Install a container runtime
Option A — Docker Desktop (recommended):

-Download and install: https://www.docker.com/products/docker-desktop/
-Verify:
```bash
docker version
docker run hello-world
```

Option B — Colima (lightweight, OSS):
-If you have homebrew:
```bash
brew install colima docker
colima start --cpu 4 --memory 8 --disk 60
docker run hello-world
```

### Install a container runtime

### Quick start
-Check prerequisites
```bash
java -version
nextflow -version
docker run hello-world
```
