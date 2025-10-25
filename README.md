# GWAS Demo Pipeline (Nextflow + Docker)

This repo demonstrates a minimal, reproducible setup for GWAS-style workflows using **Nextflow** with **containerized tools**.  
It assumes **system prerequisites** (Java, Nextflow, a container runtime) and then runs everything inside containers for reproducibility.

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

### Quick start
-Check prerequisites
```bash
java -version
nextflow -version
docker run hello-world
```
