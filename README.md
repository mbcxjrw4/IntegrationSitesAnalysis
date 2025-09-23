# IntegrationSitesAnalysis
---
## 🔬 Overview
A Nextflow-based bioinformatics pipeline to analyse lentiviral vector (LV) integration sites in T cell therapy products and post-infusion patient samples. This pipeline compares integration site patterns to public datasets and evaluates genotoxicity risk by assessing polyclonality, integration preferences, and genomic features.

---
## Features
- Insertion site sequence preference (sequence logo analysis)
- Gene-level integration patterns and insertion frequency
- Open chromatin analysis: Assess integration site overlap with open chromatin regions.
- GC content analysis: Calculate GC content in regions flanking integration sites.
- Clonality analysis: Evaluate polyclonality and potential oncogenic integration.
- Automated report: Generate an HTML report summarising all analyses.
---

## 📁 Project Structure
```bash
.
├── R/                     # R scripts for each analysis step
├── config/
│   ├── config.yaml        # Project parameters
│   └── load_config.R      # Load YAML config in R scripts
├── input/                  # Raw input data (intSites.tsv)
├── resource/               # Reference files used in analyses
├── results/                # Final outputs and plots
├── environment.yml         # Conda environment for reproducibility
└── workflows/
    └── main_workflow.nf   # Nextflow pipeline

```
---
## Installation
### Install Nextflow:
```bash
nextflow run workflows/main_workflow.nf --config config/config.yaml
```
### Create Conda environment:
```bash
conda env create -f environment.yml
conda activate integration-sites
```

---
## Usage
### Run the pipeline
```bash
nextflow run workflows/main_workflow.nf --config config/config.yaml
```
### Optional parameters (via config.yaml)
```bash
input: "input/intSites.tsv"
outdir: "results"
intermediate: "intermediate"
genome_fasta: "/path/to/GRCh38.p13.genome.fa"
report_rmd: "report.Rmd"
patient_id: "Patient001"
```
---
## Contributing
- Update config/config.yaml with new patient IDs or sample files.
- Ensure new R scripts follow the CLI argument pattern: --input, --output, --config.
- All intermediate files should use _updated suffix if modified.
---
## License
This project is licensed under the MIT License.
