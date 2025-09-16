# IntegrationSitesAnalysis
---
## 🔬 Overview
This pipeline analyses lentiviral vector (LV) integration sites from cell therapy drug products and post-infusion samples, focusing on genotoxicity risk assessment. It performs a comprehensive comparison with public datasets and evaluates safety based on: 
- Insertion site sequence preference (sequence logo analysis) 
- Gene-level integration patterns and insertion frequency 
- Overlap with open chromatin regions 
- Analysis of GC content at integration sites 
- Clonality assessment of vector-containing T cells 
- The output includes statistical summaries and a complete safety report for regulatory or scientific review.

---

## 📁 Project Structure
```bash
IntegrationSitesAnalysis/ 
├── R/                         # R scripts for each analysis module 
├── config/                    # Configuration file (YAML) 
├── input/                     # Input integration site data 
├── resource/                  # Reference files (genes, chromatin, GC, etc.) 
├── utilities/                 # Supporting scripts for data prep 
└── protocol.sh                # Driver script to execute the full pipeline
```
