# ALPINE Manuscript Analysis Scripts

This repository contains analysis scripts and visualization code for reproducibility of results presented in the ALPINE manuscript: "ALPINE: a comprehensive tool for CRISPR gene editing outcome classification from long-read sequencing data".

## Repository Structure

```
ALPINE-manuscript-analysis/
├── README.md                           # This file
├── simulation_generation/              # Scripts for generating simulated datasets
├── simulation_benchmark/              # Scripts for benchmarking tools on simulated data  
├── public_dataset_analysis/           # Scripts for analyzing PRJNA913199 public dataset
```

## Prerequisites

### Software Requirements

- **Python 3.8+** with packages: pandas, numpy, pysam, matplotlib, seaborn, biopython
- **R 4.0+** with packages: ggplot2, dplyr, tidyr, RColorBrewer, stringr, purrr, broom
- **External Tools**:
  - [ALPINE](https://github.com/Maggi-Chen/ALPINE) - CRISPR outcome classification tool
  - [knock-knock](https://github.com/jeffhussmann/knock-knock) - CRISPR outcome analysis
  - [CRISPResso2](https://github.com/pinellolab/CRISPResso2) - CRISPR analysis toolkit
  - [badread](https://github.com/rrwick/Badread) - Long-read sequencing simulator

## Script Descriptions

### Simulation Generation (`simulation_generation/`)
- **`download_ref_files.sh`** - Downloads reference genome and annotation files
- **`generate_integration_template.py`** - Creates template sequences for integration simulations
- **`generate_variant_template.py`** - Generates variant template sequences for simulation
- **`get_ha_seq.py`** - Extracts homology arm sequences from reference
- **`get_itr_seq.py`** - Extracts inverted terminal repeat (ITR) sequences
- **`merge_simulated_classes_ont.py`** - Combines simulated Oxford Nanopore data classes
- **`merge_simulated_classes_pacbio.py`** - Combines simulated PacBio HiFi data classes  
- **`simulate_hifi.sh`** - Runs PacBio HiFi read simulation using badread
- **`simulate_ont.sh`** - Runs Oxford Nanopore read simulation using badread

### Simulation Benchmarking (`simulation_benchmark/`)
- **`check_accuracy_alpine.py`** - Calculates precision, recall, and F1 scores for ALPINE
- **`check_accuracy_knock.py`** - Calculates precision, recall, and F1 scores for knock-knock
- **`plot_simulation_benchmark.R`** - Generates Figure 2 and simulation benchmark plots

### Public Dataset Analysis (`public_dataset_analysis/`)
- **`count_9del_alpine.py`** - Counts 9bp deletion events in ALPINE output
- **`count_9del_knock.py`** - Counts 9bp deletion events in knock-knock output
- **`crispresso2_statistical_test.R`** - Performs statistical analysis on CRISPResso2 results
- **`get_read_proportion_alpine_knock_crispresso2.py`** - Extracts outcome proportions from all three tools
- **`get_readname_integration.py`** - Identifies integration event read names for validation
- **`proportion_knock_from_per_sample_out.py`** - Processes per-sample knock-knock output files
- **`visualize_supp_figures.R`** - Generates Supplementary Figures S7-S11

## Data Sources

- **PRJNA913199**: Public dataset from "Interstrand crosslinking of homologous repair template DNA enhances gene editing in human cells" (Nature Biotechnology, 2023)
- **Simulated data**: Generated using badread based on reference templates

## Notes

- No raw data files are included in this repository
- Users must download the PRJNA913199 dataset independently from NCBI SRA
- Simulated data is generated on-the-fly using the provided scripts
- All analysis scripts expect data in standard formats (FASTQ, BAM, etc.)

## Citation

If you use these scripts, please cite the ALPINE preprint:

```
ALPINE: a comprehensive tool for CRISPR gene editing outcome classification from long-read sequencing data.
bioRxiv 2026.03.27.714831; doi: https://doi.org/10.64898/2026.03.27.714831
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contact

For questions about these analysis scripts, please open an issue on this repository.

---

**Note**: This repository contains analysis scripts only. For the ALPINE tool itself, please visit the [main ALPINE repository](https://github.com/Maggi-Chen/ALPINE).