# Wild Mouse TCR

> Code and data for the analyses accompanying **“The T-cell receptor repertoire of wild mice.”**  
The repository currently contains top-level folders **`Code/`** and **`Data/`**, plus a root `README.md`.

---

## Table of Contents

- [Repository structure](#repository-structure)
- [Inputs & expected formats](#inputs--expected-formats)
- [Reproducing the paper analyses](#reproducing-the-paper-analyses)
- [Outputs](#outputs)

---

## Repository structure

Wild-mouse-TCR/
├─ Code/        # R scripts / notebooks for data processing, analysis, and figure generation
├─ Data/        # Small, tracked data assets (e.g., metadata templates, lookups, toy examples)
├─ .gitignore
└─ README.md

---

## Inputs & expected formats

While raw data locations may differ by user, the analysis typically expects:

- **TCR tables** for α and/or β chains with (columns may vary):
  - sequence_id, v_call, j_call, junction_aa, duplicate_count
- **Sample metadata**:
  - mouse_id, sex, age, site

---

## Reproducing the paper analyses

The scripts in **Code/** are organized to: ingest TCR-seq tables, perform quality control and chain-specific processing, compute diversity and sharing metrics, and generate figures for the manuscript. Run them sequentially (or via a project script) to reproduce the results.

A typical run looks like:

1. 1_data_preparation.R  
   - Reads raw/processed TCR tables and sample metadata  
   - Collated raw data into summary tables for analysis

2. 2_data_exploration_visualisation.R  
   - Paper analysis, figures and supplementary plots (saved to Results/)

---

## Outputs

Running the pipeline will create all figures and analaysis presented in the manuscript.

