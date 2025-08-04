# Carbapenem‑Resistant *Acinetobacter baumannii* Analysis

This repository contains the analysis scripts and supporting files for the
study of multidrug‑resistant *Acinetobacter baumannii* isolates collected from
intensive care units in Tehran. The code reproduces the descriptive and
inferential statistics presented in the manuscript, including summary tables,
bar and pie charts of antibiotic resistance and susceptibility, and
contingency table tests comparing resistance profiles between two hospitals
(Milad and Rasul Akram) and the distribution of sequence types (STs).

## Repository structure

- `analysis.py` – executable Python script that performs the full analysis.
- `requirements.txt` – list of Python dependencies required to run the script.
- `LICENSE` – MIT licence covering the source code in this repository.
- `Resistance_BarChart.png` – bar chart of resistance percentages (created
  when running `analysis.py`).
- `Colistin_PieChart.png` – pie chart of colistin susceptibility categories
  (created when running `analysis.py`).

## Quick start

1. Ensure you have Python 3.8 or later installed.
2. Clone or download this repository.
3. Install the dependencies:

   ```bash
   pip install -r requirements.txt
   ```

4. Run the analysis script:

   ```bash
   python analysis.py
   ```

The script will print summary tables and test statistics to the terminal and
save two PNG figures in the working directory.

## Data provenance

The counts used in this analysis derive from aggregated susceptibility testing
of 118 carbapenem‑resistant *A. baumannii* isolates collected in 2024 from
patients in two Tehran hospitals (68 isolates from Milad and 50 from
Rasul Akram). Individual patient‑level data are not stored in this repository.

## Citation

If you use this code in your work, please cite the accompanying paper and
reference the specific Zenodo DOI associated with this release.
