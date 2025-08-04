"""
Analysis pipeline for carbapenem‑resistant *Acinetobacter baumannii* isolates
==========================================================================

This script reproduces the descriptive and inferential statistics reported in
the accompanying manuscript. It operates entirely on the aggregate counts of
antibiotic resistance derived from 118 carbapenem‑resistant isolates collected
from two intensive care units in Tehran. The script computes resistance
percentages, produces summary figures (a bar chart for resistance by drug and
 a pie chart for colistin susceptibility), and performs chi‑square tests to
compare resistance profiles between hospitals and sequence types.

Usage
-----
Execute the script with a Python interpreter (version 3.8 or later). All
required dependencies are listed in ``requirements.txt``. For example:

.. code:: bash

   pip install -r requirements.txt
   python analysis.py

Running the script will output summary statistics to the terminal and write
two PNG files (``Resistance_BarChart.png`` and ``Colistin_PieChart.png``) to
the current working directory.

Notes
-----
The counts and proportions used here originate from the study described in
the manuscript. No patient‑level data are included; all computations are
based on aggregated numbers of resistant, intermediate and susceptible isolates.
"""

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt


def main() -> None:
    """Run the full analysis pipeline."""
    # List of antibiotics in the order they were tested
    antibiotics = [
        "Meropenem",
        "Imipenem",
        "Ceftriaxone",
        "Piperacillin-Tazobactam",
        "Ciprofloxacin",
        "Tobramycin",
        "Cefotaxime",
        "Gentamicin",
        "Ceftazidime",
        "Amikacin",
        "Cefepime",
        "TMP-SMX",
        "Doxycycline",
    ]

    # Total counts of resistant (R), intermediate (I) and susceptible (S) isolates
    # across both hospitals (Milad and Rasul Akram). These numbers come directly
    # from the study's susceptibility testing results.
    R_counts = np.array([118, 118, 118, 118, 118, 118, 115, 114, 114, 106, 99, 99, 45])
    I_counts = np.array([0, 0, 0, 0, 0, 0, 3, 1, 1, 2, 4, 5, 6])
    S_counts = 118 - R_counts - I_counts

    # Construct a DataFrame summarising resistance categories
    df_total = pd.DataFrame(
        {
            "Resistant (R)": R_counts,
            "Intermediate (I)": I_counts,
            "Susceptible (S)": S_counts,
        },
        index=antibiotics,
    )

    print("Overall resistance profile (118 isolates):")
    print(df_total)

    # Calculate percentage resistance for each antibiotic and summary statistics
    res_percent = (R_counts / 118) * 100
    avg_resistance = res_percent.mean()
    max_res_idx = int(np.argmax(res_percent))
    min_res_idx = int(np.argmin(res_percent))
    print(f"\nMean percent resistance across drugs: {avg_resistance:.1f}%")
    print(
        f"Highest resistance: {antibiotics[max_res_idx]} ({res_percent[max_res_idx]:.1f}%)"
    )
    print(
        f"Lowest resistance: {antibiotics[min_res_idx]} ({res_percent[min_res_idx]:.1f}%)"
    )

    # Generate bar chart of resistance percentages
    plt.figure(figsize=(10, 5))
    bars = plt.bar(antibiotics, res_percent, color="cornflowerblue")
    plt.xticks(rotation=60, ha="right")
    plt.ylabel("Percent of isolates resistant (%)")
    plt.title("Resistance to each antibiotic (n=118 isolates)")
    plt.ylim(0, 105)
    plt.tight_layout()
    plt.savefig("Resistance_BarChart.png", dpi=300)
    plt.close()

    # Pie chart for colistin susceptibility (from CBDE MIC results)
    colistin_counts = np.array([108, 7, 3])
    col_labels = ["Susceptible (91.5%)", "Intermediate (5.9%)", "Resistant (2.5%)"]
    plt.figure(figsize=(4, 4))
    plt.pie(
        colistin_counts,
        labels=col_labels,
        autopct="%1.1f%%",
        colors=["mediumseagreen", "gold", "tomato"],
    )
    plt.title("Colistin susceptibility distribution (n=118 isolates)")
    plt.savefig("Colistin_PieChart.png", dpi=300)
    plt.close()

    # Hospital-specific counts for Milad (68 isolates) and Rasul Akram (50 isolates)
    res_Milad = np.array([68, 68, 68, 68, 68, 68, 66, 66, 66, 60, 56, 57, 25])
    res_Rasul = np.array([50, 50, 50, 50, 50, 50, 49, 48, 48, 46, 43, 42, 20])
    sus_Milad = np.array([0, 0, 0, 0, 0, 0, 0, 2, 2, 7, 10, 8, 39])
    sus_Rasul = np.array([0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 5, 6, 28])
    int_Milad = 68 - res_Milad - sus_Milad
    int_Rasul = 50 - res_Rasul - sus_Rasul

    print("\nChi‑square tests comparing resistance between hospitals:")
    for idx, drug in enumerate(antibiotics):
        # Construct a 2x2 contingency table: rows = [resistant, non‑resistant];
        # columns = [Milad, Rasul Akram]. Non‑resistant = susceptible + intermediate.
        resistant = np.array([res_Milad[idx], res_Rasul[idx]])
        non_resistant = np.array([
            68 - res_Milad[idx],
            50 - res_Rasul[idx],
        ])
        contingency = np.array([resistant, non_resistant])
        # If an expected count would be zero the chi‑square test is not defined;
        # skip such cases or report that the comparison cannot be performed.
        if np.any(contingency == 0):
            print(f"{drug}: cannot compute chi‑square (zero expected counts)")
            continue
        chi2, p_val, _, _ = stats.chi2_contingency(contingency, correction=False)
        print(f"{drug}: p = {p_val:.3f}")

    # Chi‑square test for colistin susceptibility distribution between hospitals
    # Rows = hospitals; columns = categories [S, I, R]
    colistin_table = np.array([[62, 5, 1], [46, 2, 2]])
    chi2_col, p_col, _, _ = stats.chi2_contingency(colistin_table, correction=False)
    print(
        f"\nColistin susceptibility (2x3 table): chi^2 = {chi2_col:.3f}, p = {p_col:.3f}"
    )

    # Sequence type (ST) distribution across hospitals
    # Seven ST categories identified (incl. one novel type)
    ST_categories = ["ST218", "ST451", "ST1417", "ST3374", "ST391", "ST1104", "ST_new"]
    counts_milad = np.array([2, 2, 1, 1, 0, 0, 1])
    counts_rasul = np.array([1, 1, 0, 0, 1, 1, 0])
    ST_table = np.vstack([counts_milad, counts_rasul])
    chi2_st, p_st, _, _ = stats.chi2_contingency(ST_table, correction=False)
    print(
        f"ST distribution (2x7 table): chi^2 = {chi2_st:.3f}, p = {p_st:.3f}\n"
    )


if __name__ == "__main__":
    main()
