# Refining DNA Sequence Alignment Accuracy in HAlign through Gradient Descent Techniques

> **PeerJ Computer Science Submission** | Manuscript: 142426  
> **Authors:** Ria K, Rhema Yashaswini O, K N Meera, T M Mamatha, Yuqing Lin  
> **Corresponding Author:** K N Meera — kn_meera@blr.amrita.edu

---

## Description

This repository contains the complete source code and simulated dataset for the paper *"Refining DNA Sequence Alignment Accuracy in HAlign through Gradient Descent Techniques"* submitted to PeerJ Computer Science.

The project proposes an improved DNA sequence alignment framework built on top of HAlign that uses:
- **Seed-based filtering** — fast k-mer anchor detection to narrow down candidate alignment regions
- **Banded Smith-Waterman local alignment** — dynamic programming restricted to high-likelihood regions for efficiency
- **Gradient descent optimization** — automated tuning of alignment parameters (gap penalty, match score, mismatch penalty) to maximize accuracy
- **Genetic algorithm** — population-based search combined with gradient descent to escape local optima

The improved aligner is benchmarked against Old HAlign and a simulated T-Coffee implementation across three species: *Homo sapiens*, *Arabidopsis thaliana*, and *Neurospora crassa*, using accuracy, precision, recall, and F1-score as evaluation metrics.

---

## Repository Structure

```
Refining-DNA-Sequence-Alignment/
├── dna_sequence.py      # Complete alignment implementation:
│                        #   - Seed-based filtering
│                        #   - Banded Smith-Waterman alignment
│                        #   - Gradient descent parameter optimization
│                        #   - Genetic algorithm
│                        #   - Heatmap, 3D hyperboloid, and helix visualization
└── README.md            # This file
```

---

## Dataset Information

All DNA sequences used in this study are **simulated/synthetic** sequences generated programmatically within `dna_sequence.py`. No external biological database downloads are required.

### Reference DNA Sequence (used in experiments)

```
ATGCGTACGTTAGCTAGCTAGTCGAATCGTAGCTAGATCGTAGCTTACGTAGCTAGCTAGCTAGCTAGGATCGTAGCTAGCTAGTCGATCGATCGTAGCTAGCTAGCTAGCTAG
CTAGCATCGTAGCATCGTAGCTAGCTAGCTAGCTAGCTAGCTTAGCTAGCATCGTAGCTAGCTAGCTAGCTAGCTAGCTAGCATCGTAGCATCGTACGATCGTAGCTAGCTAGT
CGATCGTAGCTAGCTAGCATCGTAGCTAGCTAGCTAGCTAGCTGACTGACTGACTGACTGACTGACTGACTGACTGACGATCGTAGCTAGCTAGCTAGCTAGCATCGTAGCTAG
CTAGCTAGCTAGCTAGCATCGTAGCTAGCTAGCTAGCTAGCTAGCCGATCGATCGCTAGCTAGCATCGTAGCATCGTA
```

### Query DNA Sequence

```
GCTAGTCGAATCGTAGCTAGATCGTAGCTTACGTAGCTAGCTAGCTAGCTAG
```

### Species Datasets (Simulated)

| Species | Reference Length | Query Length | Notes |
|---|---|---|---|
| *Homo sapiens* | ~500 bp | ~52 bp | Synthetic, GC-rich |
| *Arabidopsis thaliana* | ~500 bp | ~52 bp | Synthetic, AT-rich |
| *Neurospora crassa* | ~500 bp | ~52 bp | Synthetic, balanced |

All sequences are randomly generated with controlled GC-content per species to simulate realistic alignment conditions.

---

## Code Information

### File: `dna_sequence.py`

A single self-contained Python script implementing the full alignment pipeline:

| Component | Description |
|---|---|
| `seed_filter()` | Scans query for k-mer seeds; filters by frequency threshold |
| `banded_smith_waterman()` | Local alignment within a band window around each anchor |
| `gradient_descent_optimize()` | Iterative parameter tuning: gap open, gap extend, match, mismatch |
| `genetic_algorithm()` | Population-based search over parameter space with crossover + mutation |
| `evaluate_metrics()` | Computes accuracy, precision, recall, F1-score |
| `plot_hyperboloid()` | 3D surface plot of alignment scores across reference-query matrix |
| `plot_heatmap()` | Comparative heatmap of accuracy/precision/recall/F1 across methods |
| `plot_helix()` | Text-based DNA double helix visualization |

### Alignment Score Matrix (Smith-Waterman)

```
S(i, j) = max { 0,
                S(i-1, j-1) + s(q_i, r_j),
                S(i-1, j)   + g,
                S(i, j-1)   + g }
```

where `s(q_i, r_j)` is the match/mismatch score and `g` is the gap penalty.

### Gradient Descent Update Rule

```
p_i  ←  p_i + η · gradient_i
```

Parameters optimized: `p1` = gap opening, `p2` = gap extension, `p3` = match score, `p4` = mismatch penalty.

### Fitness Function

```
F(i) = AlignmentAdjustment × [0.4·Sid + 0.5·Ecomp + 0.1·Bcov]
AlignmentAdjustment = 1 − 0.1 × AlignmentRiskLevel
```

---

## Requirements

### System Requirements

| Dependency | Version | Purpose |
|---|---|---|
| Python | 3.8+ | Core runtime |
| NumPy | 1.21+ | Matrix operations, dynamic programming |
| Matplotlib | 3.4+ | Heatmaps, 3D plots, visualizations |
| SciPy | 1.7+ (optional) | Advanced optimization utilities |

### Installation

Install all dependencies with pip:

```bash
pip install numpy matplotlib scipy
```

No additional downloads or database setup required. All sequences are generated within the script.

---

## Usage Instructions

### Step 1 — Clone the repository

```bash
git clone https://github.com/Ria2826/Refining-DNA-Sequence-Alignment
cd Refining-DNA-Sequence-Alignment
```

### Step 2 — Install dependencies

```bash
pip install numpy matplotlib scipy
```

### Step 3 — Run the alignment pipeline

```bash
python dna_sequence.py
```

### Expected Output

Running the script produces:

1. **Console output** — alignment scores, best match position, coverage, and gap count for each method
2. **Figure 1** — 3D hyperboloid surface plot showing alignment score across the reference-query matrix, with peak alignment marked
3. **Figure 2** — Heatmap comparing accuracy, precision, recall, and F1-score for Improved Aligner vs. Old HAlign vs. T-Coffee
4. **Figure 3** — Text-based DNA double helix visualization of the query sequence alignment

### Sample Console Output

```
=== Improved Aligner ===
Alignment Score : 104
Start Position  : 16
Coverage        : 52
Gap Count       : 3

=== Old HAlign (Simulated) ===
Heuristic Seed Score : 400

=== T-Coffee (Simulated) ===
Iterative Score  : 104

=== Performance Metrics ===
Method              Accuracy  Precision  Recall   F1
Improved Aligner    0.92      0.91       0.93     0.92
Old HAlign          0.85      0.83       0.87     0.85
T-Coffee Sim.       0.89      0.88       0.90     0.89
```

---

## Methodology

### Pipeline Overview

```
INPUT: Reference DNA + Query DNA
       ↓
Step 1: Seed Filtering
        → Scan query for k-mer seeds (length k)
        → Reject high-frequency seeds (repetitive regions)
        → Retain low-frequency anchors
       ↓
Step 2: Banded Smith-Waterman
        → Define band window around each anchor
        → Compute local alignment score matrix within band
        → Track best score and position
       ↓
Step 3: Genetic Algorithm Initialization
        → Generate N candidate parameter sets
        → Assign random/heuristic values
       ↓
Step 4: Optimization Loop (G generations)
        → Evaluate fitness per candidate
        → Crossover + Mutation
        → Gradient descent refinement
        → Track best fitness Fmax
       ↓
Step 5: Output
        → Best parameter set → final alignment
        → Metrics: accuracy, precision, recall, F1
        → Visualizations: 3D plot, heatmap, helix
```

### Evaluation Metrics

| Metric | Formula |
|---|---|
| Accuracy | (TP + TN) / (TP + TN + FP + FN) |
| Precision | TP / (TP + FP) |
| Recall | TP / (TP + FN) |
| F1-Score | 2 × Precision × Recall / (Precision + Recall) |

### Performance Results (Three Species)

| Species | Method | Accuracy | Precision | Recall | F1 | Runtime (s) |
|---|---|---|---|---|---|---|
| *Homo sapiens* | Improved Aligner | 0.94 | 0.93 | 0.95 | 0.94 | 1.7 |
| | Old HAlign | 0.87 | 0.85 | 0.89 | 0.87 | 0.7 |
| | T-Coffee Sim. | 0.91 | 0.90 | 0.92 | 0.91 | 2.1 |
| *Arabidopsis thaliana* | Improved Aligner | 0.92 | 0.91 | 0.93 | 0.92 | 1.6 |
| | Old HAlign | 0.84 | 0.82 | 0.86 | 0.84 | 0.6 |
| | T-Coffee Sim. | 0.89 | 0.88 | 0.90 | 0.89 | 1.9 |
| *Neurospora crassa* | Improved Aligner | 0.95 | 0.94 | 0.96 | 0.95 | 1.8 |
| | Old HAlign | 0.89 | 0.86 | 0.91 | 0.88 | 0.8 |
| | T-Coffee Sim. | 0.92 | 0.91 | 0.93 | 0.92 | 2.2 |

*Note: Metrics are illustrative values based on simulated experimental conditions.*

---

## Citations

If you use this code or dataset in your research, please cite:

```
Ria K, Rhema Yashaswini O, Meera KN, Mamatha TM, Lin Y. (2024).
Refining DNA Sequence Alignment Accuracy in HAlign through Gradient Descent Techniques.
PeerJ Computer Science. Manuscript 142426.
```

### Key References

1. Wang L, et al. Interactive visualization of sequence alignment data with 3D plots. *BMC Bioinformatics*, 22(1):1–10, 2021.
2. Zhang Y, et al. GPU acceleration of banded Smith-Waterman for efficient long-read alignment. *Computational Biology and Chemistry*, 88:107392, 2021.
3. Guo Z, et al. HAlign 3: Fast multiple alignment of ultra-large numbers of similar DNA/RNA sequences. *Frontiers in Bioinformatics*, 2, 2022.
4. Guo Z, et al. HAlign 4: A new strategy for rapidly aligning millions of sequences. *Bioinformatics*, 40(12):2991–2998, 2024.
5. Smith F, Brown R, Lee S. Deep learning-based sequence alignment. *Briefings in Bioinformatics*, 23(6):bbad136, 2022.
6. Singh SN, et al. Optimization Research on DNA Sequence Alignment Algorithm: HybridAlign approach. *ACM Transactions on Bioinformatics*, 12(2):45–58, 2025.

---

## License & Contribution Guidelines

This project is made available for **academic and research purposes** in accordance with PeerJ's open data and open code policy.

- **Code license:** MIT License — free to use, modify, and distribute with attribution
- **Data:** All sequences are synthetically generated and freely available for reuse
- **Contributions:** Pull requests and issue reports are welcome via the GitHub repository

For questions regarding the methodology or results, contact:  
**K N Meera** — kn_meera@blr.amrita.edu  
Department of Mathematics, Amrita School of Engineering, Amrita Vishwa Vidyapeetham, Bengaluru 560035, India

---

*Submitted to PeerJ Computer Science | Manuscript 142426 | https://peerj.com/manuscripts/142426/edit/*
