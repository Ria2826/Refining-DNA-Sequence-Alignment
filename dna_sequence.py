import numpy as np
import matplotlib.pyplot as plt
import math
import seaborn as sns

# 1. Seed Filter
def find_seeds(ref, query, k, max_freq):
    seed_starts = []
    ref_positions = []
    seen = set()
    for i in range(len(query) - k + 1):
        s = query[i:i+k]
        if s not in seen:
            occ = [j for j in range(len(ref) - k + 1) if ref[j:j+k] == s]
            if occ and len(occ) <= max_freq:
                seed_starts.append(i)
                ref_positions.append(occ[0])
            seen.add(s)
    return seed_starts, ref_positions

# 2. Banded Smith Waterman Alignment
def smith_waterman_banded(query, ref_window):
    match_score = 2
    mismatch_penalty = -1
    gap_penalty = -2
    m = len(query)
    n = len(ref_window)
    score_matrix = np.zeros((m+1, n+1))
    max_score = 0
    peak_i, peak_j = 0, 0
    for i in range(1, m+1):
        for j in range(1, n+1):
            diag = score_matrix[i-1, j-1] + (match_score if query[i-1] == ref_window[j-1] else mismatch_penalty)
            up = score_matrix[i-1, j] + gap_penalty
            left = score_matrix[i, j-1] + gap_penalty
            score_matrix[i, j] = max(0, diag, up, left)
            if score_matrix[i, j] > max_score:
                max_score = score_matrix[i, j]
                peak_i, peak_j = i, j
    return max_score, score_matrix, peak_i, peak_j

# 3. Improved hyperboloid plot
def improved_hyperboloid_score_matrix(score_matrix):
    X = np.arange(score_matrix.shape[1])
    Y = np.arange(score_matrix.shape[0])
    X, Y = np.meshgrid(X, Y)
    Z_orig = score_matrix
    Z_min, Z_max = np.min(Z_orig), np.max(Z_orig)
    Z_norm = (Z_orig - Z_min) / (Z_max - Z_min + 1e-5)
    offset = 1e-2
    X_c = X - np.mean(X)
    Y_c = Y - np.mean(Y)
    hyperboloid_shape = np.sqrt(X_c**2 + Y_c**2 + offset)
    Z = Z_norm * hyperboloid_shape * 10
    return X, Y, Z

def plot_improved_hyperboloid(score_matrix, best_win_start, best_win_end, best_peak_i, best_peak_j, query_len, ref_len):
    X, Y, Z = improved_hyperboloid_score_matrix(score_matrix)
    fig = plt.figure(figsize=(14, 10))
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(X, Y, Z, cmap='viridis', edgecolor='none', alpha=0.8)
    ax.set_xlabel('Reference Window Index')
    ax.set_ylabel('Query Index')
    ax.set_zlabel('Transformed Score Height')
    plt.title('DNA Alignment Score Hyperboloid Structure (Improved for Longer Sequences)')
    fig.colorbar(surf, shrink=0.6, aspect=10)
    ax.text2D(0.07, 0.95, f"Reference Window: {best_win_start} - {best_win_end}", transform=ax.transAxes, color='blue')
    ax.text2D(0.07, 0.91, f"Query Range: 0 - {query_len}", transform=ax.transAxes, color='magenta')
    ax.text2D(0.07, 0.87, f"Alignment peak at: ({best_peak_j}, {best_peak_i})", transform=ax.transAxes, color='green')
    ax.scatter(best_peak_j, best_peak_i, Z[best_peak_i, best_peak_j], color='red', s=80, label='Peak Alignment')
    ax.quiver(X.min(), Y.min(), Z.max(), X.max()-X.min(), 0, 0, color='blue', arrow_length_ratio=0.1, linewidth=2)
    ax.quiver(X.min(), Y.min(), Z.max(), 0, Y.max()-Y.min(), 0, color='magenta', arrow_length_ratio=0.1, linewidth=2)
    ax.legend()
    plt.show()

# 4. DNA double helix text visualization
def print_dna_double_helix(sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    seq = sequence.upper()
    comp = ''.join(complement.get(base, 'N') for base in seq)
    n = len(seq)
    bases_per_turn = 10
    angle_step = 2 * math.pi / bases_per_turn
    print("\nDNA Double Helix Representation (Approximate):\n")
    for i in range(n):
        angle = i * angle_step
        indent = int(5 + 5 * math.sin(angle))
        left_base = seq[i]
        right_base = comp[n - 1 - i]
        connector = '-' * max(1, 6 - indent)
        print(' ' * indent + left_base.upper() + ' ' + connector + ' ' * (10 - indent) + right_base.upper())

# 5. Simulated old HAlign method scoring
def old_halign_score(ref, query):
    seed_length = 8
    max_seed_freq = 10
    seed_starts, ref_positions = find_seeds(ref, query, seed_length, max_seed_freq)
    total_score = 0
    for i, rpos in zip(seed_starts, ref_positions):
        matches = 0
        for offset in range(seed_length):
            if rpos + offset < len(ref) and i + offset < len(query) and ref[rpos + offset] == query[i + offset]:
                matches += 1
        total_score += matches * 2 - (seed_length - matches)
    return total_score

# 6. Simulated T-Coffee scoring with iterative refinement
def t_coffee_simulated_score(ref, query, iterations=5):
    base_score, _, _, _ = smith_waterman_banded(query, ref)
    match_score = 2
    mismatch_penalty = -1
    score = base_score
    for _ in range(iterations):
        mismatch_penalty += 0.1
        temp_score = 0
        for i in range(min(len(query), len(ref))):
            if query[i] == ref[i]:
                temp_score += match_score
            else:
                temp_score += mismatch_penalty
        score = max(score, temp_score)
    return score

# 7. Heatmap comparison plot for metrics
def plot_heatmap_comparison_metrics(metrics_dict):
    methods = list(metrics_dict.keys())
    metrics = ['accuracy', 'precision', 'recall', 'f1']
    data = np.array([[metrics_dict[method][metric] for metric in metrics] for method in methods])
    plt.figure(figsize=(8,6))
    ax = sns.heatmap(data, annot=True, fmt='.2f', cmap='YlGnBu', xticklabels=[m.capitalize() for m in metrics], yticklabels=methods, vmin=0, vmax=1)
    plt.title('Heatmap: DNA Alignment Methods Metric Comparison', fontsize=14)
    plt.xlabel('Metric')
    plt.ylabel('Method')
    plt.tight_layout()
    plt.show()

# 8. Main combined function WITHOUT BLAST/CLUSTALW
def combined_alignment_and_comparison():
    ref = input('Enter reference DNA sequence: ').upper().replace(' ', '').replace('\n','')
    query = input('Enter query DNA sequence: ').upper().replace(' ', '').replace('\n','')
    seed_length = 12
    max_seed_freq = 5
    band_w = 50
    seed_starts, ref_positions = find_seeds(ref, query, seed_length, max_seed_freq)
    if not seed_starts:
        print("No seeds found under given threshold. Consider lowering max_seed_freq or seed_length.")
        return
    best_score = float('-inf')
    best_pos = None
    best_win_start = None
    best_win_end = None
    best_matrix = None
    best_peak_i = None
    best_peak_j = None
    for idx in range(len(ref_positions)):
        ref_start = ref_positions[idx] - seed_starts[idx]
        win_start = max(0, ref_start - band_w)
        win_end = min(len(ref), ref_start + len(query) + band_w)
        ref_window = ref[win_start:win_end]
        score, matrix, peak_i, peak_j = smith_waterman_banded(query, ref_window)
        if score > best_score:
            best_score = score
            best_pos = ref_start
            best_win_start = win_start
            best_win_end = win_end
            best_matrix = matrix
            best_peak_i = peak_i
            best_peak_j = peak_j
    halign_score = old_halign_score(ref, query)
    tcoffee_score = t_coffee_simulated_score(ref, query)
    print(f"\nImproved Aligner Score: {int(best_score)} at position {best_pos}")
    print(f"Simulated Old HAlign Score: {halign_score}")
    print(f"Simulated T-Coffee Score (iterative): {int(tcoffee_score)}")
    plot_improved_hyperboloid(best_matrix, best_win_start, best_win_end, best_peak_i, best_peak_j, len(query), len(ref))
    print_dna_double_helix(query)

    # Mock metrics for demonstration
    metrics = {
        "Improved Aligner": {"accuracy": 0.92, "precision": 0.91, "recall": 0.93, "f1": 0.92},
        "Old HAlign": {"accuracy": 0.85, "precision": 0.83, "recall": 0.87, "f1": 0.85},
        "T-Coffee": {"accuracy": 0.89, "precision": 0.88, "recall": 0.90, "f1": 0.89}
        # "BLAST" and "ClustalW" removed since they are not available
    }
    plot_heatmap_comparison_metrics(metrics)

if __name__ == "__main__":
    combined_alignment_and_comparison()
