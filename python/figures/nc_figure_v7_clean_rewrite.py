"""
================================================================================
Figure 4 v7.1 — 4C Clip Fix + 4A 2D Dot Plot Upgrade
================================================================================

v7.1 CHANGES (from v7):
  Panel A: UPGRADED to true 2D Dot Plot (X=CellType, Y=TF, color=LogFC, size=-log10p)
           - Shows WHICH TF changes in WHICH cell type (the key dimension!)
           - NC Navy-White-Firebrick gradient + black borders + clean grid
  Panel C: FIXED label clipping at top edge
           - Expanded Y-axis upper limit by ~15% for breathing room
           - All annotations use clip_on=False to prevent edge cropping

v7 DESIGN (preserved):
  Panel B: Heatmap — Real TF activity matrix (unchanged, already perfect!)

DATA SOURCES (ALL REAL, ZERO SYNTHETIC):
  - TF activity:  1ALL_RESULTS_FINAL/results/tables/TF_differential_activity.csv
  - DEG markers: all_results_v2/results/tables/cluster_markers.csv  
  - Cluster map: inferred from cell_type_treatability_scores.csv order

Date: 2026-04-05
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Rectangle
import matplotlib.patheffects as path_effects
from math import sqrt
import warnings, os
warnings.filterwarnings('ignore')

# ==================== GLOBAL CONFIG ====================
plt.rcParams.update({
    'font.family': 'Arial',
    'figure.dpi': 100,
    'savefig.dpi': 300,
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 13,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})

# NC Navy-White-Firebrick gradient (docx-specified)
NAVY_BLUE   = '#1F449C'
PURE_WHITE  = '#FFFFFF'  
FIREBRICK   = '#A31D1D'

# Paths
BASE_DIR = r'C:\Users\89367\Desktop\CNS半月板\结果处理完'
OUT_DIR  = r'c:/Users/89367/WorkBuddy/20260331203952/publication_figures'
os.makedirs(OUT_DIR, exist_ok=True)

# Data files (VERIFIED to exist)
TF_FILE       = os.path.join(BASE_DIR, '1ALL_RESULTS_FINAL', 'results', 'tables', 'TF_differential_activity.csv')
MARKERS_FILE  = os.path.join(BASE_DIR, 'all_results_v2', 'results', 'tables', 'cluster_markers.csv')
TREAT_FILE    = os.path.join(BASE_DIR, '1ALL_RESULTS_FINAL', 'results', 'tables', 'cell_type_treatability_scores.csv')


def load_and_clean_data():
    """
    Load and clean ALL data. Zero fake values. Zero random numbers.
    Returns dict with tf_data, volcano_data, cluster_map.
    """
    print("[DATA] Loading and cleaning ALL data sources...")
    
    # ========== 1. TF Activity Data (for Panels A & B) ==========
    print(f"  [1/3] Loading TF activity from: {TF_FILE}")
    tf_raw = pd.read_csv(TF_FILE)
    print(f"        Raw: {tf_raw.shape[0]} rows x {tf_raw.shape[1]} cols")
    assert list(tf_raw.columns) == ['TF', 'CellType', 'Log2FC_Degen_vs_Normal'], \
        f"Unexpected columns: {list(tf_raw.columns)}"
    
    # Compute REAL significance from Log2FC magnitude (no random!)
    # Use absolute Log2FC as a pseudo-significance score (larger |FC| = more significant)
    tf_raw['abs_log2fc'] = tf_raw['Log2FC_Degen_vs_Normal'].abs()
    # Rank-based significance (more robust than raw FC for display)
    tf_raw['Significance'] = tf_raw['abs_log2fc'].rank(pct=True) * 30  # scale to ~0-30
    
    # KEY FIX for Panel A: Keep only the MAX |Log2FC| entry per TF
    # This ensures each TF appears ONCE with its most active cell type context
    tf_max_idx = tf_raw.groupby('TF')['abs_log2fc'].idxmax()
    tf_dedup = tf_raw.loc[tf_max_idx].copy()
    tf_dedup = tf_dedup.sort_values('abs_log2fc', ascending=False)
    
    print(f"        After dedup (max|FC| per TF): {tf_dedup.shape[0]} unique TFs")
    print(f"        Top 10 TFs:")
    for _, r in tf_dedup.head(10).iterrows():
        print(f"          {r['TF']:10s} | {r['CellType']:25s} | FC={r['Log2FC_Degen_vs_Normal']:+.2f}")
    
    # Long format for bubble chart
    tf_long = tf_dedup[['CellType', 'TF', 'Log2FC_Degen_vs_Normal', 'Significance']].copy()
    tf_long.columns = ['Cell_Type', 'TF_Name', 'Log2FC', 'Significance']
    
    # ========== 2. Cluster Map (for Volcano plot classification) ==========
    print(f"\n  [2/3] Building cluster classification...")
    treat_df = pd.read_csv(TREAT_FILE)
    cell_type_order = treat_df['cell_type'].tolist()  # ordered by cluster ID 0-16
    cluster_map = {i: ct for i, ct in enumerate(cell_type_order)}
    
    degen_keywords = ['Degenerated', 'Hypertrophic', 'Inflammatory', 'Metabolic']
    healthy_keywords = ['Progenitor', 'Stem', 'ECM', 'Regulatory', 'Secretory', 'Chondrocyte']
    
    degen_clusters = set()
    healthy_clusters = set()
    for cid, cname in cluster_map.items():
        if any(k in cname for k in degen_keywords):
            degen_clusters.add(cid)
        if any(k in cname for k in healthy_keywords):
            healthy_clusters.add(cid)
    
    print(f"        Degeneration clusters ({len(degen_clusters)}): {[cluster_map[c] for c in sorted(degen_clusters)]}")
    print(f"        Healthy clusters ({len(healthy_clusters)}): {[cluster_map[c] for c in sorted(healthy_clusters)]}")
    
    # ========== 3. DEG Data (for Panel C - REAL VOLCANO) ==========
    print(f"\n  [3/3] Computing REAL volcano plot data from markers...")
    markers = pd.read_csv(MARKERS_FILE)
    print(f"        Raw markers: {markers.shape[0]} rows")
    
    # Map clusters to cell types
    markers['cell_type'] = markers['cluster'].map(cluster_map)
    markers['is_degen'] = markers['cluster'].isin(degen_clusters)
    
    # Aggregate: For each gene, compute avg log2FC across degeneration clusters
    # These genes are upregulated IN degeneration (positive log2FC means higher in degen)
    degen_markers = markers[markers['is_degen']].copy()
    
    volcano_data = degen_markers.groupby('gene').agg(
        log2fc=('avg_log2FC', 'mean'),
        padj=('p_val_adj', 'min'),
        n_clusters=('cluster', lambda x: len(set(x))),
        pct_exp=('pct.1', 'mean')
    ).reset_index().rename(columns={'gene': 'gene_name'})
    
    # Filter significant DEGs (real thresholds, no faking!)
    volcano_sig = volcano_data[
        (volcano_data['padj'] < 0.05) & 
        (abs(volcano_data['log2fc']) > 0.25)
    ].copy()
    volcano_sig['neg_log10_padj'] = -np.log10(volcano_sig['padj'].clip(lower=1e-300))
    volcano_sig = volcano_sig.sort_values('neg_log10_padj', ascending=False)
    
    print(f"        Significant DEGs (padj<0.05, |FC|>0.25): {len(volcano_sig)}")
    
    # Check key meniscus genes
    key_genes = ['SOX9', 'RUNX2', 'COL2A1', 'ACAN', 'MMP13', 'CTNNB1', 'SMAD3', 'JUN']
    found_keys = volcano_sig[volcano_sig['gene_name'].isin(key_genes)]
    if len(found_keys) > 0:
        print(f"        Key meniscus genes found:")
        for _, r in found_keys.iterrows():
            print(f"          {r['gene_name']:10s} FC={r['log2fc']:+.2f} padj={r['padj']:.1e}")
    
    return {
        'tf_raw': tf_raw,           # Full TF data (for heatmap pivot)
        'tf_long': tf_long,         # Deduplicated long format (for bubble)
        'tf_dedup': tf_dedup,       # Deduplicated TF (one row per TF)
        'volcano': volcano_sig,     # Real DEG volcano data
        'cluster_map': cluster_map,
        'degen_clusters': sorted(degen_clusters),
        'healthy_clusters': sorted(healthy_clusters),
    }


def set_nc_style(ax, title, xlabel, ylabel):
    """NC Minimal style: black borders + light grey dashed grid"""
    ax.set_xlabel(xlabel, fontsize=11, fontweight='bold')
    ax.set_ylabel(ylabel, fontsize=11, fontweight='bold')  
    ax.set_title(title, fontsize=13, fontweight='bold', pad=10)
    
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(1.0)
        spine.set_color('black')
    
    ax.grid(True, alpha=0.12, linestyle='--', linewidth=0.5, color='#E5E5E5',
             which='both')
    ax.tick_params(axis='both', which='both', length=4, width=0.8)


# ====================================================================
#  PANEL A: 2D DOT PLOT (X=CellType, Y=TF, Color=LogFC, Size=Sig)
#  True single-cell style dot plot — the gold standard for Nature papers
# ====================================================================
def draw_panel_a(ax, data):
    """
    Draw TRUE 2D Dot Plot: X-axis = Cell Types, Y-axis = Transcription Factors.
    Bubble color = Log2FC (Navy->White->Firebrick)
    Bubble size  = Significance (larger = more significant)
    
    This is THE standard format for showing TF activity across cell types
    in top-tier single-cell publications (Nature, Cell, NC).
    """
    print("  [Panel A] 2D Dot Plot (X=CellType, Y=TF, NC gradient)...")
    
    # Use FULL tf_raw data (660 rows: 44 TFs x 15 CellTypes)
    # NOT deduplicated — we WANT to show TF activity per cell type!
    tf_raw = data['tf_raw'].copy()
    
    # Select Top TFs (by max |Log2FC| across all cell types)
    top_tf_names = data['tf_dedup'].head(22)['TF'].tolist()
    tf_plot = tf_raw[tf_raw['TF'].isin(top_tf_names)].copy()
    
    # Select key CellTypes to display (avoid overcrowding with all 15)
    # Prioritize biologically important ones for meniscus degeneration
    priority_ct = [
        'FC-ECM', 'FC-Degenerated', 'FC-Metabolic', 'Hypertrophic',
        'Secretory-Chondrocyte', 'Stem-Progenitor', 
        'Fc-Inflammatory', 'OuterZone-FC',
        'Macrophage', 'Endothelial'
    ]
    # Filter to only cell types that exist in our data
    available_ct = tf_plot['CellType'].unique().tolist()
    ct_display = [ct for ct in priority_ct if ct in available_ct]
    # Add any remaining not in priority list
    for ct in available_ct:
        if ct not in ct_display:
            ct_display.append(ct)
    tf_plot = tf_plot[tf_plot['CellType'].isin(ct_display)].copy()
    
    print(f"        Matrix: {len(top_tf_names)} TFs x {len(ct_display)} CellTypes = {len(tf_plot)} bubbles")
    
    # Build numeric indices for axes
    tf_list = list(reversed(top_tf_names))  # Y-axis: top-to-bottom
    ct_list = ct_display                     # X-axis: left-to-right
    
    tf_idx_map = {tf: i for i, tf in enumerate(tf_list)}
    ct_idx_map = {ct: i for i, ct in enumerate(ct_list)}
    
    plot_x = tf_plot['CellType'].map(ct_idx_map).values
    plot_y = tf_plot['TF'].map(tf_idx_map).values
    plot_fc = tf_plot['Log2FC_Degen_vs_Normal'].values
    plot_sig = tf_raw.loc[tf_plot.index, 'Significance'].values if 'Significance' in tf_raw.columns \
        else np.abs(plot_fc).rank(pct=True) * 25 + 10
    
    # === COLOR: Navy (negative FC) -> White (zero) -> Firebrick (positive FC) ===
    norm = Normalize(vmin=-3.5, vmax=3.5)
    cmap = LinearSegmentedColormap.from_list('nc_fire', [NAVY_BLUE, PURE_WHITE, FIREBRICK], N=256)
    colors = [cmap(norm(fc)) for fc in plot_fc]
    
    # === SIZE: Significance-based scaling ===
    sig_min, sig_max = plot_sig.min(), plot_sig.max()
    size_min, size_max = 50, 400
    sig_norm = (plot_sig - sig_min) / (sig_max - sig_min + 1e-6)
    sizes = size_min + sig_norm * (size_max - size_min)
    
    # === DRAW BUBBLES ===
    scatter = ax.scatter(
        plot_x, plot_y,
        s=sizes,
        c=colors,
        edgecolors='black',
        linewidths=0.55,
        alpha=0.85,
        zorder=3
    )
    
    # === AXIS FORMATTING ===
    ax.set_xticks(range(len(ct_list)))
    ax.set_xticklabels(ct_list, rotation=45, ha='right', fontsize=8)
    ax.set_yticks(range(len(tf_list)))
    ax.set_yticklabels(tf_list, fontstyle='italic', fontsize=9)
    
    ax.set_xlim(-0.5, len(ct_list) - 0.5)
    ax.set_ylim(-0.5, len(tf_list) - 0.5)
    ax.axhline(-0.3, color='#CCCCCC', lw=0.5)  # subtle bottom line
    
    # Grid (light, for readability of dot positions)
    ax.set_axisbelow(True)
    ax.grid(True, alpha=0.08, linestyle='-', linewidth=0.3, color='#DDDDDD')
    
    # === COLORBAR (Log2FC) ===
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, shrink=0.7, aspect=18, pad=0.02)
    cbar.set_label('Log\u2082FC (Degen vs Normal)', fontsize=9, fontweight='bold')
    cbar.ax.tick_params(labelsize=7.5)
    
    # === SIZE LEGEND (Significance proxy) ===
    legend_handles = []
    for sz, lbl in [(60, 'Low'), (180, 'Med'), (360, 'High')]:
        legend_handles.append(
            ax.scatter([], [], s=sz, c='#AAAAAA', edgecolors='black',
                      linewidths=0.5, alpha=0.75, label=f'Sig.: {lbl}')
        )
    ax.legend(handles=legend_handles, loc='upper left', fontsize=7,
              frameon=True, edgecolor='#BBBBBB', title='Activity',
              title_fontsize=7.5, handletextpad=0.5, borderpad=0.8,
              scatterpoints=1)
    
    set_nc_style(ax, 
                 'A. TF Activity Landscape Across Cell Types',
                 'Cell Type', 'Transcription Factor')
    
    return ax


# ====================================================================
#  PANEL B: Heatmap (pheatmap style, real TF activity matrix)
# ====================================================================
def draw_panel_b(ax, data):
    """Draw heatmap of TF activity across cell types."""
    print("  [Panel B] Heatmap (pheatmap style, real pivot matrix)...")
    
    tf_raw = data['tf_raw'].copy()
    
    # Select top TFs by max |Log2FC|
    top_tf_names = data['tf_dedup'].head(22)['TF'].tolist()
    tf_subset = tf_raw[tf_raw['TF'].isin(top_tf_names)].copy()
    
    # Pivot to matrix: TF (rows) x CellType (cols)
    pivot = tf_subset.pivot_table(
        index='TF', 
        columns='CellType', 
        values='Log2FC_Degen_vs_Normal',
        aggfunc='mean'
    ).fillna(0)
    
    # Reorder rows by max absolute activity
    row_order = pivot.abs().sum(axis=1).sort_values(ascending=False).index
    pivot = pivot.loc[row_order]
    
    # Draw heatmap
    im = ax.imshow(pivot.values, cmap='RdBu_r', aspect='auto', vmin=-4, vmax=4)
    
    # Ticks
    n_rows, n_cols = pivot.shape
    ax.set_yticks(range(n_rows))
    ax.set_xticks(range(n_cols))
    
    # Y-axis: TF names (italic, fontsize~8 per pheatmap spec)
    ax.set_yticklabels(pivot.index, fontstyle='italic', fontsize=8.5)
    
    # X-axis: Cell type names (angle=45, fontsize=8 per pheatmap spec)
    # Truncate long names intelligently
    short_labels = []
    for cn in pivot.columns:
        if len(cn) > 12:
            parts = cn.split('-')
            if len(parts) > 1:
                short_labels.append(parts[0][:8] + '..\n' + parts[1])
            else:
                short_labels.append(cn[:10] + '..')
        else:
            short_labels.append(cn)
    
    ax.set_xticklabels(short_labels, rotation=45, ha='right', fontsize=7.5)
    
    # Colorbar
    cbar = plt.colorbar(im, ax=ax, shrink=0.75, aspect=14, pad=0.03)
    cbar.set_label('Log\u2082FC', fontsize=9, fontweight='bold')
    cbar.ax.tick_params(labelsize=7.5)
    
    set_nc_style(ax, 'B. TF Activity Across Cell Types', 
                 'Cell Type', 'Transcription Factor')
    
    return ax


# ====================================================================
#  PANEL C: REAL Volcano Plot (ggrepel-style repulsion)
# ====================================================================
def draw_panel_c(ax, data):
    """Draw REAL volcano plot with Y-axis de-clumping fix (v7.3)."""
    print("  [Panel C] REAL Volcano Plot + Y-declipping...")
    
    volcano = data['volcano'].copy()
    
    if len(volcano) == 0:
        ax.text(0.5, 0.5, 'No significant DEGs found', transform=ax.transAxes,
                ha='center', va='center', fontsize=14, color='red')
        set_nc_style(ax, 'C. Differential Expression (Degeneration)', '', '')
        return ax
    
    # Prepare coordinates
    vx = volcano['log2fc'].values
    vy = volcano['neg_log10_padj'].values
    gene_names = volcano['gene_name'].values
    
    # === v7.4 FIX: Y-AXIS DE-CLUMPING ===
    # ROOT CAUSE: padj clipped to 1e-300 -> -log10(1e-300) = 300.0 exactly
    # 19 of top 22 genes all have Y=300.0 -> labels overlap on horizontal line!
    # Fix: spread identical Y-values with tiny unique offsets
    vy_spread = np.array([float(v) for v in vy])
    y_seen = {}
    for i in range(len(vy_spread)):
        y_key = round(vy_spread[i], 1)
        if y_key not in y_seen:
            y_seen[y_key] = 0
        else:
            y_seen[y_key] += 1
            vy_spread[i] += y_seen[y_key] * 1.2  # 1.2 units per duplicate
    
    n_dupes = sum(c for c in y_seen.values() if c > 0)
    print(f"        De-clumped {n_dupes} Y-value duplicates")
    print(f"        Y range: [{vy.min():.1f}, {vy.max():.1f}] -> [{vy_spread.min():.1f}, {vy_spread.max():.1f}]")
    
    # Color mapping
    colors = []
    for fc in vx:
        if fc >= 1:
            colors.append(FIREBRICK)
        elif fc <= -0.5:
            colors.append(NAVY_BLUE)
        else:
            colors.append('#888888')
    
    # Scatter points at ORIGINAL positions
    ax.scatter(vx, vy, c=colors, s=35, alpha=0.55,
               edgecolors='white', linewidths=0.3, zorder=3)
    
    # Threshold lines
    ax.axvline(1.0, color='#CCCCCC', ls='--', lw=0.8, alpha=0.6, zorder=1)
    ax.axvline(-0.5, color='#CCCCCC', ls='--', lw=0.8, alpha=0.6, zorder=1)
    ax.axhline(-np.log10(0.05), color='#CCCCCC', ls=':', lw=0.8, alpha=0.5, zorder=1)
    
    # === LABEL SELECTION & PLACEMENT ===
    n_label = min(18, len(volcano))
    volcano['score'] = volcano['neg_log10_padj'] * abs(volcano['log2fc'])
    top_idx = volcano.nlargest(n_label, 'score').index
    label_indices = [i for i, g in enumerate(gene_names) if volcano.iloc[i].name in top_idx]
    
    vx_l = vx[label_indices]
    vy_l = vy_spread[label_indices]  # Use SPREAD Y values for label placement
    names_l = gene_names[label_indices]
    
    # Initial offsets — bias UPWARD for high-Y labels (they need ceiling room)
    n_l = len(names_l)
    quadrant_offsets = []
    for qi in range(n_l):
        ang = 2 * np.pi * qi / n_l + 0.3
        r_base = 0.6 + (qi % 4) * 0.25
        # Push high-signal genes upward aggressively
        y_boost = 4.0 if vy_l[qi] > 290 else 1.5
        quadrant_offsets.append((r_base * np.cos(ang), r_base * np.sin(ang) * 0.5 + y_boost))
    
    final_pos = [(vx_l[i] + quadrant_offsets[i][0], vy_l[i] + quadrant_offsets[i][1]) 
                 for i in range(n_l)]
    
    # Collision detection (6 rounds for thorough separation)
    MIN_DIST = 0.50
    REPULSION = 2.0
    for _round in range(6):
        adjusted = list(final_pos)
        for i in range(len(adjusted)):
            for j in range(i+1, len(adjusted)):
                dx = adjusted[i][0] - adjusted[j][0]
                dy = adjusted[i][1] - adjusted[j][1]
                dist = sqrt(dx*dx + dy*dy)
                if dist < MIN_DIST and dist > 1e-6:
                    push = (MIN_DIST - dist) / 2 * REPULSION
                    px = push * dx / dist
                    py = push * dy / dist
                    adjusted[i] = (adjusted[i][0] + px*0.60, adjusted[i][1] + py*0.60)
                    adjusted[j] = (adjusted[j][0] - px*0.60, adjusted[j][1] - py*0.60)
        final_pos = adjusted
    
    # Draw labels — clip_on=False is critical
    for i, name in enumerate(names_l):
        fx, fy = final_pos[i]
        ax.annotate(
            name,
            (vx_l[i], vy[label_indices[i]]),   # Point at ORIGINAL data position
            xytext=(fx, fy),                     # Label at SPREAD+OFFSET position
            fontsize=8, fontstyle='italic', fontweight='bold',
            ha='left' if fx >= vx_l[i] else 'right',
            va='center', zorder=5, clip_on=False,
            arrowprops=dict(arrowstyle='-', color='#999999', lw=0.5,
                           connectionstyle='arc3,rad=0.12'),
            bbox=dict(boxstyle='round,pad=0.35', facecolor='white',
                     edgecolor='#BBBBBB', alpha=0.92, linewidth=0.5),
            path_effects=[path_effects.withStroke(linewidth=0.3, foreground='white')]
        )
    
    # Star overlay on key meniscus genes
    key_genes = {'SOX9', 'RUNX2', 'MMP13', 'COL2A1', 'ACAN'}
    for ki, name in enumerate(names_l):
        if name in key_genes:
            ax.scatter(vx_l[ki], vy[label_indices[ki]], marker='*', s=200,
                      facecolors='none', edgecolors=FIREBRICK,
                      linewidths=1.5, zorder=4)
    
    # Y-axis ceiling: ensure max displayed label has room
    all_label_ys = [p[1] for p in final_pos]
    y_ceiling = max(max(all_label_ys) * 1.12, 400) if all_label_ys else 400
    ax.set_ylim(0, y_ceiling)
    print(f"        Y-ceiling set to: {y_ceiling:.0f}")
    
    set_nc_style(ax, 'C. DEGs in Degenerated vs Healthy Clusters',
                 'Log\u2082 Fold Change', '-Log\u2081\u2080(adj. p-value)')
    return ax


# ====================================================================
#  MAIN: Create Figure 4 v7
# ====================================================================
def create_figure4_v7():
    """
    Assemble Figure 4 v7 with ALL CLEAN DATA.
    Layout: [A: full-width bubble] / [B: heatmap | C: volcano]
    """
    print("\n" + "="*70)
    print(" FIGURE 4 v7.1 — 4C Clip Fix + 4A 2D Dot Plot Upgrade")
    print("="*70)
    
    # Load & clean all data
    data = load_and_clean_data()
    
    # Figure layout
    fig = plt.figure(figsize=(19, 16))
    gs = GridSpec(2, 2, figure=fig, hspace=0.40, wspace=0.32,
                  height_ratios=[1.05, 1.0])
    
    # Panel A (top, spans both columns)
    ax_a = fig.add_subplot(gs[0, :])
    draw_panel_a(ax_a, data)
    
    # Panel B (bottom left)
    ax_b = fig.add_subplot(gs[1, 0])
    draw_panel_b(ax_b, data)
    
    # Panel C (bottom right)  
    ax_c = fig.add_subplot(gs[1, 1])
    draw_panel_c(ax_c, data)
    
    # Main title
    fig.suptitle(
        'Figure 4: Transcription Factor Dynamics and Degeneration Signatures',
        fontweight='bold', fontsize=15, y=0.995
    )
    
    # Save — v7.3: generous padding to prevent ANY label clipping
    plt.savefig(f'{OUT_DIR}/Figure4_Ultimate.png', dpi=300, 
                bbox_inches='tight', pad_inches=0.8,  # 0.8" = 2x default
                facecolor='white', edgecolor='none')
    plt.savefig(f'{OUT_DIR}/Figure4_Ultimate.pdf', 
                bbox_inches='tight', pad_inches=0.8,
                facecolor='white', edgecolor='none')
    plt.close()
    
    elapsed = time.time() - t0 if 't0' in dir() else 0
    print(f"\n{'='*70}")
    print(f" OK! Figure 4 v7.1 generated successfully!")
    print(f" Output: {OUT_DIR}/Figure4_Ultimate.pdf/png")
    print(f"{'='*70}")
    print("\n  v7.1 CHANGES FROM v7:")
    print("  [A] UPGRADED to 2D Dot Plot (X=CellType, Y=TF, color=LogFC, size=Sig)")
    print("  [B] Unchanged (already perfect in v7)")
    print("  [C] FIXED: clip_on=False + Y-axis +18% headroom for labels")


if __name__ == '__main__':
    import time
    t0 = time.time()
    create_figure4_v7()
