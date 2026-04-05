"""
拟时序-Cargo时空对应分析
目标：将外泌体cargo（蛋白/miRNA）映射到半月板退变轨迹的不同阶段
回答：什么蛋白/miRNA在退变的哪个阶段起作用？

轨迹：FC-ECM (起点) → FC-Regulatory/Hypertrophic (中间) → FC-Degenerated (终点)

数据来源:
- protein_cargo_celltype_scores.csv: 各细胞类型的蛋白cargo评分
- miRNA_cargo_celltype_scores.csv: 各细胞类型的miRNA cargo评分
- TF_differential_activity.csv: TF的退变vs正常差异表达

分析策略:
1. 定义拟时序阶段细胞类型映射
2. 计算每个蛋白/miRNA在不同轨迹阶段的活性
3. 识别阶段特异性cargo (early/mid/late phase)
4. 构建cargo-轨迹时间线
5. 预测关键cargo的干预窗口期

作者: WorkBuddy
日期: 2026-04-04
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle
import warnings
warnings.filterwarnings('ignore')

# ==================== 轨迹阶段定义 ====================
# 根据单细胞分析和记忆：FC-ECM → FC-Degenerated 为核心退变轨迹
TRAJECTORY_STAGES = {
    'Early': {
        'cell_types': ['FC-ECM', 'Endothelial', 'Macrophage'],
        'description': '早期：正常ECM维持阶段',
        'color': '#2ca02c'  # Green
    },
    'Mid': {
        'cell_types': ['FC-Regulatory', 'FC-Inflammatory', 'FC-Metabolic',
                      'Stem_Progenitor', 'Meniscus_Progenitor'],
        'description': '中期：炎症响应与代偿阶段',
        'color': '#ff7f0e'  # Orange
    },
    'Late': {
        'cell_types': ['FC-Degenerated', 'Hypertrophic_Chondrocyte',
                      'Secretory_Chondrocyte'],
        'description': '晚期：退变与肥大阶段',
        'color': '#d62728'  # Red
    }
}

# ==================== 数据加载 ====================
DATA_DIR = r'C:\Users\89367\Desktop\CNS半月板\结果处理完\1ALL_RESULTS_FINAL\results\tables'
OUTPUT_DIR = r'c:/Users/89367/WorkBuddy/20260331203952/trajectory_cargo_analysis'

import os
os.makedirs(OUTPUT_DIR, exist_ok=True)

def load_data():
    """加载所有相关数据"""
    protein_scores = pd.read_csv(f'{DATA_DIR}/protein_cargo_celltype_scores.csv', index_col=0)
    mirna_scores = pd.read_csv(f'{DATA_DIR}/miRNA_cargo_celltype_scores.csv', index_col=0)
    tf_diff = pd.read_csv(f'{DATA_DIR}/TF_differential_activity.csv')

    print(f"[OK] Protein cargo shape: {protein_scores.shape}")
    print(f"[OK] miRNA cargo shape: {mirna_scores.shape}")
    print(f"[OK] TF differential shape: {tf_diff.shape}")

    return protein_scores, mirna_scores, tf_diff


# ==================== 轨迹阶段映射函数 ====================
def map_to_stage(cell_type):
    """将细胞类型映射到轨迹阶段"""
    for stage, info in TRAJECTORY_STAGES.items():
        if cell_type in info['cell_types']:
            return stage
    return 'Other'


def calculate_stage_cargo_scores(cargo_scores_df, cargo_type='protein'):
    """
    计算每个cargo在轨迹阶段的平均评分

    Args:
        cargo_scores_df: cargo x 细胞类型 评分矩阵
        cargo_type: 'protein' 或 'mirna'

    Returns:
        DataFrame: cargo x 轨迹阶段 平均评分
    """
    # 为每个阶段计算平均分
    stage_scores = {}

    for stage, info in TRAJECTORY_STAGES.items():
        stage_cells = [ct for ct in info['cell_types'] if ct in cargo_scores_df.columns]

        if stage_cells:
            # 取平均
            stage_scores[stage] = cargo_scores_df[stage_cells].mean(axis=1)
        else:
            stage_scores[stage] = pd.Series(0, index=cargo_scores_df.index)

    stage_df = pd.DataFrame(stage_scores)

    # 归一化 (0-1)
    stage_df = (stage_df - stage_df.min().min()) / (stage_df.max().max() - stage_df.min().min())

    return stage_df


def identify_phase_specific_cargos(stage_df, threshold=0.7):
    """
    识别阶段特异性cargo (某个阶段评分显著高于其他阶段)

    Args:
        stage_df: cargo x 轨迹阶段 评分矩阵
        threshold: 特异性阈值 (某阶段评分 / 其他阶段平均评分)

    Returns:
        dict: {stage: [cargo_list]}
    """
    phase_specific = {'Early': [], 'Mid': [], 'Late': []}

    for cargo in stage_df.index:
        scores = stage_df.loc[cargo]
        max_stage = scores.idxmax()
        max_score = scores[max_stage]

        # 计算其他阶段平均
        other_stages = [s for s in scores.index if s != max_stage]
        avg_other = scores[other_stages].mean()

        # 特异性判断
        if avg_other > 0 and max_score / avg_other > threshold:
            phase_specific[max_stage].append(cargo)

    return phase_specific


# ==================== 时空可视化函数 ====================
def plot_cargo_trajectory_timeline(stage_df, cargo_type, top_n=20):
    """绘制cargo在轨迹阶段的活性时间线"""

    # 选择top cargo (按最大活跃度排序)
    stage_df['max_activity'] = stage_df.max(axis=1)
    top_cargos = stage_df.nlargest(top_n, 'max_activity').drop(columns=['max_activity'])

    # 绘图
    fig, ax = plt.subplots(figsize=(14, 10))

    # 为每个cargo绘制条带
    y_positions = np.arange(len(top_cargos))
    bar_height = 0.6

    for i, (idx, row) in enumerate(top_cargos.iterrows()):
        # Early阶段
        ax.barh(y_positions[i], row['Early'], left=0, height=bar_height,
                color=TRAJECTORY_STAGES['Early']['color'], alpha=0.8, label='Early' if i==0 else '')

        # Mid阶段
        ax.barh(y_positions[i], row['Mid'], left=row['Early'], height=bar_height,
                color=TRAJECTORY_STAGES['Mid']['color'], alpha=0.8, label='Mid' if i==0 else '')

        # Late阶段
        ax.barh(y_positions[i], row['Late'], left=row['Early']+row['Mid'], height=bar_height,
                color=TRAJECTORY_STAGES['Late']['color'], alpha=0.8, label='Late' if i==0 else '')

        # 标注cargo名称
        ax.text(0.02, y_positions[i], idx, va='center', fontsize=9, fontweight='bold')

    # 轨迹阶段标签
    ax.text(0.33, len(top_cargos) + 0.5, 'EARLY\n(ECM Maintenance)',
            ha='center', va='bottom', fontsize=12, fontweight='bold',
            color=TRAJECTORY_STAGES['Early']['color'])
    ax.text(0.66, len(top_cargos) + 0.5, 'MID\n(Inflammatory/Compensatory)',
            ha='center', va='bottom', fontsize=12, fontweight='bold',
            color=TRAJECTORY_STAGES['Mid']['color'])
    ax.text(1.0, len(top_cargos) + 0.5, 'LATE\n(Degeneration/Hypertrophy)',
            ha='center', va='bottom', fontsize=12, fontweight='bold',
            color=TRAJECTORY_STAGES['Late']['color'])

    # 阶段分隔线
    ax.axvline(x=0.33, ymin=0, ymax=0.95, color='black', linestyle='--', alpha=0.3)
    ax.axvline(x=0.66, ymin=0, ymax=0.95, color='black', linestyle='--', alpha=0.3)

    ax.set_yticks([])
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.5, len(top_cargos) + 1)

    ax.set_xlabel('Trajectory Progress', fontsize=12, fontweight='bold')
    ax.set_title(f'{cargo_type.capitalize()} Cargo Activity Across Meniscus Degeneration Trajectory',
                fontsize=14, fontweight='bold', pad=20)

    # 去掉边框
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_linewidth(1.2)
    ax.spines['left'].set_visible(False)

    plt.tight_layout()

    # 保存
    output_file = f'{OUTPUT_DIR}/trajectory_timeline_{cargo_type}.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.savefig(output_file.replace('.png', '.pdf'), bbox_inches='tight')
    plt.close()

    print(f"[OK] {cargo_type} timeline saved to {output_file}")


def plot_phase_specific_heatmap(phase_specific_cargos, cargo_type):
    """绘制阶段特异性cargo热图"""

    if not any(phase_specific_cargos.values()):
        print(f"[WARNING] No phase-specific {cargo_type} found")
        return

    # 收集所有阶段特异性cargo
    all_cargos = []
    for stage, cargos in phase_specific_cargos.items():
        for cargo in cargos:
            all_cargos.append({'cargo': cargo, 'stage': stage})

    df = pd.DataFrame(all_cargos)

    # 构建矩阵
    unique_cargos = df['cargo'].unique()
    stage_matrix = pd.DataFrame(0, index=unique_cargos, columns=['Early', 'Mid', 'Late'])

    for _, row in df.iterrows():
        stage_matrix.loc[row['cargo'], row['stage']] = 1

    # 绘图
    fig, ax = plt.subplots(figsize=(10, max(8, len(unique_cargos) * 0.3)))

    cmap = sns.color_palette(['white', '#2ca02c', '#ff7f0e', '#d62728'], 4)

    sns.heatmap(stage_matrix, cmap=cmap, cbar=False, ax=ax,
                yticklabels=True, xticklabels=True,
                linewidths=0.5, linecolor='gray')

    ax.set_xlabel('Trajectory Stage', fontsize=12, fontweight='bold')
    ax.set_ylabel(f'{cargo_type.capitalize()} Cargo', fontsize=12, fontweight='bold')
    ax.set_title(f'Phase-Specific {cargo_type.capitalize()} Cargos',
                fontsize=14, fontweight='bold', pad=15)

    plt.tight_layout()

    # 保存
    output_file = f'{OUTPUT_DIR}/phase_specific_{cargo_type}.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.savefig(output_file.replace('.png', '.pdf'), bbox_inches='tight')
    plt.close()

    print(f"[OK] Phase-specific {cargo_type} heatmap saved")


# ==================== 干预窗口期预测 ====================
def predict_intervention_window(cargo_name, stage_df, cargo_type='protein'):
    """
    预测特定cargo的最佳干预窗口期

    Args:
        cargo_name: cargo名称
        stage_df: 轨迹阶段评分矩阵
        cargo_type: cargo类型

    Returns:
        dict: 干预窗口期信息
    """
    if cargo_name not in stage_df.index:
        return None

    scores = stage_df.loc[cargo_name]
    max_stage = scores.idxmax()
    max_score = scores[max_stage]

    # 根据峰值阶段推荐干预窗口
    intervention_windows = {
        'Early': '预防期：在退变早期补充，维持ECM稳态',
        'Mid': '治疗期：在炎症代偿期补充，控制炎症进展',
        'Late': '修复期：在退变晚期补充，促进组织再生'
    }

    return {
        'cargo': cargo_name,
        'peak_stage': max_stage,
        'peak_score': max_score,
        'intervention_window': intervention_windows.get(max_stage, 'Unknown'),
        'scores': scores.to_dict()
    }


# ==================== 主分析流程 ====================
def main():
    """执行完整分析流程"""

    print("=" * 70)
    print("拟时序-Cargo时空对应分析")
    print("=" * 70)

    # 1. 加载数据
    print("\n[STEP 1] Loading data...")
    protein_scores, mirna_scores, tf_diff = load_data()

    # 2. 计算轨迹阶段评分
    print("\n[STEP 2] Calculating trajectory stage scores...")

    # Protein cargo
    protein_stage_df = calculate_stage_cargo_scores(protein_scores, 'protein')
    print(f"[OK] Protein stage scores calculated")

    # miRNA cargo
    mirna_stage_df = calculate_stage_cargo_scores(mirna_scores, 'mirna')
    print(f"[OK] miRNA stage scores calculated")

    # 3. 识别阶段特异性cargo
    print("\n[STEP 3] Identifying phase-specific cargos...")

    protein_phase_specific = identify_phase_specific_cargos(protein_stage_df, threshold=1.5)
    print(f"[OK] Protein phase-specific cargos:")
    for stage, cargos in protein_phase_specific.items():
        print(f"   {stage}: {len(cargos)} cargos")
        if cargos:
            print(f"   Examples: {', '.join(cargos[:5])}")

    mirna_phase_specific = identify_phase_specific_cargos(mirna_stage_df, threshold=1.5)
    print(f"[OK] miRNA phase-specific cargos:")
    for stage, cargos in mirna_phase_specific.items():
        print(f"   {stage}: {len(cargos)} cargos")
        if cargos:
            print(f"   Examples: {', '.join(cargos[:5])}")

    # 4. 可视化
    print("\n[STEP 4] Generating visualizations...")

    # Protein时间线
    plot_cargo_trajectory_timeline(protein_stage_df, 'protein', top_n=25)

    # miRNA时间线
    plot_cargo_trajectory_timeline(mirna_stage_df, 'mirna', top_n=20)

    # 阶段特异性热图
    plot_phase_specific_heatmap(protein_phase_specific, 'protein')
    plot_phase_specific_heatmap(mirna_phase_specific, 'mirna')

    # 5. 干预窗口期预测示例
    print("\n[STEP 5] Predicting intervention windows for key cargos...")

    # 选择关键蛋白
    key_proteins = ['FN1', 'COL1A1', 'TGFB1', 'TGFB3', 'VEGFA', 'FGF2']

    for protein in key_proteins:
        window_info = predict_intervention_window(protein, protein_stage_df, 'protein')
        if window_info:
            print(f"\n{protein}:")
            print(f"  峰值阶段: {window_info['peak_stage']}")
            print(f"  峰值评分: {window_info['peak_score']:.3f}")
            print(f"  干预窗口: {window_info['intervention_window']}")

    # 选择关键miRNA
    key_mirnas = ['hsa-miR-140-5p', 'hsa-miR-29a-3p', 'hsa-miR-146a-5p']

    for mirna in key_mirnas:
        window_info = predict_intervention_window(mirna, mirna_stage_df, 'mirna')
        if window_info:
            print(f"\n{mirna}:")
            print(f"  峰值阶段: {window_info['peak_stage']}")
            print(f"  峰值评分: {window_info['peak_score']:.3f}")
            print(f"  干预窗口: {window_info['intervention_window']}")

    # 6. 保存结果
    print("\n[STEP 6] Saving results...")

    # 保存轨迹阶段评分
    protein_stage_df.to_csv(f'{OUTPUT_DIR}/protein_trajectory_stage_scores.csv')
    mirna_stage_df.to_csv(f'{OUTPUT_DIR}/mirna_trajectory_stage_scores.csv')

    # 保存阶段特异性cargo (统一长度)
    max_len = max(len(v) for v in protein_phase_specific.values())
    protein_phase_data = {k: v + [''] * (max_len - len(v)) for k, v in protein_phase_specific.items()}
    protein_phase_df = pd.DataFrame(protein_phase_data)
    protein_phase_df.to_csv(f'{OUTPUT_DIR}/protein_phase_specific_cargos.csv', index=False)

    max_len = max(len(v) for v in mirna_phase_specific.values())
    mirna_phase_data = {k: v + [''] * (max_len - len(v)) for k, v in mirna_phase_specific.items()}
    mirna_phase_df = pd.DataFrame(mirna_phase_data)
    mirna_phase_df.to_csv(f'{OUTPUT_DIR}/mirna_phase_specific_cargos.csv', index=False)

    # 生成综合报告
    report_lines = [
        "# 拟时序-Cargo时空对应分析报告",
        "",
        "## 分析概述",
        "",
        "本分析将外泌体cargo（蛋白和miRNA）映射到半月板退变轨迹的不同阶段，",
        "识别阶段特异性cargo并预测最佳干预窗口期。",
        "",
        "## 轨迹阶段定义",
        "",
        "- **Early (早期)**: 正常ECM维持阶段",
        f"  - 细胞类型: {', '.join(TRAJECTORY_STAGES['Early']['cell_types'])}",
        "",
        "- **Mid (中期)**: 炎症响应与代偿阶段",
        f"  - 细胞类型: {', '.join(TRAJECTORY_STAGES['Mid']['cell_types'])}",
        "",
        "- **Late (晚期)**: 退变与肥大阶段",
        f"  - 细胞类型: {', '.join(TRAJECTORY_STAGES['Late']['cell_types'])}",
        "",
        "## 关键发现",
        "",
        "### 阶段特异性Protein Cargo",
        ""
    ]

    for stage, cargos in protein_phase_specific.items():
        report_lines.append(f"#### {stage}阶段 ({len(cargos)}个)")
        if cargos:
            report_lines.append(f"- 特异性蛋白: {', '.join(cargos[:10])}")
        report_lines.append("")

    report_lines.extend([
        "### 阶段特异性miRNA Cargo",
        ""
    ])

    for stage, cargos in mirna_phase_specific.items():
        report_lines.append(f"#### {stage}阶段 ({len(cargos)}个)")
        if cargos:
            report_lines.append(f"- 特异性miRNA: {', '.join(cargos[:10])}")
        report_lines.append("")

    report_lines.extend([
        "### 干预窗口期推荐",
        "",
        "根据轨迹分析，不同cargo的最佳干预窗口期如下：",
        ""
    ])

    # 添加关键cargo的干预窗口
    for protein in key_proteins:
        window_info = predict_intervention_window(protein, protein_stage_df, 'protein')
        if window_info:
            report_lines.append(f"#### {protein}")
            report_lines.append(f"- 峰值阶段: {window_info['peak_stage']}")
            report_lines.append(f"- 干预窗口: {window_info['intervention_window']}")
            report_lines.append("")

    # 保存报告
    with open(f'{OUTPUT_DIR}/trajectory_cargo_analysis_report.md', 'w', encoding='utf-8') as f:
        f.write('\n'.join(report_lines))

    print(f"\n[OK] All results saved to: {OUTPUT_DIR}")
    print(f"[OK] Analysis report: {OUTPUT_DIR}/trajectory_cargo_analysis_report.md")

    print("\n" + "=" * 70)
    print("Analysis completed successfully!")
    print("=" * 70)


if __name__ == '__main__':
    main()
