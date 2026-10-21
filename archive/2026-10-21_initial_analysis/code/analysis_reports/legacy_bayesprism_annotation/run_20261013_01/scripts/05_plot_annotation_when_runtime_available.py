# 本次图源的Python绘图入口；当前缺matplotlib，未运行渲染。
# 在已有可用绘图环境运行，输出仅限本次figures；不安装软件。
import os
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
from pathlib import Path
import subprocess
import sys
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

# Step1. 读入已保存图源与稿件合同
out = Path(__file__).resolve().parents[1]
figure_dir = out / 'figures'
qa = Path('C:/Users/马泽凯/.codex/skills/nature-figure/scripts')
sys.path.insert(0, str(qa))
from audit_panel_alignment import require_matplotlib_panel_alignment
mpl.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'], 'font.size': 7, 'pdf.fonttype': 42, 'svg.fonttype': 'none', 'axes.spines.top': False, 'axes.spines.right': False})
master = pd.read_csv(out / 'CELL_STATE_ANNOTATION_MASTER.tsv', sep='\t')
fdr = pd.read_csv(figure_dir / 'FUNCTIONAL_HEATMAP_GLOBAL_FDR_SOURCE.tsv', sep='\t', index_col=0)
composition = pd.read_csv(figure_dir / 'ECOTYPE_MEMBER_CONFIDENCE_SOURCE.tsv', sep='\t', index_col=0)
markers = pd.read_csv(figure_dir / 'MARKER_DETECTION_HEATMAP_SOURCE.tsv', sep='\t')
plots = []

# Step2. 所有状态和Hallmark保留，颜色是全局FDR而非通路活性
values = -np.log10(fdr.astype(float).clip(lower=1e-12))
fig, ax = plt.subplots(figsize=(12, 12))
mesh = ax.imshow(values.to_numpy(), cmap='Blues', aspect='auto', vmin=0, vmax=4)
ax.set_yticks(range(len(values)), values.index, fontsize=7)
ax.set_xticks(range(len(values.columns)), [c.removeprefix('HALLMARK_') for c in values.columns], rotation=90, rotation_mode='anchor', ha='right', fontsize=7)
ax.set_title('Conditional Hallmark ORA in frozen NMF feature space', fontsize=9)
fig.colorbar(mesh, ax=ax, fraction=.02, pad=.02, label='-log10(global BH FDR), capped at 4')
fig.subplots_adjust(left=.24, bottom=.31, right=.94, top=.96)
plots.append((fig, 'Functional_ORA_global_FDR'))

# Step3. 生态型构成表示冻结成员数，不是组织细胞比例
levels = ['High confidence', 'Moderate confidence', 'Low confidence', 'Unresolved']
composition = composition.reindex(columns=levels, fill_value=0)
fig, ax = plt.subplots(figsize=(7.2, 3.5))
bottom = np.zeros(len(composition))
for level, color in zip(levels, ['#287D8E', '#81B6BA', '#D5B98D', '#B6B6B6']):
    counts = composition[level].to_numpy()
    ax.bar(composition.index, counts, bottom=bottom, color=color, label=level, width=.7)
    bottom += counts
ax.set_ylabel('Frozen member states')
ax.set_title('Candidate-program confidence within legacy E1-E7', fontsize=9)
ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1), fontsize=7, frameon=False)
ax.set_ylim(0, 13)
ax.set_yticks(range(0, 14, 2))
fig.subplots_adjust(left=.10, right=.68, bottom=.15, top=.86)
plots.append((fig, 'Ecotype_member_confidence'))

# Step4. 每父类型一个检出率热图，使用原作者审计供者等权值；不把NA补0
for cell_type in master.CellType.unique():
    data = markers[markers.FrozenParentCellType.eq(cell_type)].copy()
    data['RowLabel'] = data.Dataset.str.replace('Colon_GSE214695', 'Colon', regex=False).str.replace('Ileum_IBDverse', 'Ileum', regex=False) + ' / ' + data.CellState
    measured = data[data.SourceMeasured.eq(True)]
    matrix = measured.pivot(index='RowLabel', columns='Gene', values='DetectionFraction')
    all_rows = data.RowLabel.drop_duplicates()
    all_genes = data.Gene.drop_duplicates()
    matrix = matrix.reindex(index=all_rows, columns=all_genes)
    fig, ax = plt.subplots(figsize=(max(7.2, len(all_genes) * .14 + 2), max(3.0, len(all_rows) * .24 + 1.2)))
    cmap = mpl.colormaps['Blues'].copy()
    cmap.set_bad('#D9D9D9')
    mesh = ax.imshow(np.ma.masked_invalid(matrix.to_numpy()), cmap=cmap, aspect='auto', vmin=0, vmax=1)
    ax.set_xticks(range(len(all_genes)), all_genes, rotation=90, rotation_mode='anchor', ha='right', fontsize=7)
    ax.set_yticks(range(len(all_rows)), all_rows, fontsize=7)
    ax.set_title(f'{cell_type}: donor-balanced marker detection; grey = not measured', fontsize=9)
    fig.colorbar(mesh, ax=ax, fraction=.02, pad=.02, label='Mean donor detection fraction')
    fig.subplots_adjust(left=.24, bottom=.30, right=.94, top=.87)
    plots.append((fig, f'Marker_detection_{cell_type}'))

# Step5. 导出前做几何检查，每次PDF导出后做文字与碰撞QA
for fig, name in plots:
    target = figure_dir / name
    fig.canvas.draw()
    require_matplotlib_panel_alignment(fig, json_out=str(target) + '.alignment.json', tolerance_pt=1.5, gutter_tolerance_pt=1.5, strict=True)
    fig.savefig(str(target) + '.pdf')
    fig.savefig(str(target) + '.svg')
    fig.savefig(str(target) + '.png', dpi=600)
    plt.close(fig)
    subprocess.run([sys.executable, str(qa / 'audit_pdf_text.py'), str(target) + '.pdf'], check=True)
    subprocess.run([sys.executable, str(qa / 'audit_figure_collisions.py'), str(target) + '.pdf', '--json-out', str(target) + '.collision-audit.json', '--overlay-pdf', str(target) + '.collision-audit.pdf'], check=True)
    print('Rendered; final visual inspection still required:', name)
