# 对最终导出的 recovery PDF 执行 nature-figure 字号与碰撞审计，并生成目视预览。
# 在本地 .tools/cibersortx_review 环境运行；读取 figures，审计 JSON、叠加 PDF、日志写入 figure_qa。
# 这些检查针对排版，不证明统计有效性；出现失败时必须修图后重新审计。
from pathlib import Path
import json
import subprocess
import sys
import pymupdf

# Step1. 检查绘图脚本配置
project = Path(__file__).resolve().parents[3]
root = project / 'IBD_ecotype_analysis/06_validation/GSE214695_scRNA_recovery_20261009'
skill = Path.home() / '.codex/skills/nature-figure/scripts'
qa = root / 'figure_qa'
qa.mkdir(exist_ok=True)
preflight = subprocess.run([sys.executable, '-X', 'utf8', str(skill / 'validate_figure.py'), str(Path(__file__).parent / 'plot_gse214695_recovery_nature.py'), '--backend', 'python', '--json'], capture_output=True, text=True, encoding='utf-8')
(qa / 'source_preflight.json').write_text(preflight.stdout, encoding='utf-8')
assert preflight.returncode == 0, preflight.stdout
comparison_preflight = subprocess.run([sys.executable, '-X', 'utf8', str(skill / 'validate_figure.py'), str(Path(__file__).parent / 'plot_gse214695_ileum_comparison.py'), '--backend', 'python', '--json'], capture_output=True, text=True, encoding='utf-8')
(qa / 'comparison_source_preflight.json').write_text(comparison_preflight.stdout, encoding='utf-8')
assert comparison_preflight.returncode == 0, comparison_preflight.stdout

# Step2. 每张最终 PDF 独立检查；叠加文件只作 QA，不混入汇报图件
rows = []
for pdf in sorted((root / 'figures').glob('*.pdf')):
    text_audit = subprocess.run([sys.executable, '-X', 'utf8', str(skill / 'audit_pdf_text.py'), str(pdf), '--min-pt', '5', '--json'], capture_output=True, text=True, encoding='utf-8')
    (qa / f'{pdf.stem}_text.json').write_text(text_audit.stdout, encoding='utf-8')
    collision = subprocess.run([sys.executable, '-X', 'utf8', str(skill / 'audit_figure_collisions.py'), str(pdf), '--json-out', str(qa / f'{pdf.stem}_collision.json'), '--overlay-pdf', str(qa / f'{pdf.stem}_collision_overlay.pdf')], capture_output=True, text=True, encoding='utf-8')
    (qa / f'{pdf.stem}_collision.log').write_text(collision.stdout + collision.stderr, encoding='utf-8')
    rows.append({'figure': pdf.name, 'text_exit_code': text_audit.returncode, 'collision_exit_code': collision.returncode})
    print(rows[-1], flush=True)
    doc = pymupdf.open(pdf)
    page = doc[0]
    page.get_pixmap(matrix=pymupdf.Matrix(2, 2), alpha=False).save(qa / f'{pdf.stem}_preview.png')
    doc.close()
(qa / 'audit_exit_codes.json').write_text(json.dumps(rows, indent=2), encoding='utf-8')
assert all(x['text_exit_code'] == 0 and x['collision_exit_code'] == 0 for x in rows), '存在排版审计失败，须查看 JSON/叠加图后修正'
