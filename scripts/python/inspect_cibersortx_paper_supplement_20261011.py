# 读取 CIBERSORTx 原论文的出版社补充 PDF，保留逐页文本、来源摘要与公式核对入口。
# 本机 bundled Python（pdfplumber）运行；论文原文仅作本地阅读，索引写入静态审计目录。
from pathlib import Path
import urllib.request
import hashlib
import json
import pdfplumber

# Step1. 固定官方公开补充材料来源，已有 PDF 不重复下载
project = Path(__file__).resolve().parents[2]
outdir = project / 'CIBERSORTx_server_plan_20261011'
papers = outdir / 'paper_methods'
papers.mkdir(exist_ok=True)
url = 'https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41587-019-0114-2/MediaObjects/41587_2019_114_MOESM1_ESM.pdf'
source = papers / 'CIBERSORTx_2019_Supplement.pdf'
if not source.exists():
    with urllib.request.urlopen(url, timeout=45) as response:
        payload = response.read()
    assert payload.startswith(b'%PDF'), '来源没有返回 PDF，不能将 HTML 错误页当作补充材料'
    source.write_bytes(payload)
payload = source.read_bytes()
assert payload.startswith(b'%PDF')

# Step2. 保存真实 PDF 页码；正文页码比 PDF 页码少一页，公式需看渲染图
with pdfplumber.open(source) as document:
    texts = [page.extract_text() or '' for page in document.pages]
    metadata = document.metadata
page_text = '\n\n'.join('## PDF page ' + str(i + 1) + '\n' + text for i, text in enumerate(texts))
(papers / 'CIBERSORTx_2019_Supplement_pages.txt').write_text(page_text + '\n', encoding='utf-8')
assert len(texts) == 43
assert '12.3.6.' in texts[35] and '12.3.7.' in texts[36]
print('Supplement pages:', len(texts), '; refinement: PDF 36 / printed 35; smoothing: PDF 37 / printed 36')

# Step3. 只把来源与章节定位纳入项目审计，原文不随代码重新发布
sections = [{'section': '12.2', 'pdf_pages': [34], 'topic': 'Group-mode NNLS, BH and adaptive geometric CV'}, {'section': '12.3.1–12.3.4', 'pdf_pages': [34, 35, 36], 'topic': 'Global detection, sorting, split fit and differential tests'}, {'section': '12.3.5–12.3.6', 'pdf_pages': [36, 37], 'topic': 'Coefficient refinement and weak fold-change replacement'}, {'section': '12.3.7–12.3.8', 'pdf_pages': [37], 'topic': 'Smoothing, OLS anchors and sample order'}, {'section': '12.4.3', 'pdf_pages': [38, 39], 'topic': 'S-mode synthetic mixtures and adjusted signature'}]
evidence = {'doi': '10.1038/s41587-019-0114-2', 'publisher_url': url, 'local_pdf': str(source), 'sha256': hashlib.sha256(payload).hexdigest(), 'bytes': len(payload), 'pages': len(texts), 'pdf_metadata': metadata, 'sections': sections, 'visually_checked_pdf_pages': [36, 37], 'visual_check_note': 'Rendered with Poppler and inspected: squared confidence term and absolute half-difference; edge padding notation contains row c although Q has N-W+1 rows. Static code uses the last valid window row.', 'ecotyper_primary_url': 'https://pmc.ncbi.nlm.nih.gov/articles/PMC8526411/', 'ecotyper_xml_status': 'Europe PMC fullTextXML returned HTTP 500; webpage indexed primary text used for downstream context only.'}
(outdir / 'paper_supplement_provenance.json').write_text(json.dumps(evidence, ensure_ascii=False, indent=2, default=str) + '\n', encoding='utf-8')
print('Source SHA256:', evidence['sha256'])
print('Evidence:', outdir / 'paper_supplement_provenance.json')
