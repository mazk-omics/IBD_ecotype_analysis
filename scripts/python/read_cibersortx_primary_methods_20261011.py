# 从 Europe PMC 的官方全文接口读取 CIBERSORTx 与 EcoTyper 原始论文，保存 XML 来源及方法段落。
# 本机 cibersortx_review 环境运行；只整理公开文献，不连接服务器、不修改算法或研究数据。
from pathlib import Path
import urllib.request
import urllib.error
import xml.etree.ElementTree as ET
import hashlib
import json
import csv

# Step1. 两篇原始论文分别记录，避免把 EcoTyper 的下游筛选归给 HiRes
project = Path(__file__).resolve().parents[2]
outdir = project / 'CIBERSORTx_server_plan_20261011/paper_methods'
outdir.mkdir(exist_ok=True)
papers = [('PMC6610714', 'CIBERSORTx_2019'), ('PMC8526411', 'EcoTyper_Cell_2021')]
records = []
sections = []
supplements = []
for accession, label in papers: # 两篇论文使用相同的来源与段落格式，分别保留身份
    url = 'https://www.ebi.ac.uk/europepmc/webservices/rest/' + accession + '/fullTextXML'
    destination = outdir / (accession + '.xml')
    if not destination.exists():
        request = urllib.request.Request(url, headers={'User-Agent': 'IBD-method-audit/1.0'})
        try:
            with urllib.request.urlopen(request, timeout=40) as response:
                payload = response.read()
        except urllib.error.HTTPError as error:
            records.append({'accession': accession, 'label': label, 'url': url, 'status': 'HTTP_' + str(error.code)})
            continue
        ET.fromstring(payload) # 先确认返回 XML，而非网页验证页面
        destination.write_bytes(payload)
    payload = destination.read_bytes()
    article = ET.fromstring(payload)
    title = ' '.join(article.find('.//article-title').itertext())
    doi = next((node.text for node in article.findall('.//article-id') if node.attrib.get('pub-id-type') == 'doi'), '')
    records.append({'accession': accession, 'label': label, 'url': url, 'status': 'XML_saved', 'title': title, 'doi': doi, 'sha256': hashlib.sha256(payload).hexdigest()})

    # Step2. 保留标题与段落次序；数学原式需回到 XML 的 MathML 核对，不能只按平文本猜公式
    section_number = 0
    for section in article.findall('.//sec'):
        heading = section.find('title')
        if heading is None:
            continue
        section_number += 1
        heading_text = ' '.join(heading.itertext())
        direct_paragraphs = section.findall('p')
        paragraph_text = '\n\n'.join(' '.join(node.itertext()) for node in direct_paragraphs)
        if not paragraph_text:
            continue
        sections.append({'paper': label, 'section_number': section_number, 'id': section.attrib.get('id', ''), 'heading': heading_text, 'paragraphs': paragraph_text})
    for supplement in article.findall('.//supplementary-material'):
        hrefs = [node.attrib.get('{http://www.w3.org/1999/xlink}href', '') for node in supplement.iter() if '{http://www.w3.org/1999/xlink}href' in node.attrib]
        supplements.append({'paper': label, 'id': supplement.attrib.get('id', ''), 'description': ' '.join(supplement.itertext()), 'hrefs': hrefs})

# Step3. 原文留作本地阅读材料，索引明确引用来源；不把原文重新发布为项目方法
(outdir / 'primary_paper_sources.json').write_text(json.dumps(records, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
(outdir / 'supplement_links.json').write_text(json.dumps(supplements, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
with (outdir / 'primary_paper_section_index.tsv').open('w', encoding='utf-8', newline='') as stream:
    writer = csv.DictWriter(stream, fieldnames=['paper', 'section_number', 'id', 'heading', 'paragraphs'], delimiter='\t')
    writer.writeheader()
    writer.writerows(sections)
for accession, label in papers:
    selected = [row for row in sections if row['paper'] == label]
    text = '\n\n'.join('## ' + row['heading'] + '\n\n' + row['paragraphs'] for row in selected)
    (outdir / (label + '_sections.txt')).write_text(text + '\n', encoding='utf-8')
print('Paper retrieval:', [(row['label'], row['status']) for row in records])
print('Sections with paragraphs:', len(sections))
print('Supplement records:', supplements)
