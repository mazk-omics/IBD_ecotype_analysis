# 从已核对ZIP目录按范围读取两个小型补充图Source Data，不下载完整归档。
# 仅用于已知巨噬状态替代解释；每成员<10 MB且核对ZIP长度/CRC。
from pathlib import Path
import urllib.request
import struct
import zlib
import hashlib
import pandas as pd

# Step1. 限定有明确用途的两个小表
out = Path(__file__).resolve().parents[1]
directory = pd.read_csv(out / 'audit/GarridoTrigo2023_RANGE_ZIP_MEMBER_INVENTORY.tsv', sep='\t')
selected = directory[directory.Member.str.strip().isin(['Source_data_file_supplementary_figure_4.xlsx', 'Source_data_file_supplementary_figure_5.xlsx'])]
url = 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-023-40156-6/MediaObjects/41467_2023_40156_MOESM4_ESM.zip'
logs = []
schemas = []
for row in selected.itertuples():
    assert row.CompressedBytes <= 10 * 1024 * 1024 and row.UncompressedBytes <= 15 * 1024 * 1024
    offset = row.LocalHeaderOffset
    header_request = urllib.request.Request(url, headers={'Range': f'bytes={offset}-{offset + 4095}'})
    with urllib.request.urlopen(header_request, timeout=20) as response:
        assert response.status == 206
        header = response.read(4096)
    info = struct.unpack_from('<4s5H3L2H', header)
    assert info[0] == b'PK\x03\x04'
    start = offset + 30 + info[9] + info[10]
    request = urllib.request.Request(url, headers={'Range': f'bytes={start}-{start + row.CompressedBytes - 1}'})
    with urllib.request.urlopen(request, timeout=30) as response:
        assert response.status == 206
        content = response.read(row.CompressedBytes + 1)
    assert len(content) == row.CompressedBytes
    if row.CompressionMethod == 8:
        content = zlib.decompress(content, -15)
    # 使用central directory的CRC；data descriptor模式的本地header CRC可为0。
    assert len(content) == row.UncompressedBytes and zlib.crc32(content) == row.CRC32
    path = out / 'literature_resources' / ('GarridoTrigo2023__' + row.Member.strip())
    path.write_bytes(content)
    workbook = pd.ExcelFile(path, engine='openpyxl')
    for sheet in workbook.sheet_names:
        preview = pd.read_excel(workbook, sheet_name=sheet, header=None, nrows=5)
        schemas.append({'File': path.name, 'Sheet': sheet, 'PreviewJSON': preview.fillna('').to_json(orient='values', force_ascii=False)})
    logs.append({'URL': url, 'Member': row.Member, 'Bytes': len(content), 'SHA256': hashlib.sha256(content).hexdigest(), 'Transfer': 'HTTP206 exact member range; no whole ZIP', 'ZIP_CRC_OK': True})
    pd.DataFrame(logs).to_csv(out / 'audit/TARGETED_SOURCE_TABLE_TRANSFER_LOG.tsv', sep='\t', index=False)
    pd.DataFrame(schemas).to_csv(out / 'audit/GARRIDO_SOURCE_TABLE_SCHEMA.tsv', sep='\t', index=False)
    print(path.name, 'sheets:', workbook.sheet_names, flush=True)
