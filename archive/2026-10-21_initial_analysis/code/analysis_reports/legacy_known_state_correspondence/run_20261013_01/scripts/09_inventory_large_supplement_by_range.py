# 仅请求大Source Data ZIP的末尾目录，不下载整个大归档。
# 若服务器不支持HTTP206范围响应则停止，避免大文件传输。
from pathlib import Path
import urllib.request
import struct
import pandas as pd

out = Path(__file__).resolve().parents[1]
url = 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-023-40156-6/MediaObjects/41467_2023_40156_MOESM4_ESM.zip'
rows = []
status = 'NOT_RETRIEVED'
try:
    request = urllib.request.Request(url, headers={'Range': 'bytes=-131072'})
    with urllib.request.urlopen(request, timeout=20) as response:
        if response.status != 206:
            raise ValueError('Server did not return206; body not read')
        data = response.read(131072 + 1)
        content_range = response.headers.get('Content-Range', '')
    assert len(data) <= 131072
    position = data.rfind(b'PK\x05\x06')
    assert position >= 0, 'EOCD not present in permitted tail'
    eocd = struct.unpack_from('<4s4H2LH', data, position)
    total_entries, central_size, central_offset = eocd[4], eocd[5], eocd[6]
    start_offset = int(content_range.split(' ')[1].split('-')[0])
    cursor = central_offset - start_offset
    assert cursor >= 0 and central_size <= len(data), 'Directory outside small permitted tail'
    for i in range(total_entries):
        assert data[cursor:cursor + 4] == b'PK\x01\x02'
        entry = struct.unpack_from('<4s6H3L5H2L', data, cursor)
        name_n, extra_n, comment_n = entry[10], entry[11], entry[12]
        name = data[cursor + 46:cursor + 46 + name_n].decode('utf-8', errors='replace')
        rows.append({'Member': name, 'CompressedBytes': entry[8], 'UncompressedBytes': entry[9], 'LocalHeaderOffset': entry[16], 'CompressionMethod': entry[4], 'CRC32': entry[7], 'GeneralPurposeFlags': entry[3]})
        cursor += 46 + name_n + extra_n + comment_n
    status = 'TAIL_DIRECTORY_ONLY'
    pd.DataFrame(rows).to_csv(out / 'audit/GarridoTrigo2023_RANGE_ZIP_MEMBER_INVENTORY.tsv', sep='\t', index=False)
    print(pd.DataFrame(rows).to_string(index=False))
except Exception as error:
    print('Range inventory stopped:', error)
    status = 'NOT_RETRIEVED: ' + str(error)
(out / 'audit/GarridoTrigo2023_RANGE_STATUS.txt').write_text(status, encoding='utf-8')
