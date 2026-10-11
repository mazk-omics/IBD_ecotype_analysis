# 本机核对 runHiRes 参数存储与机器指令序言，读固定 ELF 和补参输出，写独立证据。
# 使用 cibersortx_review 环境；只读取该数值函数的开头，不执行程序或扫描账户字符串。
from pathlib import Path
import hashlib
import json
import csv
from elftools.elf.elffile import ELFFile
from capstone import Cs, CS_ARCH_X86, CS_MODE_64

# Step1. 固定文件摘要，保证寄存器证据属于同一个编译版本
project = Path(__file__).resolve().parents[2]
outdir = project / 'CIBERSORTx_server_plan_20261011'
source = outdir / 'vendor_reference/hires/src/CIBERSORTxHiRes'
source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
assert source_hash == '87511e7db27a973757e12d401b5184d306a9e85d24ae0748ab7199a640ae5f89'
typed = outdir / 'decompiled_numeric/CIBERSORTxHiRes_typed'
with (typed / 'runHiRes_parameter_storage.tsv').open(encoding='utf-8', newline='') as stream:
    parameters = list(csv.DictReader(stream, delimiter='\t'))
assert len(parameters) == 14
assert parameters[6]['Storage'] == 'XMM0_Qa:8' and parameters[7]['Storage'] == 'XMM1_Qa:8'

# Step2. 只保存数值函数开头 40 条指令，独立检查整数和浮点寄存器的入栈位置
address = 0x10cc9b
with source.open('rb') as stream:
    elf = ELFFile(stream)
    section = next(item for item in elf.iter_sections() if item['sh_addr'] <= address - 0x100000 < item['sh_addr'] + item['sh_size'])
    offset = address - 0x100000 - section['sh_addr']
    code = section.data()[offset:offset + 400]
    disassembler = Cs(CS_ARCH_X86, CS_MODE_64)
    instructions = list(disassembler.disasm(code, address))[:40]
prologue = [{'address': hex(item.address), 'instruction': item.mnemonic + ' ' + item.op_str} for item in instructions]

# Step3. 区分符号类型、分析器存储假设和实际指令，不把补参命名当作原作者命名
text = (typed / 'runHiRes_typed.c').read_text(encoding='utf-8')
usage = {row['Name']: text.count(row['Name']) - 1 for row in parameters}
evidence = {'binary_sha256': source_hash, 'function_address': hex(address), 'program_executed': False, 'analyst_types_not_original_names': True, 'parameters': parameters, 'prologue': prologue, 'decompiled_body_identifier_occurrences': usage, 'typed_pseudocode_sha256': hashlib.sha256(text.encode('utf-8')).hexdigest()}
(outdir / 'hires_parameter_abi_evidence.json').write_text(json.dumps(evidence, indent=2) + '\n', encoding='utf-8')
print('Parameters:', len(parameters), '; floating-point slots:', parameters[6]['Storage'], parameters[7]['Storage'])
print('Only declaration, no recovered body use:', [name for name, count in usage.items() if count == 0])
print('Entry register stores:', [row for row in prologue if row['instruction'].startswith('mov') and any(reg in row['instruction'] for reg in ['edi', 'rsi', 'rdx', 'rcx', 'r8d', 'r9d', 'xmm0', 'xmm1'])])
