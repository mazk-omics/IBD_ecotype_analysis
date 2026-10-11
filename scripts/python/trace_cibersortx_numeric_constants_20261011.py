# 核对 HiRes 数值调度、抽样、排序和回归函数引用的常量，只反汇编已定位的算法函数。
# 本机 cibersortx_review 环境运行；读已固定摘要的 ELF，写常量与指令证据，不执行供应商程序。
from pathlib import Path
import hashlib
import json
import struct
import csv
import math
from elftools.elf.elffile import ELFFile
from capstone import Cs, CS_ARCH_X86, CS_MODE_64, x86_const

# Step1. 二进制摘要与 Ghidra 映射地址固定，防止把不同镜像的地址混用
project = Path(__file__).resolve().parents[2]
outdir = project / 'CIBERSORTx_server_plan_20261011'
source = outdir / 'vendor_reference/hires/src/CIBERSORTxHiRes'
source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
assert source_hash == '87511e7db27a973757e12d401b5184d306a9e85d24ae0748ab7199a640ae5f89'
image_base = 0x100000
with (outdir / 'decompiled_numeric/CIBERSORTxHiRes/numerical_functions.tsv').open(encoding='utf-8', newline='') as stream:
    function_rows = list(csv.DictReader(stream, delimiter='\t'))
selected_names = ['CIBERSORTxHiRes', 'doBootstrap', 'runHiRes', 'findBestRSSThres', 'utility::sequence', 'utility::order<double>', 'LinearRegression::calcSlope', 'LinearRegression::calcIntercept', 'LinearRegression::predict', 'stats::adjustPVals', 'stats::fisherComb', 'stats::rowpvals']
functions = {row['Name']: (int(row['Address'], 16), int(row['BodyBytes'])) for row in function_rows if row['Name'] in selected_names}
observed = []

# Step2. 浮点读指令给出真实地址；不能按 Ghidra 的重叠全局符号名称猜常量偏移
with source.open('rb') as stream:
    elf = ELFFile(stream)
    sections = list(elf.iter_sections())
    disassembler = Cs(CS_ARCH_X86, CS_MODE_64)
    disassembler.detail = True
    for name, (address, size) in functions.items():
        section = next(item for item in sections if item['sh_addr'] <= address - image_base < item['sh_addr'] + item['sh_size'])
        offset = address - image_base - section['sh_addr']
        code = section.data()[offset:offset + size]
        for instruction in disassembler.disasm(code, address):
            for operand in instruction.operands:
                if operand.type != x86_const.X86_OP_MEM or operand.mem.base != x86_const.X86_REG_RIP or operand.size != 8 or instruction.mnemonic not in ['movsd', 'mulsd', 'divsd', 'addsd', 'subsd', 'comisd', 'ucomisd']:
                    continue
                target = instruction.address + instruction.size + operand.mem.disp
                data_section = next((item for item in sections if item['sh_addr'] <= target - image_base < item['sh_addr'] + item['sh_size']), None)
                if data_section is None or data_section.name != '.rodata':
                    continue # 全局可变对象的初始值不当作算法常量
                data_offset = target - image_base - data_section['sh_addr']
                raw = data_section.data()[data_offset:data_offset + 8]
                value = struct.unpack('<d', raw)[0]
                observed.append({'function': name, 'instruction_address': hex(instruction.address), 'instruction': instruction.mnemonic + ' ' + instruction.op_str, 'constant_address': hex(target), 'section': data_section.name, 'raw_little_endian_hex': raw.hex(), 'float64_value': value if math.isfinite(value) else str(value)})

# Step3. 保留地址、原始字节及解释值，后续实现可逐项核对
result = {'binary_sha256': source_hash, 'ghidra_image_base': hex(image_base), 'program_executed': False, 'constants': observed}
(outdir / 'hires_numeric_constant_evidence.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
print('Float64 constant references:', len(observed))
for row in observed:
    if row['function'] != 'doBootstrap':
        continue
    print(row['function'], row['constant_address'], row['float64_value'], row['instruction_address'])
