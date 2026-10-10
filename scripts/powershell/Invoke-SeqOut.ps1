# 从本机独立 Python 3.13 环境调用 SeqOut，后面的参数原样传给官方 CLI。
# 例如：& .\IBD_ecotype_analysis\scripts\powershell\Invoke-SeqOut.ps1 search 'Crohn spatial transcriptomics' --organism 'Homo sapiens'
$seqoutProject = Split-Path (Split-Path (Split-Path $PSScriptRoot -Parent) -Parent) -Parent
$seqoutExecutable = Join-Path $seqoutProject '.tools/seqout/environment/Scripts/seqout.exe'
$env:PYTHONUTF8 = '1'
& $seqoutExecutable @args
