# 本机使用 Ghidra 对 HiRes 和两个 Fractions 编译版本进行数值函数审计，不运行供应商程序。
# 工具位于 .tools/cibersortx_static；输出静态项目、函数目录及伪代码，认证相关内容不导出。
$taskProject = Split-Path (Split-Path $PSScriptRoot -Parent) -Parent
$taskRoot = Split-Path $taskProject -Parent
$taskTools = Join-Path $taskRoot '.tools\cibersortx_static'
$taskAudit = Join-Path $taskProject 'CIBERSORTx_server_plan_20261011'
$taskWorkspace = Join-Path $taskAudit 'decompiler_workspace'
$taskHeadless = Join-Path $taskTools 'ghidra_12.1.4_PUBLIC\support\analyzeHeadless.bat'

# Step1. Java 与缓存仅针对当前进程设置，避免改动系统环境
$env:JAVA_HOME = Join-Path $taskTools 'jdk-25.0.4.1+1'
$env:GHIDRA_HEADLESS_MAXMEM = '2G'
$env:GHIDRA_HEADLESS_JAVA_OPTIONS = "-Dapplication.settingsdir=$taskTools\settings -Dapplication.cachedir=$taskTools\cache -Dapplication.tempdir=$taskTools\temp"
$taskScriptPath = Join-Path $taskProject 'scripts\java'
$taskOutput = Join-Path $taskAudit 'decompiled_numeric'
New-Item -ItemType Directory -Force -Path $taskWorkspace | Out-Null

# Step2. 已建立的分析项目直接复用；首次运行则导入，三个不同编译文件分别保存
$taskHiResInput = Join-Path $taskAudit 'vendor_reference\hires\src\CIBERSORTxHiRes'
$taskHiResMode = @('-import', $taskHiResInput)
if (Test-Path -LiteralPath (Join-Path $taskWorkspace 'CibersortxNumericAudit.gpr')) { $taskHiResMode = @('-process', 'CIBERSORTxHiRes', '-noanalysis') }
& $taskHeadless $taskWorkspace CibersortxNumericAudit @taskHiResMode -scriptPath $taskScriptPath -postScript AuditCibersortxNumerics.java $taskOutput -max-cpu 2 -analysisTimeoutPerFile 300 *> (Join-Path $taskAudit 'ghidra_hires_numerical_pass.log')
if ($LASTEXITCODE -ne 0) { throw 'HiRes static inspection failed; read local log.' }
Write-Output 'HiRes numerical pass completed.'

$taskFractions = @('hires', 'fractions')
foreach ($taskImage in $taskFractions) { # 镜像内的两个不同 Fractions 哈希需分别分析
    $taskInput = Join-Path $taskAudit "vendor_reference\$taskImage\src\CIBERSORTxFractions"
    $taskName = "CibersortxFractions_$taskImage"
    $taskMode = @('-import', $taskInput)
    if (Test-Path -LiteralPath (Join-Path $taskWorkspace "$taskName.gpr")) { $taskMode = @('-process', 'CIBERSORTxFractions', '-noanalysis') }
    & $taskHeadless $taskWorkspace $taskName @taskMode -scriptPath $taskScriptPath -postScript AuditCibersortxNumerics.java $taskOutput "CIBERSORTxFractions_$taskImage" -max-cpu 2 -analysisTimeoutPerFile 300 -log (Join-Path $taskAudit "ghidra_fractions_$taskImage.log") *> (Join-Path $taskAudit "ghidra_fractions_$taskImage.console.log")
    if ($LASTEXITCODE -ne 0) { throw "Fractions static inspection failed: $taskImage" }
    Write-Output "Fractions numerical pass completed: $taskImage"
}
