param(
    [int]$LocalTransferPort = 0,
    [int]$RemoteTransferPort = 0
)

$ErrorActionPreference = 'Stop'
$xiyouSshPath = Join-Path $env:WINDIR 'System32\OpenSSH\ssh.exe'
if (-not (Test-Path -LiteralPath $xiyouSshPath)) {
    throw "Windows OpenSSH client is unavailable: $xiyouSshPath"
}

# Prevent simultaneous sessions launched through this helper. There is no retry
# loop: after a lost session, the agent must obtain explicit reauthorization.
$xiyouSessionMutex = [System.Threading.Mutex]::new($false, 'Local\Codex_IBD_Xiyou_PersistentSession')
$xiyouMutexOwned = $false
try {
    try {
        $xiyouMutexOwned = $xiyouSessionMutex.WaitOne(0)
    } catch [System.Threading.AbandonedMutexException] {
        $xiyouMutexOwned = $true
    }
    if (-not $xiyouMutexOwned) {
        throw 'A Xiyou session launched by this helper is already active. Reuse its exec session ID.'
    }
    Write-Output 'Opening ONE persistent SSH shell; no automatic reconnect or separate SCP calls.'
    $xiyouForwardArgs = @()
    if ($LocalTransferPort -gt 0 -and $RemoteTransferPort -gt 0) {
        if ($LocalTransferPort -gt 65535 -or $RemoteTransferPort -gt 65535) { throw 'Invalid transfer port.' }
        $xiyouForwardArgs = @('-L', "127.0.0.1:${LocalTransferPort}:127.0.0.1:${RemoteTransferPort}", '-o', 'ExitOnForwardFailure=yes')
    }
    & $xiyouSshPath @xiyouForwardArgs -tt `
        -o BatchMode=yes `
        -o IdentitiesOnly=yes `
        -o PreferredAuthentications=publickey `
        -o ConnectionAttempts=1 `
        -o ConnectTimeout=15 `
        -o StrictHostKeyChecking=yes `
        -o ServerAliveInterval=60 `
        -o ServerAliveCountMax=3 `
        Xiyou
    $xiyouExitCode = $LASTEXITCODE
    if ($xiyouExitCode -ne 0) {
        Write-Error "SSH ended with code $xiyouExitCode. Stop; do not reconnect automatically."
    }
} finally {
    if ($xiyouMutexOwned) { $xiyouSessionMutex.ReleaseMutex() }
    $xiyouSessionMutex.Dispose()
}
