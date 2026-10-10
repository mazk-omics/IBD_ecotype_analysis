# IBD EcoTyper project instructions

# Required research code style

Follow `Profile/CODE_STYLE.md` for every new or modified analysis script. The user
requires linear imperative research scripting: readable top-to-bottom execution,
short Chinese Step comments, compact direct commands, visible intermediate
objects and outputs, and parameters near their actual use. Avoid unnecessary
wrappers, classes, pipeline/configuration frameworks, mechanical line wrapping,
and semicolon-chained independent operations. Keep only algorithmically necessary
loops and checks preventing leakage, mismatches, or data loss. Style-only edits
must preserve methods, parameters, filtering, and results. Do not automatically
rewrite existing scripts outside the authorized task. This rule applies throughout
the EcoTyper project and its eventual publication-oriented code review.

# Server access rules

The user is concerned that repeated short SSH connections to Xiyou can trigger
the provider's IP blocking. The user initially withdrew permission for autonomous server login, then on
2026-10-06 explicitly authorized completing the corrected IBDverse pilot and
colon candidate, then separately authorized GSE282122 joint and ileum pilots and
a SCP259 pilot, with local result retrieval. Each scoped run must use one
persistent SSH session. Outside explicitly authorized server tasks, work locally and provide
artifacts for the user to upload through their existing connection.

If the user explicitly reauthorizes server access for a task:

- Establish at most one persistent SSH shell using
  `IBD_ecotype_analysis/scripts/powershell/Open-XiyouPersistentSession.ps1`.
  Keep the returned exec session ID and use `write_stdin` for subsequent remote
  commands. Do not launch a new ssh process for every tool call.
- Upload small scripts through that same session rather than opening separate
  SCP/SFTP connections. For large transfers, agree on an appropriate persistent
  transfer method first.
- There is no automatic reconnection loop. If authentication, host verification,
  connection setup, or the established session fails, stop remote work and report
  the failure. Do not try more keys, passwords, another port, or repeated probes.
- Use the configured Xiyou key, strict host verification, and a single connection
  attempt. Never print private keys or authentication secrets.
- A Python/R/script error after successful login is not a reason to reconnect.
  Diagnose it locally where possible or through the same live session.

The Windows native OpenSSH client does not supply the Unix ControlMaster sharing
used by many Linux SSH workflows. Keepalive options do not themselves share a
connection between multiple ssh processes.
