## 2025-05-14 - Command Injection in Perl scripts via samtools and Rscript
**Vulnerability:** Command injection was possible through file paths and command-line arguments passed to external tools (`samtools`, `Rscript`) using backticks, piped `open`, and single-argument `system` calls.
**Learning:** Legacy Perl code often uses 2-argument `open` and backticks for convenience, which implicitly invoke the shell. This is particularly dangerous in bioinformatics tools that process user-provided filenames.
**Prevention:**
1. Always use the 3-argument form of `open` for files: `open $fh, "<", $filename`.
2. Use the multi-argument form of `open` for pipes: `open $fh, "-|", "command", @args`.
3. Use the list form of `system`: `system($cmd, @args)`.
4. Avoid backticks when the command contains variables; use multi-argument `open` and read from the filehandle instead.

## 2025-01-24 - Robust Pipe Error Handling and Binary Validation
**Vulnerability:** Silent failures in external command pipes and missing validation of external binaries could lead to undefined behavior or processing of incomplete data.
**Learning:** In Perl, `open` for a pipe might succeed even if the command eventually fails. The true exit status of the piped process is only reliably caught when the filehandle is closed.
**Prevention:**
1. Always check the return value of `close` on piped filehandles: `close($fh) or die "Pipe failed: $!"`.
2. Explicitly validate the existence and executability of external binaries in the `PATH` at startup, and use their absolute paths to ensure the intended tools are being invoked.
