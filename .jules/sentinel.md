## 2025-05-14 - Command Injection in Perl scripts via samtools and Rscript
**Vulnerability:** Command injection was possible through file paths and command-line arguments passed to external tools (`samtools`, `Rscript`) using backticks, piped `open`, and single-argument `system` calls.
**Learning:** Legacy Perl code often uses 2-argument `open` and backticks for convenience, which implicitly invoke the shell. This is particularly dangerous in bioinformatics tools that process user-provided filenames.
**Prevention:**
1. Always use the 3-argument form of `open` for files: `open $fh, "<", $filename`.
2. Use the multi-argument form of `open` for pipes: `open $fh, "-|", "command", @args`.
3. Use the list form of `system`: `system($cmd, @args)`.
4. Avoid backticks when the command contains variables; use multi-argument `open` and read from the filehandle instead.
