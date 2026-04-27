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

## 2025-05-14 - Binary Planting in PATH Discovery
**Vulnerability:** Binary discovery logic that iterates through the `PATH` was susceptible to binary planting if the current directory (`.`) or an empty string was present in the `PATH`.
**Learning:** Even when manually searching the `PATH` instead of relying on the shell, one must explicitly ignore untrusted directories like the current working directory to prevent execution of malicious binaries placed there.
**Prevention:** Explicitly skip empty entries and `.` when iterating through `File::Spec->path()`.

## 2025-05-24 - Denial of Service via Division-by-Zero in Binning Logic
**Vulnerability:** The application was susceptible to crashes (DoS) when processing input that resulted in zero-valued divisors (e.g., `bootstrap`, `ywin`, `genome_size`) during data binning or progress reporting.
**Learning:** In bioinformatics tools, user-provided parameters or genomic dimensions often dictate binning logic. If these values are not strictly validated as non-zero, they can lead to runtime exceptions.
**Prevention:** Explicitly check if global counters or calculated divisors are non-zero before performing division or modulo operations. Use `List::Util::max(1, ...)` for safe defaults in progress indicators.
## 2026-04-24 - Division-by-Zero in Numerical Calculations and Progress Reporting
**Vulnerability:** The application was susceptible to script termination (DoS) due to division-by-zero when processing small or invalid genomic data, or when using low bootstrap values.
**Learning:** In Perl, division by zero is a fatal error. Genomic tools often perform divisions based on derived values (like genome size or median distance) which can be zero if input data is malformed or filter-heavy.
**Prevention:**
1. Always validate numeric command-line arguments and derived values (e.g., `$ywin`, `$genome_size`) before using them as divisors.
2. In loop-based progress reporting, ensure the modulo divisor is at least 1 using `List::Util::max(1, ...)` to handle small iteration counts safely.

## 2026-04-26 - Infinite Loops and Division-by-Zero in SNR Calculation
**Vulnerability:** The `sv::snr` function was susceptible to infinite loops and division-by-zero when processing uniform or insufficient relative entropy data.
**Learning:** Signal-to-Noise Ratio (SNR) calculations often involve partitioning data based on range. If the range is zero (all values identical) or data points are insufficient, calculations for resolution or weights can lead to zero steps in loops or division-by-zero errors.
**Prevention:** Explicitly validate that the data set has at least two elements and a non-zero range before proceeding with partitioning-based calculations.

## 2026-04-27 - Algorithmic Complexity (DoS) in SNR Calculations
**Vulnerability:** The `sv::snr` function exhibited $O(N^2)$ memory growth and unbounded CPU usage. It stored array slices in hashes for every iteration of a loop whose step size was fixed, leading to memory exhaustion and CPU hangs on large genomic datasets.
**Learning:** Partitioning algorithms that iterate through a range must have a bounded number of steps. Storing data subsets (like array slices) in long-lived structures within loops can lead to rapid memory depletion.
**Prevention:**
1. Always bound the number of iterations in search or optimization loops by making the step size proportional to the total range.
2. Use transient data references or direct indexing instead of storing redundant copies of data subsets in hashes or arrays within loops.
