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

## 2026-04-28 - Infinite Loop DoS in Binary Search via NaN values
**Vulnerability:** The `sv::binary_search` function was susceptible to infinite loops when encountering `NaN` values in the input array.
**Learning:** In Perl, all numeric comparisons (`==`, `<`, `>`) with `NaN` return false. If a binary search loop depends solely on these comparisons to update its bounds, it will never terminate if it lands on a `NaN`.
**Prevention:** Ensure search loops have a terminal `else` or fallback condition that breaks the loop if no comparison is met, and sanitize numeric input data to exclude `NaN` or non-numeric values before processing.

## 2026-04-29 - Denial of Service in Genomic Range Processing
**Vulnerability:** The `sv::range` function used a point-expansion algorithm that converted genomic intervals into individual points. For large ranges, this caused massive memory allocation and CPU spikes ($O(\text{TotalRangeSize})$).
**Learning:** In Perl, lexical variables `$a` and `$b` can shadow the global variables used by the `sort` built-in, leading to broken comparisons in algorithms that rely on sorting. Genomic data often contains extremely large intervals that must be processed as intervals rather than collections of points.
**Prevention:** Use interval-merging algorithms with $O(N \log N)$ complexity for range operations. Avoid using `$a` and `$b` as lexical argument names in subroutines to prevent interference with the `sort` built-in.

## 2026-05-15 - Algorithmic Complexity DoS in Chromosome Ordering
**Vulnerability:** The chromosome ordering logic in `svre.pl` used nested loops ($O(N^2)$) to group chromosomes by size. When processing reference genomes with thousands of chromosomes of identical size (e.g., fragmented assemblies), this led to extreme CPU usage and memory exhaustion due to the list size exploding.
**Learning:** Naive grouping and sorting algorithms using nested `grep` calls can easily exhibit quadratic or worse complexity. Additionally, incorrect logic for filtering chromosomes (comparing a hash reference to a string) caused valid data to be dropped or misaligned.
**Prevention:** Always use efficient $O(N \log N)$ sorting functions (like Perl's `sort`) with appropriate tie-breakers for stability. Ensure filtering logic correctly uses `exists` or value checks instead of comparing references to strings.

## 2026-05-18 - Algorithmic Complexity DoS in Relative Entropy Calculation (sv::ric)
**Vulnerability:** The `sv::ric` function in `sv.pm` performed $O(N)$ iteration over the entire genome length (e.g., 250M iterations for human chr1) to bin genomic data, even when only a few data points were present. This led to extreme CPU usage and execution times (e.g., ~27 seconds for 100MB genome).
**Learning:** Naive data binning that iterates through every possible coordinate in a reference is extremely inefficient for sparse genomic data. High-performance genomic tools must iterate over observed data points and calculate spans/distances between them.
**Prevention:** Replace genome-wide loops with iterations over sorted data-containing keys from hashes. Use coordinate deltas between consecutive data points to maintain identical binning and bin-size tracking logic.

## 2026-05-20 - Algorithmic Complexity DoS in Coordinate Normalization
**Vulnerability:** The `sv::refmod` and `sv::refadd` functions used `while` loops to normalize genomic coordinates against a reference length. For extremely large coordinates or small reference lengths, this resulted in $O(N)$ complexity, leading to CPU exhaustion and script hangs (DoS).
**Learning:** Naive normalization using iterative subtraction/addition is dangerous when the input range is untrusted or unbounded.
**Prevention:** Always use $O(1)$ arithmetic (like the modulo operator or direct division) for coordinate normalization instead of loops.

## 2026-05-22 - Algorithmic Complexity DoS in sv::unwind_distance
**Vulnerability:** The `sv::unwind_distance` function performed a linear scan ($O(N)$) over an entire genomic bin range (which could be extremely large depending on user-provided `binsize`). This led to CPU exhaustion (DoS) when processing bins with large coordinate spans.
**Learning:** Iterating through every integer in a genomic range is a common DoS vector in bioinformatics tools. Replacing linear scans with iterations over observed data points (hash keys) protects against large ranges, but can introduce performance regressions if the entire dataset is scanned repeatedly.
**Prevention:** Use a combination of coordinate-based iteration and caching. Store sorted genomic coordinates for each chromosome and use binary search to locate starting positions within bins. This achieves $O(\log N + \text{DataPointsInBin})$ complexity, providing both security and high performance.

## 2026-05-24 - DoS via log(0) and Division-by-Zero in Entropy Calculations
**Vulnerability:** The application was susceptible to runtime crashes (DoS) when calculating Relative Information Content (RIC) or bootstrapping if input data resulted in zero-valued probabilities or missing bins in the global distribution.
**Learning:** In Perl, `log(0)` is a fatal error. Genomic data with extremely low coverage or skewed distributions can lead to zero-valued bins. Skipping these bins in summation is mathematically correct for entropy ($0 \log 0 = 0$) and prevents script termination.
**Prevention:** Always guard `log()` and division operations with checks for positive arguments and non-zero denominators, especially when derived from data-dependent distributions.

## 2026-05-30 - Denial of Service via Infinite Loops in Random Event Simulation
**Vulnerability:** The structural variation simulation functions (`sv::randomsize`, `sv::deletion`, `sv::inversion`, `sv::tandem`, `sv::duplication`) used `while` loops to find random sizes and positions that satisfied specific constraints. If the input parameters (e.g., sequence length vs. required event size) made these constraints impossible to meet, the script would enter an infinite loop, causing a DoS.
**Learning:** "Trial and error" random generation logic is dangerous when the search space can become empty based on user-provided or data-derived parameters.
**Prevention:** Refactor random generation to be deterministic in its termination. Calculate the valid range of values first, and then pick from that range in a single step (e.g., `pos = rand(max_pos)`). Ensure that minimum requirements are clamped to the maximum possible value to prevent invalid ranges.

## 2026-05-26 - Denial of Service via Infinite Loops and Division-by-Zero in Utility Functions
**Vulnerability:** The `sv::null_distribution`, `sv::null_distribution_from_file`, `sv::rms_variance`, and `sv::refmod` functions were susceptible to infinite loops or division-by-zero crashes when provided with non-positive window sizes, reference lengths, or zero-valued bin counts.
**Learning:** Utility functions in libraries often lack the strict validation present in the main application entry points. This can lead to vulnerabilities if the library is used by other tools or if validation is bypassed.
**Prevention:** Implement strict guard clauses in all library functions that perform iterative calculations or division to ensure that divisors and loop increments (like `$ywin` or `$ref`) are strictly positive.

## 2026-06-01 - Algorithmic Complexity DoS in Poisson Probability Calculation
**Vulnerability:** The `sv::poisson` function used `Math::BigFloat`'s `bfac()` (factorial) and `bpow()` (power) methods, which exhibit $O(N)$ or worse complexity. For large inputs (e.g., $x=1,000,000$), this caused extreme CPU usage and process hangs, leading to a Denial of Service.
**Learning:** Even when using arbitrary-precision libraries, certain operations like factorials can be computationally prohibitive. For statistical purposes, approximations (like Ramanujam's or Stirling's) often provide sufficient precision with $O(1)$ complexity.
**Prevention:** Implement threshold-based logic in statistical functions to switch to efficient approximations for large input values. Always validate and bound inputs to prevent expensive operations from being triggered by untrusted data.

## 2026-06-03 - Denial of Service via Infinite Loops in Structural Variation Selection
**Vulnerability:** The `sv::createsv` subroutine used `while` loops with `redo` to find non-overlapping genomic positions for multiple structural variations. In cases where the genome was small or saturated with mutations, these loops could become infinite, leading to a Denial of Service (CPU exhaustion).
**Learning:** "Trial and error" search algorithms for valid positions must have an upper bound on attempts. If a valid configuration cannot be found within a reasonable number of retries, it is safer to terminate with an error than to continue indefinitely.
**Prevention:** Implement explicit retry counters and limits (e.g., 1000 attempts) for all iterative search loops that depend on random chance or data-driven constraints to satisfy overlap requirements.

## 2026-06-05 - Robust SAM Header Parsing and Command-Line Validation
**Vulnerability:** The application was susceptible to incorrect genomic data processing if SAM header tags (`SN`, `LN`, `ID`, `VN`, `CL`) were not in a specific order. Additionally, missing validation for command-line arguments and `GetOptions` return values could lead to unexpected behavior or DoS when provided with malformed or out-of-range inputs.
**Learning:** Manual parsing of tab-separated fields in bioinformatics headers should always be position-independent to handle different software outputs. Relying on library-level type checking (e.g., `Getopt::Long`) is preferred over manual string parsing for numeric inputs, combined with explicit bounds checking.
**Prevention:**
1. Hardened SAM header parsing in `svre.pl` to iterate through all fields and extract tags using regex.
2. Implemented strict range validation and error handling for all critical numeric command-line parameters.

## 2026-06-08 - Division-by-Zero in RIC Proportions and Input Validation in Poisson
**Vulnerability:** The application was susceptible to Denial of Service (DoS) via script termination due to division-by-zero when calculating proportions of relative information content (RIC) for structural variations if the total entropy sum was zero. Additionally, the `sv::poisson` function lacked input validation, potentially leading to runtime errors with malformed data.
**Learning:** Even after hardening main binning loops, downstream reporting and typing logic can still be vulnerable to zero-valued aggregates derived from edge-case data. Library functions also require internal validation to handle untrusted or malformed inputs securely.
**Prevention:**
1. Always guard division operations in reporting logic with checks for non-zero denominators, especially when the divisor is an aggregate sum of data-dependent values like entropy.
2. Implement strict input validation (`defined` and `isfloat` checks) at the beginning of library subroutines to ensure they fail gracefully when provided with invalid data.

## 2026-06-10 - Denial of Service via Division-by-Zero in rcount processing
**Vulnerability:** The application crashed (DoS) when calculating probability distributions in genomic bins that contained only translocations, resulting in a standard read count (`rcount`) of zero.
**Learning:** Genomic bins can have zero standard reads if they only contain translocation data (stored under the "pair" key). Binning logic must not assume `rcount > 0` even if the bin exists in the data structure.
**Prevention:** Explicitly guard all divisions using `rcount` as a divisor in `sv::ric` and output generation loops.

## 2026-06-12 - Denial of Service via Missing Input Validation in Utility Functions
**Vulnerability:** The `sv::average` and `sv::stdev` subroutines in `sv.pm` were susceptible to script termination (DoS) when provided with empty arrays (via `die`) or non-numeric data (via fatal warnings/errors during arithmetic).
**Learning:** General-purpose utility functions in shared modules often lack the strict validation found in main entry points. When these utilities are used on data derived from untrusted or sparse genomic sources, they can cause unexpected application failure.
**Prevention:** Always harden utility functions that perform arithmetic or statistical operations. Filter input data for numeric types using `isfloat` or `looks_like_number` and return safe defaults (e.g., 0) instead of using `die` for empty input sets.

## 2026-06-14 - Bypass of Numeric Range Checks via NaN/Inf in Command-Line Arguments
**Vulnerability:** Command-line arguments (`bootstrap`, `fdr`, `ywindow`, `cov`, `mapq`) were vulnerable to `NaN` or `Inf` values bypassing numeric range checks. For example, `NaN <= 0` and `NaN > 100` both evaluate to false in Perl, allowing `NaN` to pass a `die "Error" if $x <= 0 or $x > 100` check.
**Learning:** In Perl, numeric comparisons with `NaN` always return false. This can lead to security bypasses if the application assumes that failing a range check implies the value is valid.
**Prevention:** Always validate that a numeric string is a valid number using a helper like `sv::isfloat()` (which uses a strict regex) before performing numeric range comparisons.

## 2025-06-20 - Metadata Corruption via Negative Coordinate Ranges
**Vulnerability:** Regex-based metadata extraction failed to handle negative numbers in genomic ranges, leading to coordinate data leaking into reference name variables.
**Learning:** In bioinformatics, coordinates can be negative (e.g., relative to a feature or on the reverse strand). Sanitization regexes must explicitly account for optional signs in all numeric parts of a range.
**Prevention:** Use robust regexes like `s/^-?\d+(\.\.-?\d+)?___//` when stripping coordinate prefixes that might contain negative values.

## 2026-06-16 - Algorithmic Complexity DoS in sv::bootstrap
**Vulnerability:** The `sv::bootstrap` function in `sv.pm` was susceptible to an Algorithmic Complexity DoS. For each bootstrapping iteration, it performed a linear scan of a slice of the Cumulative Distribution Function (CDF) array to map random samples to bins. When processing datasets with a large number of bins (e.g., 1,000,000), this led to $O(N)$ complexity per sample, causing extreme CPU usage and multi-minute execution times.
**Learning:** Sampling from a large discrete distribution using a linear CDF sweep is a common performance bottleneck and security risk. Sorting samples and sweeping can help, but direct binary search into the CDF is more robust and efficient for independent samples.
**Prevention:** Replace linear scans or slice-based sweeps in statistical sampling functions with $O(\log N)$ binary search. Ensure that the search logic correctly handles boundary conditions (like `rand()` returning 1.0) to maintain statistical integrity.

## 2026-06-18 - Denial of Service via Undefined Hash Dereferencing in ric
**Vulnerability:** The `sv::ric` function was susceptible to runtime crashes (DoS) when processing a reference genome where some chromosomes had no mapped reads. Attempting to call `keys` on an undefined hash reference (e.g., `%{$precount->{$ref}}`) caused the script to terminate.
**Learning:** Genomic data is often sparse, and some chromosomes in a reference might not have any observations in a specific dataset. Code that iterates over chromosomes based on a reference list must gracefully handle missing data for those chromosomes in data-dependent hashes.
**Prevention:** Use the defined-or operator (`// {}`) or explicit `defined` checks before dereferencing hash or array references that might be missing for specific genomic regions or chromosomes.

## 2026-05-23 - Logic Errors in Utility Functions and Redundant Qualification
**Vulnerability:** The `sv::sortu` utility function incorrectly included undefined elements in the sorting process, leading to potential warnings and non-deterministic results. The `sv::range` function lacked a definedness check for its range parameter, causing uninitialized value warnings.
**Learning:** General-purpose utility functions in shared modules often contain edge-case bugs that can manifest as runtime warnings or logic failures when processing sparse genomic data. Redundant qualification of internal functions with package prefixes (e.g., `sv::`) within the same package is non-idiomatic in Perl and should be avoided to maintain code clarity.
**Prevention:**
1. Implement strict definedness checks at the beginning of utility functions that accept optional or data-dependent parameters.
2. Ensure that filtering logic for sorting (e.g., creating `@defined` lists) correctly excludes `undef` values to prevent sorting errors.
3. Follow idiomatic coding standards for the language; in Perl, avoid redundant package prefixes for internal calls within the same namespace.

## 2026-06-21 - Path Traversal via Output Filenames
**Vulnerability:** User-provided output filename prefixes were susceptible to path traversal via `..` components, allowing arbitrary file writes outside the intended output directory.
**Learning:** `File::Spec->rel2abs()` resolves relative paths to absolute ones but preserves `..` components. If these absolute paths containing `..` are used in file operations or passed to external tools (like `Rscript`), they can still be exploited for path traversal.
**Prevention:** Use `File::Spec->canonpath()` to normalize the path and explicitly check for `..` components (e.g., using regex `/\.\.(\/|\\|$)/`) before proceeding with file operations.

## 2025-06-25 - R Command Injection via read.table and DoS via ARG_MAX
**Vulnerability:** R's `read.table` function can execute arbitrary shell commands if the provided filename starts with a pipe character (`|`). Additionally, passing thousands of chromosome names as command-line arguments to `Rscript` can trigger "Argument list too long" (E2BIG) errors.
**Learning:** External tools often have language-specific "features" like R's pipe-opening behavior that can lead to command injection even when the main application is otherwise secure. Command-line length limits are a real DoS vector for tools processing large genomic datasets with many scaffolds.
**Prevention:**
1. Explicitly block input strings that start with `|` when they will be used as filenames in R.
2. Use temporary files to pass large metadata collections between processes instead of relying on command-line arguments.
3. Use robust `read.table` parameters (`quote=""`, `comment.char=""`) to handle untrusted data within files.

## 2025-06-25 - Regex Line Anchor Bypass (Newline Injection) in Perl
**Vulnerability:** Use of the `$` anchor in validation regexes allowed strings with trailing newlines to pass as valid numeric or coordinate data. This could lead to log injection or logic bypasses when the data is later used.
**Learning:** In Perl, the `$` anchor matches at the end of the string OR before a newline at the end of the string. The `\z` anchor is the strict end-of-string anchor.
**Prevention:** Always use `\z` instead of `$` for strict end-of-string validation in Perl. Avoid `chomp` within validation functions if they are intended to be pure validators, as `chomp` modifies the input (if not used on a copy) and can mask malformed data.

## 2025-05-25 - Data Integrity and DoS via Negative Array Indexing in Bin Merging
**Vulnerability:** The `sv::ric` function in `sv.pm` was susceptible to data corruption and script crashes (DoS) when processing chromosomes with exactly one small genomic bin. Due to Perl's negative array indexing, a check for the "previous" bin (`$sortbin[$#sortbin - 1]`) would return the first bin itself if only one existed, causing the bin to merge with itself (doubling its count) and then delete itself.
**Learning:** In languages that support negative array indexing (like Perl or Python), boundary checks using offsets from the array length must explicitly verify that the array has sufficient elements to prevent unintentional self-reference or wrap-around behavior.
**Prevention:** Always use `scalar(@array) >= N` or explicit index bounds checks before accessing array elements using negative offsets or calculated indices near the boundaries.
