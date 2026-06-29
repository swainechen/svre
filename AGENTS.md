# Sentinel 🛡️ Instructions

## Genomic Coordinates
- **Never support scientific notation for genomic coordinates.** Genomic coordinates must always be treated as integers.
- Any attempt to introduce scientific notation support for coordinates should be rejected as it is non-standard in this domain and can lead to data integrity issues.
- Always use `RE_INT` for coordinate-related regex matches.
