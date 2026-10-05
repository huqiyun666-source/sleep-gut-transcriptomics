# Repository privacy and secret audit

Audit date: 2026-10-05  
Scope: every file under this repository package, including scripts, documentation, tables, provenance copies, and final figure PDFs where text extraction was applicable.

## Automated checks

| Check | Result | Notes |
|---|---|---|
| Local user-home, desktop, mounted-volume, and temporary-system absolute paths | PASS | No matches. |
| Email-address review | PASS | One corresponding-author email is intentionally included from author-provided publication metadata. No other email address was detected. |
| Credential patterns (API key, secret key, access/auth token, password assignment, cookies, private-key headers) | PASS | No matches. |
| Credential-bearing/private URL patterns | PASS | No token, key, session, authentication, signature, or embedded user-info URL parameters detected. |
| Raw sequencing and large binary analysis extensions | PASS | No FASTQ, SRA, BAM, SAM, RDS, or RData files included. |
| Shell history, SSH private keys, cookies, or hidden account files | PASS | None included. |
| macOS metadata | PASS | No `.DS_Store` files included. |

## Manual review

- Public accession, PubMed, DOI, GEO, SRA, BioProject, ENA, and journal URLs were retained because they are public provenance, not private endpoints.
- Names appearing in scholarly citations or dataset titles are publication metadata. No private contact names, phone numbers, or unpublished author identities were added.
- Author names use the English publication spellings confirmed by the authors. The corresponding-author email and postal information were supplied specifically for publication metadata.
- No suspected credential was found, so no secret value is reproduced in this report.

## Conclusion

PASS. All included files are classified as `public_safe = YES` in `file_manifest.csv`. This is a content audit, not a substitute for the authors' final pre-publication review.
