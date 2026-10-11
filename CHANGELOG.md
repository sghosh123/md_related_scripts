# Changelog

## 0.4.0a1 (integration draft)

- Integrated an interface feature extractor for eight Cα distances and four orientation angles into the existing `copper_trafficking` namespace.
- Added shared periodic-boundary geometry helpers and explicit feature-selection validation.
- Added streaming interfacial heavy-atom contact extraction with trajectory and segment provenance.
- Extended multi-distance analysis with optional cutoff-based presence flags and provenance IDs.
- Added synthetic regression/CLI tests and updated scientific conventions, examples, workflow catalog, and migration notes.
- This is a pre-release integration draft. Real-system regression, source/data licensing, and final citation metadata remain outstanding.

## 0.3.0

See `docs/migration-log.md` for the prior simulation-template and audit-correction release notes.

## 0.4.0a2

- Added explicit salt-bridge/charged-group distance screens with selection validation and per-frame provenance.
- Added complete frame-manifest export to streaming interface-contact analysis.
- Added replicate-aware residue-pair contact consistency summaries that include zero-contact frames and pool file segments within independent trajectory IDs.
