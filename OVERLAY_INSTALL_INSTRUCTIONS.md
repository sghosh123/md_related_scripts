# Applying this repository update

This ZIP contains the **complete merged repository** in a folder named `copper-trafficking-mechanism/`. It does not make changes to GitHub merely by being uploaded as a ZIP file.

## Recommended: apply to your existing local Git clone

1. Back up or commit your current working tree.
2. Extract the ZIP.
3. Copy the *contents* of the extracted `copper-trafficking-mechanism/` folder into the root of your existing local clone, allowing matching files to be overwritten. Keep the existing `.git/` directory; do not replace or delete it.
4. Inspect `git status` and `git diff`, then run:

   ```bash
   python -m pip install -e '.[test]'
   python -m pytest -q
   ```

5. Commit the changes on a branch and push that branch to GitHub. Review the diff before merging into the default branch.

The archive intentionally does not include the repository's `.git/` history. Never delete `.git/` to apply this update.

## If you only want to apply the overlay files

Use the separate `copper-trafficking-mechanism-v0.4.0a1-overlay.zip`, which contains only changed and new paths. Extract it into the root of your existing local clone with overwrite enabled, then inspect the diff and run the tests before committing.

## Release boundary

This build is an integration draft, not an assertion that all methods are validated on real trajectories. Check `docs/migration-log.md`, `docs/analysis-conventions.md`, and `LICENSE-STATUS.md`. The existing license status is unresolved; confirm rights and choose an appropriate license before making a public release.
