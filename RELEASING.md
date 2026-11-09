# Releasing

Releases are made by hand. The public repository
(`AKubanski/rlds-paired`) has no GitHub Actions workflows, and the Zenodo
record of the package (concept DOI `10.5281/zenodo.20801099`) is maintained by
manual uploads. Every new version is therefore added with **New version** on
that record, which keeps all versions under one concept DOI.

## Every release

1. Update `CHANGELOG.md` with what changed.
2. Sync the version number into `rlds_paired/_version.py`, `.zenodo.json` and
   `CITATION.cff` (which also gets today's `date-released`):

   ```bash
   python scripts/bump_version.py 1.2.0
   ```

3. Install and run the tests:

   ```bash
   python -m pip install -e ".[dev,mcmc]" h5py
   python -m pytest tests/ -q
   ```

   With the Pantheon+ files in `RLDS_DATA_DIR` all tests run; without them the
   six tests that need these files are skipped.
4. Commit, push and tag with git. Uploading through the browser does not work
   for this repository: it accepts files up to 25 MB and at most 100 files at
   a time, and the MCMC chains in `validation/chains/` are larger.

   ```bash
   git add -A
   git commit -m "v1.2.0: <summary>"
   git push
   git tag v1.2.0
   git push origin v1.2.0
   ```

5. On Zenodo, open the rlds-paired record → **New version** → replace the
   file with the archive of the tag (GitHub → Tags → v1.2.0 → zip) → set
   Version to `1.2.0` and the publication date → take the description,
   keywords and related works from `.zenodo.json` → **Publish**. Note the new
   version DOI.
6. Optional: create a GitHub Release from the tag (Releases → Draft a new
   release) with the CHANGELOG entry as release notes.

Keep the Zenodo–GitHub integration switched **off** for this repository. If it
is on, publishing a GitHub Release makes Zenodo create a separate record with
a new concept DOI.

## Why the version is bumped before tagging

`.zenodo.json` and `CITATION.cff` are part of the archived tree. A stale
`version` there would label the new archive with the old version number, so
`bump_version.py` runs before the commit that is tagged. `CITATION.cff` cites
the concept DOI, which stays valid for every version.

## Workflows

Earlier versions of this file described a CI workflow and a release workflow
in `.github/workflows/` (tests, version bump, tag, GitHub Release, PyPI via
Trusted Publishing). These files are not part of the public repository, and
the procedure above does not use them.
