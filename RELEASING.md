# Releasing

Releasing is **one manual step**. You never edit version numbers by hand.

## One-time setup (do once)

1. **PyPI Trusted Publisher** — on pypi.org, project `rlds-paired` →
   *Publishing* → add a GitHub publisher:
   - Repository: `Askuba77/rlds-paired`
   - Workflow: `release.yml`
   - Environment: `pypi`
2. **Zenodo–GitHub** — log in to Zenodo with GitHub, flip the
   `Askuba77/rlds-paired` repository switch **on**. Zenodo will archive every
   published GitHub Release and mint a versioned DOI.
3. After the **first** release, copy the software DOI Zenodo mints into
   `CITATION.cff` (uncomment the `doi:` line) and the README, then commit.

## Every release

1. Update `CHANGELOG.md` with what changed.
2. GitHub → **Actions → Release → Run workflow** → type the version (e.g. `1.0.1`).

That's it. The workflow then:

1. runs the full test suite — nothing irreversible happens unless it passes;
2. runs `scripts/bump_version.py`, syncing the version into
   `rlds_paired/_version.py`, `.zenodo.json` and `CITATION.cff`
   (plus today's `date-released`);
3. builds the wheel + sdist from the bumped working tree;
4. commits the version bump and pushes it;
5. creates tag `v<version>` and a GitHub Release — **Zenodo archives the
   tagged tree, so the Zenodo record carries the correct version**;
6. publishes the already-built artifacts to PyPI via Trusted Publishing
   (OIDC, no stored token).

## Why a script + tag, not just a tag

Zenodo reads `.zenodo.json` from the **committed tree of the release tag**. A
stale `version` there makes every Zenodo record share one version label. The
workflow bumps and commits *before* tagging, so the archived tree is always
correct — and you only ever type the number once.

## Manual fallback

You can also bump locally without the workflow:

```bash
python scripts/bump_version.py 1.0.1
```

then commit, tag `v1.0.1`, and create the GitHub Release yourself.
