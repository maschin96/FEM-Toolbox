# GitHub-hosted verification

Since the repository became public on 19 October 2026, both workflows use
GitHub-hosted `ubuntu-24.04` runners. Each job receives an ephemeral machine;
untrusted pull requests never execute on the persistent self-hosted runner.
Python 3.11–3.13 retain the pinned Debian Bookworm container images. Each Python
job creates a fresh environment. Artifact checks install transferred wheel/sdist
archives outside the checkout and verify their checksums, metadata and examples.

The test and latest-dependency jobs use upstream NumPy wheels on the hosted CPU.
Oldest-dependency jobs retain `oldest-constraints.txt`. The CPU-specific NumPy
source build is no longer a prerequisite; `build_runner_numpy.py` remains an
explicit utility for the historical generic KVM runner and is not run by CI.
PETSc real/complex and the five handbook PDFs retain pinned Docker toolchains.
The handbook workflow checks negative fixtures and source-archive reproduction;
on main it additionally requires clean release provenance with `--release`.

All actions remain pinned to full commit SHAs and the token has read-only
contents permission. The public PR event uses no secrets or privileged token.
