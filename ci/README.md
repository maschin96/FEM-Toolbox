# Self-hosted verification

The workflows target the registered Linux X64 runner. Python 3.11–3.13 use
immutable Debian Bookworm container images; container UID/GID 1001 matches
`machine-1` so native and container jobs can share and clean their workspace.
Each Python job creates a fresh environment. Artifact installation tests use
fresh directories outside the persistent source checkout.

Run the runner as a system service, with its account in the Docker group.
The service must be restarted after adding group membership. A runner attached
to a closed SSH terminal can report online but fail before spawning its worker.
The host needs Docker and Python 3.11 or newer for the handbook orchestration.
PETSc and TeX retain their existing pinned native container toolchains.

This VM exposes a generic KVM CPU without x86-64-v2. Current NumPy wheels cannot
import there. The `runner-numpy` matrix therefore builds the latest resolvable
NumPy source for each supported Python, with `cpu-baseline=none`,
`cpu-dispatch=none` and NumPy's bundled linear-algebra fallback. It checks import,
matrix multiplication and eigenvalues and publishes wheel/source hashes and
build settings. This follows the [NumPy CPU build options](https://numpy.org/doc/stable/reference/simd/build-options.html).
Consumers verify the transferred hashes and install that exact wheel before
installing the project. Oldest-dependency jobs continue using the original
constraints and wheels. Versions are not downgraded to accommodate the runner.

These wheels are CI dependencies, never FEM release assets. Successful runs
demonstrate source-built latest dependency compatibility on this CPU; they do
not certify the upstream x86-64-v2 binaries on this VM or establish performance
parity with an optimized BLAS build. Numerical FEM tests and isolated package
installation gates remain enabled.
