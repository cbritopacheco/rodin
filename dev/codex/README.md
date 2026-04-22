# Codex setup + maintenance scripts

This folder provides local scripts that mirror the CI workflow conventions in `.github/workflows/` so Codex (and humans) can bootstrap and maintain a consistent Rodin dev environment.

## Repository audit highlights

The scripts were derived from these workflow patterns:

- **Compiler/toolchain setup** from `.github/workflows/copilot-setup-steps.yml`.
- **Ubuntu dependency package set** from Build/Tests/Installation workflows.
- **macOS dependency package set** from Build/Installation workflows.
- **CMake flags and test labels** from `.github/workflows/Build.yml` and `.github/workflows/Tests.yml`.

## Scripts

### `setup_codex_env.sh`

Bootstraps a machine for Codex work:

- Initializes submodules.
- Optionally installs OS dependencies (`--install-deps`).
- Configures a dedicated build directory (`build/codex` by default) with unit/manufactured tests enabled.

Example:

```bash
bash dev/codex/setup_codex_env.sh --install-deps
```

### `maintenance.sh`

Runs periodic maintenance tasks aligned with CI:

- Submodule refresh.
- CMake reconfigure.
- Build.
- Unit/manufactured test labels (excluding `slow`).
- Optional installation smoke test (`--installation-check`).

Example:

```bash
bash dev/codex/maintenance.sh
```

## Notes

- These scripts are intentionally conservative and keep docs/examples/benchmarks disabled by default for faster local feedback.
- Override `BUILD_TYPE`, `CTEST_PARALLEL_LEVEL`, or `BUILD_DIR` via environment variables or flags where supported.
