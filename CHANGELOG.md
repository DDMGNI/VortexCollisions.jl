# Release Notes

All notable changes to VortexCollisions.jl.

This package is pre-1.0, so *every* minor release is potentially breaking in the sense of
[SemVer](https://semver.org) for `0.x` versions. The sections below name what actually
changed, so that a compat-only bump can be told apart from a rename or a change in results.

This file was started on 2026-08-31 and deliberately holds no entries. Nothing has been
released yet — there are no tags, and `Project.toml` stands at `0.1.0` — so the development
history that predates this file is in `git log` alone, reaching back to 2017. It is named as a
gap rather than reconstructed, because a changelog assembled after the fact loses exactly the
reasoning that makes it worth keeping.

## [Unreleased] — targeting 0.1.0

### New Features

- **Continuous integration.** The repository had no `.github` directory at all; it now carries the
  same `CI.yml`, `CompatHelper.yml`, `Documenter.yml` and `TagBot.yml` as every other repository in
  the tree — the test suite over two Julia versions and three operating systems, a doctest job, and
  coverage reporting.
- **A documentation build.** `docs/` previously held one untracked PDF; it now has `make.jl`,
  `Project.toml` and a two-page manual, so `Documenter` has something to build and to deploy.

### Bug Fixes

### Breaking Changes

- **A `[compat]` section, where there was none.** `julia = "1.10"` — the LTS and the floor across
  the tree, and the field the CI matrix resolves its lower entry from — plus bounds for the four
  non-stdlib dependencies: `AbstractFFTs = "1"`, `FFTW = "1"`, `HDF5 = "0.17"`,
  `ProgressMeter = "1"`, matching what the manifest resolves today. An unbounded dependency is free
  to break the package on its next breaking release.

## Open Issues
