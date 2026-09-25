# Contributing to MEDS

MEDS is open source under the [Apache License 2.0](LICENSE).

## Your contribution and its license

- **Contributions are licensed under Apache 2.0.** When you submit a pull request, your
  contribution is licensed under the Apache License 2.0, as section 5 of that license provides.
  You keep the copyright. There is no contributor agreement to sign and no copyright transfer.
- **Submit only work you have the right to license.** Many employers own what their staff write
  at work or under a research grant. If yours might, check with them before contributing.
- **Add yourself to [`AUTHORS`](AUTHORS) in your first pull request.** If your employer holds
  the rights to your work, add your employer instead. `AUTHORS` is the list that the copyright
  line in [`NOTICE`](NOTICE) refers to.
- **Start each new source file with the license identifier.** Use
  `! SPDX-License-Identifier: Apache-2.0` in Fortran, and
  `# SPDX-License-Identifier: Apache-2.0` in Python, shell and CMake files. If the file starts
  with a `#!` line, put the identifier on the line after it.

## Pull requests

- **Branch from `beta`, and open pull requests against `beta`, not `main`.** `beta` collects
  merged work that is not yet in a release. It merges into `main` when a release is cut.
- A closing keyword such as `Fixes #N` in a pull request into `beta` does not close the issue.
  The issue is closed when the release reaches `main`.
- Record what you changed in [`CHANGELOG.md`](CHANGELOG.md), in the same pull request.

## Where things are

- Building MEDS and running the tests: [`docs/building.md`](docs/building.md).
- The source tree, and where a new file goes: [`src/README.md`](src/README.md).
- Coding and documentation conventions: [`CLAUDE.md`](CLAUDE.md).
