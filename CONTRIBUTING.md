# Contributing to LibEphemeris

Thank you for your interest in contributing to LibEphemeris! Contributions of
all kinds are welcome: bug reports, feature requests, documentation
improvements, and code changes.

## Getting Started

1. Fork the repository and clone it locally:

   ```bash
   git clone https://github.com/<your-username>/libephemeris.git
   cd libephemeris
   ```

2. Install the development dependencies (requires [uv](https://docs.astral.sh/uv/)):

   ```bash
   uv pip install -e ".[dev]"
   ```

3. Create a new branch for your changes:

   ```bash
   git checkout -b my-feature
   ```

4. Make your changes, add tests where applicable, and run the targeted suites.
   Do **not** run the full suite during development — it is very long. Use the
   fast task suites instead:

   ```bash
   poe test:leb:fast                  # Recommended fast unit suite (~1 min)
   poe test:skyfield:fast             # Skyfield backend unit suite
   poe test:leb:core                  # Essential subset (~20 s)
   pytest tests/test_file.py -v       # A single file while iterating
   ```

5. Run the linter, the formatter, the type checker, and the provenance gates:

   ```bash
   poe lint
   poe format
   poe typecheck
   uv run python scripts/check_provenance.py
   uv run python scripts/check_algorithm_provenance.py
   ```

6. Push your branch and open a Pull Request against the `main` branch.

## Reporting Issues

- Use the [GitHub Issues](https://github.com/g-battaglia/libephemeris/issues) tracker.
- Include a clear description of the problem, steps to reproduce, and the
  expected vs. actual behaviour.
- If relevant, include the Python version, OS, LibEphemeris version, calculation
  mode (`auto` / `skyfield` / `leb` / `horizons`), and the data tier in use.

## Code Style

- `from __future__ import annotations` at the top of every module.
- Imports grouped: stdlib, third-party, local (relative imports).
- Line length 88, Python 3.12+, double quotes, [Ruff](https://docs.astral.sh/ruff/)
  formatter and linter.
- Google-style docstrings with `Args` / `Returns` / `Raises`.
- Naming: `snake_case` functions, `PascalCase` classes,
  `SCREAMING_SNAKE_CASE` constants, `_underscore` private.
- Always return native Python floats, never NumPy scalars.
- The public API keeps 1:1 signature compatibility with the reference
  ephemeris API; changing a public signature or a flag semantic is a breaking
  change and needs to be discussed in an issue first.

## Provenance Requirements

LibEphemeris is a scientific library, and every number it produces must be
traceable to a public source. For any new or materially changed numerical
implementation:

- Add an in-code `Provenance` section and a matching record in
  `docs/methodology/provenance-registry.toml`.
- Cite the defining public source (JPL/IAU documentation, peer-reviewed
  literature, or primary publication), and state units, reference frame, time
  scale, and validity limits.
- Keep source values separate from project choices, and identify the generator
  for any checked-in data file.
- Agreement with a compatibility target is **never** an acceptable source for a
  value.

### Clean-room policy

This project is developed independently of the reference implementation. By
contributing, you confirm that you have not inspected, retrieved, translated,
adapted, or copied the reference implementation's source code, source comments,
documentation prose, algorithms, or data files while preparing your
contribution, and that no reference-distribution file is part of it. The
reference API may be used only for ephemeral behavioural comparison — call its
public API and compare outputs — and those outputs must never be persisted,
fitted, or committed to this repository.

## License

This project is licensed under the **GNU Affero General Public License v3.0**
(`AGPL-3.0-only`). See the [LICENSE](LICENSE) file for the full text and
[LICENSING.md](LICENSING.md) for details.

## Contributor License Agreement (CLA) — Copyright Assignment

By submitting a pull request or any other contribution to this repository, you
agree to the following terms in respect of every contribution you have made and
will make to the project:

1. **Copyright Assignment.** You hereby irrevocably assign to the project
   maintainer, **Giacomo Battaglia** ("Maintainer"), all right, title, and
   interest in the copyright of your past, present, and future contributions to
   this repository. To the extent any such assignment is not effective under
   applicable law, you instead grant the Maintainer a worldwide, perpetual,
   irrevocable, royalty-free, exclusive licence to use, reproduce, modify,
   distribute, sublicense, and re-license those contributions, with the right to
   grant sublicenses through multiple tiers. The Maintainer in turn grants you a
   perpetual, worldwide, non-exclusive, royalty-free licence to use and
   re-license your own contribution for any purpose, so this assignment does not
   stop you from reusing your own work.

2. **Re-licensing.** The Maintainer may re-license the project — or any part of
   it, including your contribution — under any other license, whether
   open-source or proprietary, at their sole discretion.

3. **AGPL Availability.** The project will continue to be publicly available
   under the AGPL-3.0 license. The copyright assignment enables dual-licensing
   and commercial offerings that help sustain long-term development.

4. **Moral rights.** To the fullest extent permitted by applicable law, you
   agree not to assert or enforce, against the Maintainer or its licensees, any
   moral rights you may hold in your contributions; where such rights are
   waivable, you waive them. Some jurisdictions (including Italy) treat certain
   moral rights as inalienable; this clause applies only as far as the law
   allows.

5. **Attribution.** Your authorship is acknowledged in the Git commit history,
   in the [AUTHORS](AUTHORS) file, and, where appropriate, in release notes.
   Copyright assignment does not erase your credit as the original author of
   your contribution.

6. **Originality & authority.** You represent that each contribution is your
   original work, that you have the right to assign its copyright (and, where
   you contribute in the course of employment, that you have your employer's
   authorization to do so), and that the contribution does not knowingly
   infringe any third party's rights. If any part of your contribution is
   subject to a third-party license, you must clearly state this in the pull
   request.

7. **Governing law.** This agreement is governed by the laws of Italy, without
   regard to its conflict-of-laws rules, and the courts of the Maintainer's
   place of residence shall have jurisdiction, without prejudice to any
   mandatory protections available to you under your local law.

Your agreement is given by your own act of submitting a contribution (as stated
above). On GitHub, a CLA-bot additionally checks each pull request and comments
with instructions if your account is not yet recorded as having agreed; this
automated check runs on GitHub only (it does not run on any mirror).

## Questions?

If you have any questions about contributing or the CLA, feel free to reach out
at [kerykeion.astrology@gmail.com](mailto:kerykeion.astrology@gmail.com?subject=Contributing%20to%20LibEphemeris).
