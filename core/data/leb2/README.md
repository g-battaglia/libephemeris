# Bundled base coefficients

`base_core.leb2` is the project-native coefficient image distributed with
LibEphemeris v3.2.2 / data-v4 at commit `7c77b70717cf8420b85c19d3aaad128060998007`,
path `libephemeris/data/leb2/base_core.leb2`.

Size: 10,232,271 bytes. SHA-256:
`b891de02f82b68346046da67389fd3713aa81ab299287c632cb20cc396397964`.
Stored core coverage: JD 2396757.5–2506332.5, one additional day at each end
compared with data-v3. The unchanged [upstream manifest](data-v4-manifest.json)
records identities and source-bounded intervals for all twelve data-v4 LEB2
assets; only this base core is embedded.

The optional core feature `bundled-data`, enabled by the Python package,
provides the baseline local calculation path. Its nominal date range is
1850–2150; the reader enforces each channel's actual coverage. Explicit
files and discovered tier files take precedence. Wider DE/SPK and LEB data
remain external. See the project licence and `THIRD_PARTY_NOTICES.md`.
