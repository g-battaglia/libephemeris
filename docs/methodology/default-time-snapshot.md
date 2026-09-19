# Private default Delta-T byte snapshot

`libephemeris._time_snapshot` is a disconnected, private implementation of
the bounded contract in validation
`golden/specs/occultation_contact_residuals/default-time-three-asset-snapshot-contract.md`,
commit `19db2e17a803ce6a10975be8477caa28d152c10e` (SHA-256
`bc161b1f48c832afd3a6c464ab668b1aa338487f83f8adb22e590f62274846bd`).
It constructs an independent Timescale for the successful enhanced default
SMH-2016 Delta-T branch. It does not change the ordinary public time path.

## Source and ownership

The numerical model uses the published Delta-T work recorded under
`DELTA_T` in the provenance registry and the installed MIT-licensed
Skyfield 1.54 `Timescale`, `DeltaT` and spline primitives. The project
owns the file admission, array validation, explicit curve construction,
branch declaration and final native-float day conversion. NumPy 2.3.5
supplies the version-pinned array and interpolation operations. The
operation order is reviewed for numerical equivalence to the installed
Skyfield 1.54 builder; the code does not import an ambient bundled-data
loader into the private construction path.

The accepted production inputs are fixed independently of their paths:

| Role | Exact bytes | SHA-256 |
| --- | ---: | --- |
| `iers.npz` | 62,966 | `c7d7536d898dfa9f8cd43e8044ff51e108cc8289675a13fee9822010a1c4935c` |
| `historic_deltat.npy` | 10,576 | `f5346b780b36a0325b1847dc6c0083d66edc7e88b7f648b4c98a67bbd02b5d3f` |
| `delta_t.npz` | 1,547 | `2d12bd3e789543b78a1f53c8b76ed7fecffdf7e5149cfb6a0aed21a8b3db5ff6` |

Each descriptor is closed before parsing. The private factory checks exact
byte count and SHA-256 on copied immutable bytes, then parses only those
bytes with pickle loading disabled. The fixed schema and relational checks
reject malformed epochs, leap tables and spline intervals. A separate
synthetic factory exercises the same path but cannot issue a production
asset identity. A live production tag is bound to the still-owned byte
objects and is rechecked on retrieval; closing the evaluator drops it.
SHA-256 identity relies on the usual collision-resistance assumption.

For this first environment, the implementation also guards the exact
Skyfield and NumPy versions and the directly reviewed Skyfield source
files `timelib.py`, `curvelib.py` and `functions.py`. Their SHA-256 digests
are recorded in the validation contract. These guards do not constitute a
complete transitive dependency or runtime-code manifest.

## Computational boundary

The private evaluator accepts a native binary64 UT1 Julian date and an
explicit declaration of the default time configuration. It rejects an
override, alternate model, optional IERS branch or ordinary fallback. It
uses only its fresh Timescale, converts Delta-T seconds to a native Python
float and applies the project default day conversion. The initial private
numerical guard accepts only `|UT1 JD| <= 10^12`; any finite date outside
that bound, or one whose intermediate or final value becomes nonfinite, yields
a distinct typed unsupported result. The bound is a computational guard,
not a claim of physical validity over that span.
No numeric Delta-T output is stored in this documentation or its tests.

The declaration is caller supplied: it is not proof that a live ordinary
`deltat()` call selected this branch. Existing `_TS` and cached Skyfield
`Time` objects have no retroactive asset identity. A future call-local
receipt must capture the selected configuration and equality of the
ordinary and private results, then bind time to the exact LEB state graph
and both Sun/Moon returned words. This component alone certifies neither
an ordinary astronomical state nor physical accuracy, contact geometry,
public routing or release readiness.
