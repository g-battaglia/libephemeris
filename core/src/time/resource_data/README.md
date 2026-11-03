# Independent resource clock data

These exact already-admitted Skyfield 1.54 numeric clock inputs are included
by the autonomous Rust Time curve. They are not ephemeris kernels. The source
identities are the same as src/tables/smh2016_tables.rs; no numerical word was
changed. Runtime NumPy only decodes containers and preserves caller operators;
curve construction/evaluation is Rust. Skyfield data is MIT licensed under the
existing third-party license admission.

[
  {
    "name": "iers.npz",
    "bytes": 62966,
    "sha256": "c7d7536d898dfa9f8cd43e8044ff51e108cc8289675a13fee9822010a1c4935c",
    "source": "/opt/homebrew/lib/python3.14/site-packages/skyfield/data/iers.npz"
  },
  {
    "name": "historic_deltat.npy",
    "bytes": 10576,
    "sha256": "f5346b780b36a0325b1847dc6c0083d66edc7e88b7f648b4c98a67bbd02b5d3f",
    "source": "/opt/homebrew/lib/python3.14/site-packages/skyfield/data/historic_deltat.npy"
  },
  {
    "name": "delta_t.npz",
    "bytes": 1547,
    "sha256": "2d12bd3e789543b78a1f53c8b76ed7fecffdf7e5149cfb6a0aed21a8b3db5ff6",
    "source": "/opt/homebrew/lib/python3.14/site-packages/skyfield/data/delta_t.npz"
  }
]
