# Reporting security issues

Report a suspected security vulnerability to the maintainer at
[kerykeion.astrology@gmail.com](mailto:kerykeion.astrology@gmail.com).
Include the affected package/version, operating system, a minimal reproduction
and whether untrusted files, network input or credentials are involved.
Do not include secrets, production data or files from another ephemeris
distribution. No response-time or support-lifetime guarantee is implied.

For ordinary scientific discrepancies, use the canonical
[issue tracker](https://github.com/g-battaglia/libephemeris/issues) with your own
reproduction inputs, units, clocks, model settings and data identities.

Treat acquired data inodes as immutable while calculation owners use them.
Use sealed network policy for calculations that must not acquire data, and
bounded worker admission when exposing synchronous calculations to a service.
Neither setting is a complete sandbox for hostile files or arbitrary code.
