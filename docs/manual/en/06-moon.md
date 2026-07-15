# Chapter 6 — The Moon: nodes, apogee, perigee, and Lilith

## What you will learn

In this chapter, you will discover why the Moon is the most complex celestial body to calculate, what the lunar nodes are and why they are linked to eclipses, how apogee and perigee work (and what a "supermoon" is), and what is meant by Black Moon Lilith in astrology.

---

## 6.1 Why the Moon is special

The Moon is the closest celestial body to us and also the one with the most complex motion. This is why it deserves an entire chapter.

**It is very fast.** While the Sun moves by about 1° per day along the zodiac, the Moon travels about 12°–15° per day — it crosses an entire zodiac sign in just over two days. This means that the time of birth matters a lot for the Moon's position: in just 2 hours, the Moon moves by more than 1°.

**Its orbit is unstable.** The lunar orbit is an ellipse, but an ellipse that continuously changes shape, inclination, and orientation. The Sun "pulls" the Moon, altering its orbit in complex ways. Mathematicians have cataloged hundreds of lunar perturbations — it is the most studied three-body problem in the history of astronomy.

**Parallax matters.** For all other celestial bodies, the difference between the geocentric position (from the center of the Earth) and the topocentric position (from your location on the surface) is negligible — less than 1 arcsecond. For the Moon, this difference can reach almost 1°. This is why a solar eclipse is total in one place and partial in another just a few hundred kilometers away.

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 18.0)

# Position of the Moon — note the speed
pos, _ = ephem.calc_ut(jd, ephem.MOON, 0)

print(f"Moon: {pos[0]:.4f}° of ecliptic longitude")
print(f"Latitude: {pos[1]:+.4f}°")
print(f"Distance: {pos[2]:.6f} AU ({pos[2] * 149597870.7:.0f} km)")
print(f"Speed: {pos[3]:.4f}°/day")
```

```
Moon: 19.1832° of ecliptic longitude
Latitude: +0.3292°
Distance: 0.002405 AU (359781 km)
Speed: 14.9963°/day
```

---

## 6.2 The lunar nodes: where the Moon crosses the ecliptic

The Moon's orbit does not lie on the same plane as the ecliptic (the plane of the Earth's orbit around the Sun). It is inclined by about 5.1° relative to it. The two points where the lunar orbit intersects the ecliptic are called **nodes**.

Imagine the ecliptic as the floor of a room. The Moon's orbit is a tilted circle that crosses the floor at two points:

- The **ascending node** (or North Node) is the point where the Moon "rises" — it passes from below to above the ecliptic, from negative to positive latitude.
- The **descending node** (or South Node) is the opposite point — the Moon "descends" below the ecliptic. It is always 180° from the ascending node.

The nodes are not fixed: they slowly move backwards along the zodiac, completing a revolution in about **18.6 years**. This retrograde motion of the nodes is caused by the gravitational pull of the Sun on the lunar orbit.

### Mean node and true node

LibEphemeris offers two versions of the lunar nodes:

- The **mean node** (`MEAN_NODE`) is calculated using a regular mathematical polynomial. It moves in a uniform and predictable way — no oscillations, no jumps. It is like the perfect clock of nodal motion.

- The **true node** (`TRUE_NODE`) takes into account all the real perturbations. The true node oscillates back and forth around the mean position with an amplitude of about ±1.5°. In detailed ephemerides, the true node can briefly move in direct motion (forward), something the mean node never does.

In astrology, the choice between the mean and true node is a subject of debate. Many astrologers use the mean node because of its regularity; others prefer the true node because it reflects the astronomical reality.

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)

# Mean node
pos_mean, _ = ephem.calc_ut(jd, ephem.MEAN_NODE, 0)

# True node
pos_true, _ = ephem.calc_ut(jd, ephem.TRUE_NODE, 0)

diff = pos_true[0] - pos_mean[0]
print(f"Mean node: {pos_mean[0]:.4f}°")
print(f"True node:  {pos_true[0]:.4f}°")
print(f"Difference: {diff:.4f}° ({diff * 60:.1f} arcminutes)")

# Dedicated functions (accept JD in TT)
jd_tt = jd + ephem.deltat(jd)
lon_mean = ephem.calc_mean_lunar_node(jd_tt)
lon_true, lat_true, dist_true = ephem.calc_true_lunar_node(jd_tt)

print(f"\nMean node (dedicated function): {lon_mean:.4f}°")
print(f"True node  (dedicated function): {lon_true:.4f}°")
```

```
Mean node: 15.6610°
True node:  15.6269°
Difference: -0.0340° (-2.0 arcminutes)

Mean node (dedicated function): 15.6624°
True node  (dedicated function): 15.6269°
```

### Why the nodes matter: eclipses

Eclipses occur **only** when the Sun is near a lunar node. Why? A solar eclipse requires the Moon to pass in front of the Sun — but since the lunar orbit is tilted by 5.1°, the Moon usually passes above or below the Sun, missing it. Only when the new Moon occurs near a node is the Moon aligned enough with the ecliptic to cover the Sun (we will discuss this in detail in Chapter 9).

---

## 6.3 Apogee and perigee: the near and far Moon

The Moon's orbit is not a perfect circle but an **ellipse** — an oval shape. This means that the distance between the Earth and the Moon is constantly changing:

- At the **perigee** (closest point), the Moon is about 356,000 km from the Earth.
- At the **apogee** (farthest point), the distance rises to about 407,000 km.

The difference is remarkable: at perigee, the Moon appears about 14% larger and 30% brighter than at apogee. When a full Moon coincides with the perigee, the media calls it a **"supermoon"** — although the difference to the naked eye is difficult to perceive.

But there is a complication: due to the perturbations of the Sun, the position of the apogee and perigee **changes continuously**. The ellipse of the lunar orbit slowly rotates in space, completing a full revolution in about 8.85 years. Moreover, the ellipse itself changes shape from month to month.

Because of this, LibEphemeris offers three different versions of the apogee (and similarly for the perigee):

**Mean apogee** — a smooth mean curve built from ERFA/IERS Delaunay arguments
and conventional inclined-orbit geometry. It is widely used in
astrology under the name "Black Moon Lilith".

**Osculating apogee** — the instantaneous position of the apogee, calculated from the true orbit of the Moon at that exact moment. It includes all perturbations and can oscillate by 20°–30° relative to the mean position. "Osculating" comes from the Latin *osculari* (to kiss) — the osculating orbit is the ellipse that "kisses" the true trajectory at a given instant.

**Interpolated apogee** — a compatibility-oriented analytical curve formed
from lunar Delaunay arguments and periodic perturbation terms. It is not an
observed body or a claim of greater physical accuracy.

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)

# Mean apogee (via calc_ut)
mean, _ = ephem.calc_ut(jd, ephem.MEAN_APOG, 0)

# Osculating apogee (via calc_ut)
oscu, _ = ephem.calc_ut(jd, ephem.OSCU_APOG, 0)

# Interpolated apogee (dedicated function, requires JD in TT)
jd_tt = jd + ephem.deltat(jd)
lon_int, lat_int, dist_int = ephem.calc_interpolated_apogee(jd_tt)

print(f"Mean apogee:       {mean[0]:.4f}°")
print(f"Osculating apogee:   {oscu[0]:.4f}°")
print(f"Interpolated apogee: {lon_int:.4f}°")
print(f"Mean-osculating difference: {oscu[0] - mean[0]:.1f}°")
```

---

## 6.4 Lilith: the Black Moon

In astrology, **Black Moon Lilith** (or simply "Lilith") is the apogee of the lunar orbit — the point where the Moon is farthest from the Earth. It is not a physical body: it is a geometric point in space, the empty focus of the Moon's orbital ellipse.

The name comes from mythology: Lilith is the first wife of Adam according to some Jewish traditions, associated with the dark, wild, and untamable side of feminine nature. In astrology, it represents deep instincts, what is hidden or repressed.

There are three versions of Lilith, which correspond to the three versions of the apogee described above:

**Mean Lilith** is the most widely used in astrology. It corresponds to the mean apogee and moves in a regular manner, completing a revolution of the zodiac in about 8.85 years (just under 9 years). It is the version you find in most astrological software and printed ephemerides.

**True (osculating) Lilith** is the true instantaneous position of the apogee. It can differ from the mean Lilith by 20°–30° and has irregular movements, including brief retrogradations. Some astrologers prefer it for its adherence to astronomical reality.

**Interpolated Lilith** is the versioned compatibility curve built from lunar
Delaunay arguments and a reviewed, hash-pinned refinement. For physical
Moon-distance events, calculate the Moon's distance directly.

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)
jd_tt = jd + ephem.deltat(jd)

# Mean Lilith (= mean apogee)
mean_lilith = ephem.calc_mean_lilith(jd_tt)

# True/osculating Lilith
true_lon, true_lat, true_dist = ephem.calc_true_lilith(jd_tt)

signs = ["Ari", "Tau", "Gem", "Cnc", "Leo", "Vir",
         "Lib", "Sco", "Sgr", "Cap", "Aqr", "Psc"]

def sign(lon):
    return f"{lon % 30:.1f}° {signs[int(lon / 30)]}"

print(f"Mean Lilith: {sign(mean_lilith)}")
print(f"True Lilith:  {sign(true_lon)}")
print(f"Difference:   {true_lon - mean_lilith:.1f}°")
```

### The White Moon (Selena)

Some astrological traditions use a distinct point called the **White Moon** or
**Selena**. ID 56 implements the uniform seven-year convention stated by
Velichko and Larin (2007), pp. 17–18. Its uniform rate and phase are unwrapped
from their January 1800 and January 2000 rows (pp. 20 and 45), then checked
against three unused rows through 2007. LibEphemeris declares the otherwise
unstated time convention explicitly: first day of the month at 00:00 TT. It is
a symbolic longitude, not a physical lunar apsis or satellite; the returned
radius is an IAU-standards-derived display convention.

```python
selena = ephem.calc_hypothetical_position(ephem.WHITE_MOON, jd_tt)
print(f"Selena longitude: {selena[0]:.6f}°")
```

The complete page-level derivation and arithmetic are in the
[hypothetical-body provenance record](../../methodology/hypothetical-bodies.md#id-56-white-moon--selena).

Use the documented mean, osculating, or interpolated lunar-apsis APIs when an
astronomical apsis is required.

---

## 6.5 Interpolated perigee

Just as there is an interpolated apogee, there is also an **interpolated
perigee**: a smoothed coordinate defined by LibEphemeris. Do not use
this abstract point as a shortcut for tides or "supermoon" timing. For those
questions, search the actual geocentric Moon-distance curve for a minimum.

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 8, 12.0)
jd_tt = jd + ephem.deltat(jd)

# Interpolated perigee
lon, lat, dist = ephem.calc_interpolated_perigee(jd_tt)

# Distance in km (1 AU = 149597870.7 km)
dist_km = dist * 149597870.7

print(f"Interpolated perigee: {lon:.4f}°")
print(f"Distance: {dist_km:.0f} km")
```

LibEphemeris evaluates interpolated apogee and perigee with separate analytical
Delaunay series and a versioned compatibility refinement. The refinement is
hash-pinned and cannot be modified without an explicit review. The returned
distances are conventional radii, not measurements of the Moon's instantaneous
distance. See
[Lunar Nodes and Apsides](../../methodology/lunar-apsides.md).

---

## 6.6 Node crossings

In certain cases, you need to know **exactly when** the Moon crosses the plane of the ecliptic — that is, when its ecliptic latitude passes through zero. This is important for calculating eclipses and for some astrological techniques.

The `mooncross_node_ut` function finds the next time the Moon crosses the ecliptic, starting from an initial date:

```python
import libephemeris as ephem

jd = ephem.julday(2024, 4, 1, 0.0)

# Find the next node crossing
jd_cross, lon_cross, lat_cross = ephem.mooncross_node_ut(jd)

year, month, day, hours = ephem.revjul(jd_cross)
h = int(hours)
m = int((hours - h) * 60)

print(f"Next crossing: {year}-{month:02d}-{day:02d} {h:02d}:{m:02d} UT")
print(f"Longitude at node: {lon_cross:.4f}°")

# lat_cross is ~0 by definition (it is the moment
# when the Moon is exactly on the ecliptic)
```

```
Next crossing: 2024-04-08 12:18 UT
Longitude at node: 15.6269°
```

The Moon crosses a node about **twice a month** — once at the ascending node (latitude from negative to positive) and once at the descending node (from positive to negative). You can distinguish the two cases by the sign of the velocity in latitude at the time of crossing: positive = ascending node, negative = descending node.

---

## Summary

- The Moon is the most complex celestial body to calculate: it is very fast (~13°/day), its orbit is continuously perturbed by the Sun, and the parallax can reach almost 1°.
- The **lunar nodes** are the points where the lunar orbit intersects the ecliptic. They exist in mean (regular) and true (with oscillations) versions. Eclipses occur only near the nodes.
- **Apogee** (farthest point) and **perigee** (closest point) exist in three versions: mean, osculating, and interpolated. The mean apogee is the "Black Moon Lilith" of astrology.
- **Lilith** (Black Moon) is the apogee of the lunar orbit — a geometric point, not a physical body. Mean Lilith is the most used; the true version can differ by 20°–30°.
- `mooncross_node_ut` finds the exact moment when the Moon crosses the ecliptic.

### Functions introduced

- `calc_mean_lunar_node(jd_tt)` — longitude of the mean lunar node
- `calc_true_lunar_node(jd_tt)` — longitude, latitude, and distance of the true lunar node
- `calc_mean_lilith(jd_tt)` — longitude of mean Lilith (mean apogee)
- `calc_true_lilith(jd_tt)` — longitude, latitude, and distance of true Lilith (osculating apogee)
- `calc_interpolated_apogee(jd_tt)` — interpolated (smoothed) lunar apogee
- `calc_interpolated_perigee(jd_tt)` — interpolated lunar perigee
- `mooncross_node_ut(jd_ut)` — next lunar node crossing
- `MEAN_NODE`, `TRUE_NODE` — lunar nodes via `calc_ut`
- `MEAN_APOG`, `OSCU_APOG` — lunar apogee via `calc_ut`
