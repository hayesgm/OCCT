# Rational-spline roundtrip regression

The two BREP solids are a reduction of an ordinary Trinket Product's warm
geometry and reopened STEP export. Each has 44 faces; a same-basis rational
B-spline corner at diagnostic face index 8 reproduces an order-dependent
intersection. Baseline forward subtraction is empty, while reverse subtraction
leaves a valid solid with 690 occurrence entities. This is diagnostic fixture
geometry, not a replacement or accepted Product artifact.

The regression requires empty directional cuts of both complete solids and the
corner pair. It preserves genuine planar crossings, a spherical-pole crossing,
a cone-apex crossing, distinct regular trims of a looping spline, and actual
face self-intersection. Baseline passes the controls but fails the original
corner-pair difference; the corrected walker must pass both.

At the collapsed corner, nearly zero first derivatives produce unstable normal
directions. The existing local tangency check can consequently seed spurious
walking lines. The correction uses CSLib second-derivative normal recovery when
a tangent magnitude is below the existing tolerance; unavailable or ambiguous
higher derivatives leave the check inconclusive. The existing neighborhood
projection checks and distance/angle thresholds remain unchanged. No global
surface-coincidence shortcut is used.

SHA256:
- warm.brep: e8373907add0301a75b7e08af37cb8d2c98243d836a3bde71b7a09033d6f9957
- reopened.brep: 4a13568e99553b0fd73400b7afd4488e362aa6c19b332de8a7720f70414a5a0b
