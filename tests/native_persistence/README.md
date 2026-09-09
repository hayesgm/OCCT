# Native persistence and projection conformance

These standalone C++ tests exercise the installed OCCT libraries with generated
binary records and mathematical geometry. No external CAD fixtures are needed.
Build the library normally, then configure this directory against that exact
installation. The tests are independent of DRAW and do not alter the OCCT build.

```sh
cmake -S tests/native_persistence -B /absolute/test-build -DOCCT_ROOT=/absolute/occt-install
cmake --build /absolute/test-build
ctest --test-dir /absolute/test-build --output-on-failure
```

The selected installation must provide the stored analytic and rigid-transform
reader capabilities. The build uses only its explicit include/library paths and
normal runtime search path. It must match this checkout's production source.

* `stored_transforms` checks exact elementary bits, both location readers,
  reconstruction controls, malformed/truncated input, stream scopes, rounding,
  reflection/inverse/power, compound/reference structure and finite arithmetic.
* `analytic_readers` checks line/circle/plane/sphere binary table readers,
  repeated exact payloads, reconstruction and malformed-input controls.
* `projection_controls` checks generated classifier projection cases, normal
  polynomial values/derivatives, nonsingular/point/periodic projected curves,
  and the two one-sided singular-boundary normals.

The existing interior singular-normal refusal is a separate known defect in
Release builds with `No_Exception`: the fallback can return non-finite values.
Run `projection_controls --interior-normal` to exercise the unchanged refusal
assertion separately. It fails on that configuration and is not classified as a
passing control. This change does not repair it or alter exception build flags.
The upstream reflection `gp_Trsf::Power(3)` issue is also outside this change;
the supported inverse/power controls here do not claim general Power correctness.
