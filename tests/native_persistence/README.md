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

* `acquired_edges` retains legitimate same-argument coincidences after canonical
  endpoint merging when zero or one endpoint was originally shared; separated
  edges remain separate. Distinct open arcs already sharing both original
  endpoints do not acquire a spurious common block when those vertices
  canonicalize.
* `closed_edges` distinguishes an originally shared closing endpoint from
  distinct original endpoints that become coincident during the operation.
* `sliver_faces` uses independently generated narrow cylindrical faces whose
  closing boundary is within vertex tolerance. Perpendicular hatching supplies
  a valid interior point; same-domain comparison succeeds in both directions
  and rejects a translated face.
* `trimmed_hatches` verifies the perpendicular retry against a generated hole:
  the returned witness lies inside the actual trimmed domain, the hole is
  excluded, and an ordinary rectangular face retains its original behavior.

The new controls keep the Boolean fuzzy value at 1e-7 and retain all existing
hatcher-domain validity checks. No supplier model, importer repair, tolerance
relaxation or external framework is used.

The closed-edge control also uses initially open circular arcs sharing only one
original vertex. Ordinary VE processing merges their endpoints after both pave
ranges have been initialized. The unchanged long ranges must form a common
block through the existing coincidence comparison; a separated-tool control
retains open ranges and no common block. Read-only stage observations establish
the late transition without injecting paves or changing execution order.
