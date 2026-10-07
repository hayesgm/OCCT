// Binary64 plane frame with both rounded-zero and subnormal squared components.
// All serialized coefficients are normal or zero; retain every original bit.
static const double THE_UNDERFLOW_FRAME[] = {
  0x1.3333333333333p-3, 0x1.4325c8f0bf238p+1, -0x1.f3ffffffffffdp-108,
  -0x1.a7e810355e07fp-313, 0x1.c37186d10a9a5p-209, 0x1p+0,
  0x1.30c90d827725p-578, -0x1p+0, 0x1.c37186d10a9a5p-209,
  0x1p+0, 0x1.75c4e89d655b2p-521, 0x1.a7e810355e07fp-313
};
