# Test on 2D and realspace

This is a graphene for testing 2D disentanglement and real space WF

- I put atoms in the middle of unit cell, to test that the WF center is correctly computed
- the `unk` files are reduced by 3, as shown in `p2w.in`, to save space
- the `A` matrices from Wannier90 `chk` are written in `graphene.w90.amn`
- I didn't use Wannier90 output xsf, because it is shifted by 1 grid point,
    so I use WJL output xsf as reference
