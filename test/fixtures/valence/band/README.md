# Steps to run w90

The `chk.fmt` file is generated with

```bash
wannier90.x silicon
w90chk2chk.x -export silicon
```

The unitary matrices in `chk` is almost the same as `.amn`, since I set `num_iter = 0; dis_num_iter = 0`.
