# A concise representation for adaptive hexagonal meshes

This is the code accompanying the paper [A concise representation for adaptive hexagonal meshes](https://lhf.impa.br/ftp/papers/hex.pdf), submitted for publication (2024).

Adaptive hexagonal meshes were introduced in the paper [Hexagonal LOD for interactive terrain rendering](https://www.researchgate.net/publication/228673405) by Sussner et al. (2005).

Running `make` outputs a sample adaptive hexagonal mesh in eps and csv. The eps file is converted to pdf using pstopdf.

The file `in.csv` is a sample adaptive hexagonal mesh generated by random refinement. It can be loaded with `make I=in.csv`.

The code should work in both Python 2 and Python 3.

