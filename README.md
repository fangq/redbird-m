![](./doc/images/redbird_banner.png)

# Redbird - a model-based diffuse optical imaging toolbox

- Copyright: (C) Qianqian Fang (2005-2025) <q.fang at neu.edu>
- License: GNU Public License V3 or later
- Version: 0.5.0
- Github: https://github.com/fangq/redbird-m

Redbird is a compact, portable, yet feature-rich diffuse optical imaging (DOI) and diffuse optical tomography (DOT) 
toolbox for MATLAB and GNU Octave. It provides a fast and experimentally-validated diffusion equation forward solver
based on the finite-element (FE) method, as well as advanced non-linear image reconstruction algorithms.

Redbird is a result of over two-decades of active research in DOT and image reconstruction methods from Dr. Fang's
lab. It has been extensively used in over a dozen journal publications related to optical breast imaging, prior-guided
image reconstruction techniques, multi-modal imaging, and wide-field ultra-high-density DOT system development, 
including analyses of clinical measurements from over 400 human subjects. The forward solver has been exensively
validated against our in-house Monte Carlo (MC) solvers including MCX and mesh-based Monte Carlo (MMC). The inverse
solvers are packed with highly advanced reconstruction techniques, including dual-mesh modeling, multi-spectral
chromorphore estimation, adjoint method, structural prior-guided reconstrucitons, multi-right-hand-side (block)
iterative linear solvers, simultaneous frequency-domain (FD) and continuous wave (CW) data reconstructions, 
log-amplitude-phase reconstructions and many more, many of which were directly derived from Dr. Fang's PhD works.

