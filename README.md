# RF shimming for combined brain and cervical-spine MRI at 9.4 T

MATLAB code accompanying a manuscript on static RF shimming for combined brain and cervical-spine MRI at 9.4 T. The study investigates SAR-constrained RF shimming using simulated and measured transmit field maps, including subject-tailored and universal RF shims. The investigations are based on a 16-channel double-row coaxial-end dipole array previously introduced by [Solomakha et al., NMR in Biomedicine (2026)](https://doi.org/10.1002/nbm.70228).

## Reproducing the figures

| Script | Figure |
| --- | --- |
| `rf_shimming_simulations.m` | Figure 2: simulation-based RF shimming |
| `array_row_contributions.m` | Figure 3: contribution of individual array rows |
| `universal_rf_shimming.m` | Figure 4: universal RF shimming |

The code was developed using MATLAB R2022b. Optimization Toolbox and Image Processing Toolbox are required. Parallel Computing Toolbox is optional and accelerates the combinatorial calculation for Figure 4.

The public [`tdm_ptx`](https://github.com/fglang/tdm_ptx) repository is an external dependency. By default, it is expected as a sibling of this repository:

```text
parent/
├── head-neck-rf-shimming/
└── tdm_ptx/
```

Each script calls `setup_paths` automatically and can be run from any working directory. For a different checkout location, the path to `tdm_ptx` can be specified in the script's `setup_paths` call.

Input files are expected in `data/`. See [data/README.md](data/README.md) for details. Figure 4 stores its computationally expensive intermediate results in `results/US_all_combinations.mat` and reuses them on subsequent runs.

## Licenses and provenance

Original code in this repository is licensed under the MIT License. The included third-party packages `Panel` and `find_pareto_frontier` retain their respective BSD licenses. See `THIRD_PARTY_NOTICES.md` for details.