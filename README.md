# Benchmark Paper
[![DOI](https://zenodo.org/badge/674374498.svg)](https://zenodo.org/doi/10.5281/zenodo.12607316)

This repository contains the code and files to reproduce the analyses presented in our paper: 

Luo, S., Germain, PL., Robinson, M.D. et al. Benchmarking computational methods for single-cell chromatin data analysis. Genome Biol 25, 225 (2024). https://doi.org/10.1186/s13059-024-03356-x

The benchmarks compare multiple computational workflows for single-cell
chromatin accessibility data, evaluate performance at different processing
stages, and provide guidelines for method selection.

## Usage
For downloading and processing all the datasets we used in our benchmark, see `./data`.

For generating the visualizations in our benchmark paper, see `./analysis` and `./result_files`.

For the reusable snakemake pipeline, see the repository: [sc_chromatin_benchmark](https://github.com/RoseYuan/sc_chromatin_benchmark).


## Citation
If you use this code or find it useful, please cite:
```bibtex
@article{Luo2024Benchmark,
  title   = {Benchmarking computational methods for single-cell chromatin data analysis},
  author  = {Luo, Siyuan and Germain, Pierre-Luc and Robinson, Mark D. and von Meyenn, Ferdinand},
  journal = {Genome Biology},
  volume  = {25},
  pages   = {225},
  year    = {2024},
  doi     = {10.1186/s13059-024-03356-x}
}
```
