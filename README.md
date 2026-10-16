## Installation Instructions

To install the package, run

```
remotes::install_github("eweine/log1pNMF")
```

## Using OpenMP

The package's main function, `fit_poisson_log1p_nmf`, has multithreading capabilities via [OpenMP](https://www.openmp.org/). 
On most Windows and Linux operating systems, `OpenMP` should be installed already. However, on newer versions
of MacOS, `OpenMP` needs to be installed manually. For instructions on how to set this up please refer to
`data.table`'s [helpful guide](https://github.com/Rdatatable/data.table/wiki/Installation#Enable-openmp-for-macos).

## Manuscript figures

For information on the code used to run the experiments and produce the
figures in the accompanying manuscript, see
[`inst/paper_figures`](inst/paper_figures) and the README in that
directory.
