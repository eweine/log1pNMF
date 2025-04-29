## Installation Instructions

To install the package, run

```
remotes::install_github("eweine/passPCA")
```

## Using OpenMP

The package's main function, `fit_poisson_log1p_nmf`, has multithreading capabilities via `OpenMP`. 
On most Windows and Linux operating systems, `OpenMP` should be installed already. However, on newer versions
of MacOS, `OpenMP` needs to be installed manually. For instructions on how to set this up please refer to
`data.table`'s [helpful guide](https://github.com/Rdatatable/data.table/wiki/Installation#Enable-openmp-for-macos).
