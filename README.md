<h1> <img 
src="http://www.combine.rwth-aachen.de/files/cbio/content/PhysioSpaceLogo2.png" 
width=100> PhysioSpace Methods</h1>
PhysioSpace is a robust statistical method for relating high dimensional 
omics data.^[Lenz, Michael, et al. "PhysioSpace: relating gene expression 
experiments from heterogeneous sources using shared physiological processes."
PLoS One 8.10 (2013): e77627]. It is designed to take advantage of the vast
availability of public omics data, which in combination with statistical
approaches make a potent tool capable of analyzing heterogenious biological
data sets.

PhysioSpaceMethods is a R package which provides an implementation of 
PhysioSpace method alongside other handy functions for making PhysioSpace an 
accessible tool for R users.


### Installation Instructions
You can install this package by:
```r
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_gitlab("jrc-combine/PhysioSpaceMethods", host = "git.rwth-aachen.de", build_manual = TRUE, build_vignettes = TRUE)
```

### Usage Instructions
PhysioSpaceMethods can map user samples inside a physiological space, 
calculated prior from a compendium of known samples. Explanation of how to use 
this package is provided in a vignette accessible by:
```r
browseVignettes(package = "PhysioSpaceMethods")
```

### Test Environments
The package was tested with R 3.5.1 on Windows 10, Mac OS X (10.12.6) 
and Linux (CentOS 7.4).

