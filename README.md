<h1> <img src="http://www.combine.rwth-aachen.de/files/cbio/content/PhysioSpaceLogo2.png" width=100> PhysioSpace Methods</h1>
PhysioSpace is a robust statistical method for relating high dimensional omics data.^[Lenz, Michael, et al. "PhysioSpace: relating gene expression experiments from heterogeneous sources using shared physiological processes." PLoS One 8.10 (2013): e77627]. It is designed to take advantage of the vast availability of public omics data, which in combination with statistical approaches make a potent tool capable of analyzing heterogenious biological data sets.

PhysioSpaceMethods is a R package which provides an implementation of PhysioSpace method alongside other handy functions for making PhysioSpace an accessible tool for R users.



#### Table of Contents
**[Installation Instructions](#installation-instructions)**<br>
**[Usage Instructions](#usage-instructions)**<br>

### Installation Instructions
#### Installing via Devtools (Recommended method):
Easiest way to install PhysioSpaceMethods is via <a href="https://cran.r-project.org/web/packages/devtools/">Devtools</a>.
After installing Devtools from cran, you can install PhysioSpaceMethods by:
```r
devtools::install_github(repo = "JRC-COMBINE/PhysioSpaceMethods", build_vignettes = TRUE)
```

#### Alternative installation methods (Manual download):
In case you encountered any problem while installing PhysioSpaceMethods, you can download the repository first and 
install the package locally.
In your terminal, first clone the repository in your desired directory:
```Shell
cd [Your desired directory]
git clone https://github.com/JRC-COMBINE/PhysioSpaceMethods.git
```
Then install the downloaded package using <a href="https://cran.r-project.org/web/packages/devtools/">Devtools</a>:
```Shell
R -e "devtools::install_local('./PhysioSpaceMethods/', build_vignettes = TRUE)"
```

### Usage Instructions
PhysioSpaceMethods can map user samples inside a physiological space, calculated prior from a compendium of known samples. Explanation of how to use this package is provided in a vignette accessible by:
```r
browseVignettes(package = "PhysioSpaceMethods")
```

### Test Environments
The package was tested with R 3.5.1 on Windows 10, Mac OS X (10.12.6) and Linux (CentOS 7.4).

