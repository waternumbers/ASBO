# anomalous

A native R implementation of the [PELT](https://doi.org/10.1080/01621459.2012.737745) and [CAPA](https://doi.org/10.1002/sam.11586) algorithms for changepoint and anomaly detection.

The implementation is build around an `S3` class `anomalous_partition` for storing inforamtion about the changepoints and anomalies and `R6` classes for the cost functions relating to fitting different distributions to the data. 
The CAPA and PELT implementations make use of the methods of the `anomalous_partition` class. This allows for the replacement of the current `anomalous_partition` class with a new class with similarly names methods.
The `anomalous_partition` class makes us of the cost classes via named methods, again allowing for generalisation and extension with further cost functions.

## Using the code

This package is not available on [CRAN](https://cran.r-project.org/). The
latest development version can be installed from
the git repository from within R using the devtools package: 

```
devtools::install_github("waternumbers/anomalous")
```

Prebuild packages of the development version available from the [r-universe](https://waternumbers.r-universe.dev/anomalous)

## Acknowledgements

Development of this code was supported by UK Research and Innovation (UKRI)
through the Engineering and Physical Sciences Research Council (EPSRC)
project "Reducing End Use Energy Demand in Commercial Settings Through
Digital Innovation" (EP/T025964/1).
