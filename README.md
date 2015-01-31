## Simplifying Trip Data into Networks via Coupled Matrix Factorization
This directory contains MATLAB code for summarizing spatio-temporal trip records into a network of movements based on these trips' starting location/time and destination location/time.
The code builds upon and utilizes routines from existing public software packages as listed below, and it runs on MATLAB version 7 or higher.

A typical use case consists of, in the given order, "preprocess", "convert2matrices" and "couplenmf". Some postprocessing is necessary to obtain the exact output (e.g. figures) from the paper but this step is more subject to how the user prefers the visual presentation.

File(s) in this directory: 

+ extract_trips.m: routine used in "preprocess.m" to extract trip record from a single file (e.g. one taxi in the data)
+ preprocess.m: extract trip records from raw GPS traces from the data set below 
+ latlong2ind.m: convert latitude and longitude coordinates to our spatial index; used in "convert2matrices.m"
+ ind2latlong.m: convert our spatial index format to latitude and longitude coordinates; used in "convert2matrices.m"
+ convert2matrices.m: transform the trip records resulting from "preprocess.m" to the matrices containing spatial temporal pickup and dropoff information to be decomposed
+ couplenmf.m: the major algorithm that performs coupled nonnegative matrix factorization to obtain the latent factors

Dependence: 

+ [NNLS solver](http://www.cc.gatech.edu/~hpark/nmfsoftware.php): a fast solver for nonnegative least squares problem implemented in MATLAB by Park & Kim. It is free for download at the authors' website.

The methods (including preprocessing and postprocessing) are described in detail in the paper [Simplifying Trip Data into Networks via Coupled Matrix Factorization](https://sites.google.com/site/chiatungkuo/publication/) and its supplementary materials.
The taxi data set used in the experiments is a public data set downloadable at [CRAWDAD Archive](http://crawdad.org/epfl/mobility/) after registration; More detail description of the data set can be found at their website.
