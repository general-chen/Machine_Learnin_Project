## This code is used for cross-validation
  - reads data from .csv and converts to pandas dataframe.
  - retrieves CP and CL, CD.
  - use Bessel filter to smooth pressure data, no filtering for CL and CD.
## This code uses 8-fold cross-validation for hypermarameters optimization of neurons and layers.
  - neurons: [16, 24, 32, 40, 48, 56, 64, 72]
  - layers:  [1,2,3,4,5,6,7,8]
## The best learning rate for this delta wing datasets is:
  - $10^{-5}$
