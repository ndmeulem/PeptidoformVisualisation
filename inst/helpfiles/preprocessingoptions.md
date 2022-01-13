Preprocessing options
======================

Logarithmic transformation
----------------------------
Oftentimes, intensity data show a very skewed distribution. Logarithmic transformation can help reduce this skewness, resulting in a distribution a bit closer to the Gaussian one. This can allow us to then fit a linear model to the data.
There are four options to choose from:
- **none**: no logarithmic transformation will be performed. For example, if the data have already been log transformed previously.
- **log2, log10 or natural**: which base with respect to which logarithms are computed. natural: e = exp(1)

Minimum number of non zero columns
------------------------------------
How many intensity values there should at least be present for each peptidoform, i.e. in at least how many samples a quantitative value has been picked up.
Defaults to at least two samples. Any lower than this value might lead to too many missing values in the computations.

Normalisation
---------------
Normalisation to be performed on the entire dataset.  
There are three options to choose from:
- **none**: no normalisation will be performed.
- **center.mean**: mean centering will be performed, i.e. the mean of a variable will be subtracted from all observations on that variable such that the new mean of the variable is zero.
- **center.median**: median centering will be performed, i.e. the median of a variable will be subtracted from all observations on that variable such that the new median of the variable is zero.

When all options are filled out, the user can click the preprocess button and all preprocessing functions will be carries out, as well as the calculation of a density plot and boxplot of the data after preprocessing. There is also a boxplot and density plot present of the original data as it was uploaded, for comparison.
Do note that the preprocessing should be carried out before the user can go on to the next tab.