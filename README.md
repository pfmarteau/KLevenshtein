# kLevenshtein: a Regularized String Kernel derived from the Levenshtein's distance.

**kLevenshtein_cdist.py** is a python3.* implementation of kLevenshtein, a similarity measure dedicated to string matching. 

**kLevenshtein_c** directory contains a C biding of kLevenshtein, nearly 40 times faster that kLevenshtein_cdist implementation. To install the C version, just execute inside the K_Levenshtein subdirectory:
> python3 setup.py install --user

and to test:
> python3 test.py

KLevenshtein is derived from Levenshtein's distance while ensuring the property that kLevenshtein(.,.) is a positive definite kernel (homogeneous to an inner product in the so-called Reproducing Kernel Hilbert Space). The procedure to derive this kernelization of Levenshtein's distance is detailed in Marteau & Gibet 2014 [3]. kLevenshtein is a R-Convolution kernel as defined in [2]. 

kLevenshtein comes with two smothing meta-parameters entering into the computation of the local kernel: 

*k(x,y) = 1/3.(exp{-d(x,y)/sigma}+epsilon)/(1+epsilon)*, where d(x,y) is the hamming distance.
* *sigma* is used to scale the local distance computations (default value is set to 1.0).
* *epsilon* is used to avoid the vanishing of products of local kernel values (assimilated to probabilities) evaluated along the alignment paths (default value is set to 1e-3).

## Difference with Levenshtein distance
KLevenshtein is more resilient to sub-sequence flipping as shown in the following example, thanks to a high epsilon value (sigma=5, epsilon=1):

| *String1*	          |	*String2*	             | *Levenshtein distance*|	*Klevenshtein similarity*	|
|---------------------|------------------------|:---------------------:|:--------------------------:|
| the little big man	|	the big little man	   | **8**	               |	**0.0071**	              |
| the little big man	|	the little mairmaid	   | **6**	               |	**0.0065**	              |
| the little mairmaid	|	the big little man	   | 10	                   |	0.0053	                  |
| the little mairmaid	|	an old sailboat at sea |	20	                 |	0.0025	                  |


## Algorithmic complexity
The algorithmic complexity of kLevenshtein is O(N^2), where N is the length of the longest of the two time series given in entry, and when no corridor size is specified. This complexity drops down to O(C.N) when a symmetric corridor of size C is exploited. 

## Implementations
* kLevenshtein_cdist.py uses numpy and scipy.spatial.distance.cdist to accelerate the computation of the local kernel matrix
* The C/Python binding implementation in **kLevenshtein_c** directory is 40 times faster.

## Interpretation
The main idea leading to the construction of positive definite kernels from a given elastic distance used to match two time series is simply the following: instead of keeping only one of the best alignment paths, the new kernel sums up the costs of all the existing sub-sequence alignment paths with some smoothing factor that will favor good alignments while penalizing bad alignments. In addition, this smoothing factor can be optimized. Thus, basically we replace the min (or max) operator by a summation operator and introduce a symmetric corridor function (the function h in the recursive equation above) to optionally limit the summation and the computational cost. Then we add a regularizing recursive term (Kxx) such that the proof of the positive definiteness property can be understood as a direct consequence of the Haussler’s convolution theorem  [2]. Please see  [3] for the details.


## References

[2]   David Haussler. Convolution kernels on discrete structures. Technical Report UCS-CRL-99-10, University of California at Santa Cruz, Santa Cruz, CA, USA, 1999.

[3]   Pierre-François Marteau and Sylvie Gibet. On Recursive Edit Distance Kernels with Application to Time Series Classification. IEEE Transactions on Neural Networks and Learning Systems, pages 1–14, June 2014. 14 pages.

## see also 
[KDTW](https://github.com/pfmarteau/KDTW), a regularized Dynamic Time Warping (DTW) Kernel
