# Workflow evaluations

This folder contains code and example data used to compare results of the same input data processed by original implementations and our Galaxy worflows. Strong correlations and minimal differences are expected, the contrary would indicate errors in the Galaxy implementation.

However, perfect match is not realistic -- the complex workflows simply cannot be reimplemented to every detail (including reference inputs). On the contrary, the differences we see here support our claim in the paper (TODO: cite) that the results are sensitive to any minor change, therefore keeping frozen reference implementations is essential to reproduce results, as well as to make comparisons among studies etc.

The folders contain both original and our outputs run on an experimental sample set. The raw input files are not provided, those are not public.

The evaluation is implemented in Jupyter notebooks which require common Python packages only (NumPy, Pandas, Matplotlib). They can be easily run in the pre-canned ``jupyter/scipy-notebook`` Docker container ([run.sh](run.sh)).

## JAX RNASeq pipeline

Jupyter notebook [jax-rnaseq.ipynb](jax-rnaseq/jax-rnaseq.ipynb)

Raw inputs were kindly provided by JAX and we processed them with our implementation. 

Snapshot of the JAX results is stored in the repository, or it can be downloaded fresh from their website (code is included in the notebook). These data are available as z-scores computed on larger sample set, we invert the z-score computation (using JAX's mean and standard deviation tables) to arrive back to log2 unnormalized expression counts we are able to compare. 


## UNITO RNASeq pipeline

Jupyter notebook [unito-rnaseq.ipynb](unito-rnaseq/unito-rnaseq.ipynb)

Both UNITO original and our results are provided on selected samples. Comparison is done on log2 expression counts again.

