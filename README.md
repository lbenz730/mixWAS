# mixWAS

This repository contains code for the paper (INSERT PAPER + Link Here). An R package for the `mixWAS` algorithm is available for download as follows

```
# install.packages('devtools')
devtools::install_github('lbenz730/mixWAS')
```

`mixWAS` exports the following 4 functions, with the following arguments


* `mixWAS`: mixWAS algorithm for all sites (if all data can be supplied at once). Returns either p-value for SNP or score/variance components.
    * `snps`: list of snps (one for each site), each a vector of SNPs $\in \{0,1,2\}$
    * `phenotypes`: list of phenotypes, each a matrix of phenotypes (one per column), with names of phenotypes specified as column names
    * `covariates`: list of covariates, each a matrix or data frame of covariates
    * `covariate_map`: Default = `NULL`. If `NULL`, function assumes by all covariates are to be used for each phenotype. If this is not desired behavior, user can supply a data frame with two columns one called `variable` and a second called `phenotype`. In the variable column is the name of covariates, with phenotypes being specified as 'all' (to use the variable for all phenotypes) or the name of a phenotype. Variables can be entered multiple times if they go to multiple phenotypes (but not all). A phenotype specific data set of covariates with use 'all' covariate + phenotype specific covariates.
    * `phenotype_index`: list of vectors giving the index (numeric) of which phenotypes are in each site's matrix. If `NULL` (default), will be inferred from matrix colnames.
    * `types`: optional vector specifying data types e.g. ('continuous', 'binary', 'count'). Default = `NULL` (phenotype data types will be inferred), Note that 'count' will never be inferred, only 'binary' or 'continuous'.
    * `parallel_sites`: logical, if score/variance component computations should be parallelized over sites. Default = `FALSE`.
    * `return_p`: logical, if `TRUE` return P-values, else return components like score/variance. Default = `TRUE`.

* `mixWAS_single_site`: Compute score vector and covariance matrix for a single site via mixWAS algorithm.
    * `snps`: matrix of phenotypes (one per column), with names of phenotypes specified as column names.
    * `phenotypes`: matrix or data frame of covariates.
    * `covariates`: list of covariates, each a matrix or data frame of covariates.
    * `covariate_map`: Default = `NULL`. If `NULL`, function assumes by all covariates are to be used for each phenotype. If this is not desired behavior, user can supply a data frame with two columns one called `variable` and a second called `phenotype`. In the variable column is the name of covariates, with phenotypes being specified as 'all' (to use the variable for all phenotypes) or the name of a phenotype. Variables can be entered multiple times if they go to multiple phenotypes (but not all). A phenotype specific data set of covariates with use 'all' covariate + phenotype specific covariates.
    * `types`: optional vector specifying data types ('continuous', 'binary', 'count'). Default = NULL (phenotype data types will be inferred). Note that 'count' will never be inferred, only 'binary' or 'continuous'.
    
* `run_hypothesis_test`: run hypothesis tests from intermediate mixWAS components
  * `score`: score vector from mixWAS intermediate output
  * `V_inv`: Inverse Varariance Matrix from mixWAS intermediate output
  * `z`: Standardized Z-scores from mixWAS intermediate output
  * `q`: # of phenotypes}

* `combine_site_results`: Combine results from running mixWAS on each site individually
    * `mixWAS_components`:	list of  `mixWAS_single_site` output of length = # of sites
    * `q`: # of phenotypes
    * `phenotypes`:	Optional list of phenotypes, each a matrix of phenotypes (one per column), with names of phenotypes specified as column names. If `phenotype_index` is `NULL`, will be used to infer phenotypes. One of `phenotypes` and `phenotype_index` must be specified
    * `phenotype_index`:	list of vectors giving the index (numeric) of which phenotypes are in each site's matrix. If `NULL` (default), will be inferred from `phenotypes`.
    * `return_p`: logical, if TRUE return P-values, else return components like score/variance. Default = `TRUE`.

