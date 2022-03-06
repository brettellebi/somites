########################
# GridLMM functions
########################

my.invnorm = function(x) {
  res = rank(x)
  res = qnorm(res/(length(res)+0.5))
  return(res)
}

run_gwas <- function(d,m,p,invers_norm=F, covariates = NULL) {
  # Creates a vector of the prefix "snp" combined with the row number of positions
  ids = paste("snp", 1:nrow(m), sep="")
  # Adds this vector as first column of the DF of positions: now comprising `ids`, `CHROM`, `BIN_START`, `BIN_END`
  mm = data.frame(ids, m)
  # Also add these ids as column names to the sample x locus DF of genotypes 
  colnames(d) = ids
  # Rename columns of positions DF
  colnames(mm) = c("snp", "Chr", "pos")
  # Then pull just the first three columns
  mm = mm[,1:3]
  # Rename it as `map`
  map = mm
  # Get number of samples (rows) from genotype DF
  n_geno=nrow(d)
  # Create matrix of 1200 elements from 1 to 1200, with 10 columns
  Plots = matrix(1:1200,nc = 10)
  # Create matrix of row numbers for each element of the `Plots` matrix (1:1200)
  Plot_Row = row(Plots)
  # Create matrix of col numbers for each element of the `Plots` matrix (1:10)
  Plot_Col = col(Plots)
  # Create a data frame with "Genox" for each sample in the first column
  # And randomly sample `n_geno` times from `Plots`
  data = data.frame(Geno = paste0('Geno',1: n_geno), Plot = sample(Plots)[1:nrow(d)])
  # Add a column `Row` with the row number of each randomly-sampled element
  data$Row = Plot_Row[data$Plot]
  # Add a column `Col` with the row number of each randomly-sampled element
  data$Col = Plot_Col[data$Plot]
  # Add phenotype
  if(invers_norm) {
    ## Inverse-normalise
    data = cbind(data,
                 y = my.invnorm(p %>% 
                              dplyr::select(y = dplyr::all_of(TARGET_PHENO))), # inverse-normalise phenotype
                 p %>% 
                   dplyr::select(dplyr::all_of(covariates))) # pull covariates
  } else {
    data = cbind(data,
                 p %>% 
                   dplyr::select(y = dplyr::all_of(TARGET_PHENO), # pull phenotype
                                 dplyr::all_of(covariates))) # and covariates
  }
  # Convert genotypes DF to matrix
  X = as.matrix(d)
  # Convert NAs to 0
  X[is.na(X)]=0
  # Add sample names to `X` as row names
  row.names(X) = data[,1]
  # Center genotypes by subtracting the mean for each marker (column) from the genotype value
  X_centered = sweep(X,2,colMeans(X),'-') # center marker genotypes
  # Use centered genotypes to create `K` kinship matrix
  K = tcrossprod(X_centered) / ncol(X_centered)
  # Make rownames and colnames of `K` the randomly sampled elements from `Plots`
  rownames(K) = colnames(K) = data$Plot
  # Calculate the spatial kernel
  field = data[,c('Row','Col')] # Pull out `Row` and `Col` columns
  dists = as.matrix(dist(field)) # Compute the distances 
  h = median(dists) # Calculate default tuning parameter `h`
  K_plot = KRLS::gausskernel(field,h^2/2); diag(K_plot)=1 # Compute the NxN distance matrix with pairwise with pairwise differences between rows as measured by a Gaussian Kernel
  rownames(K_plot) = colnames(K_plot) = data$Plot
  # Create `test_formula`
  if (is.null(covariates)){
    TEST_FORMULA = "~1"
  } else if (!is.null(covariates)){
    TEST_FORMULA = paste("~1 +", paste(covariates, collapse = " + "))
  }
  TEST_FORMULA = as.formula(TEST_FORMULA)
  # Run gwas
  gwas = GridLMM::GridLMM_GWAS(
    formula = y~1 + (1|Geno), # the same error model is used for each marker. It is specified similarly to lmer
    test_formula = TEST_FORMULA, # this is the model for each marker. ~1 means an intercept for the marker. ~1 + cov means an intercept plus a slope on `cov` for each marker
    reduced_formula = ~0, # This is the null model for each test. ~0 means just the error model. ~1 means the null includes the intercept of the marker, but not anything additional
    data = data, # The dataframe to look for terms from the 3 models
    weights = NULL, # optional observation-specific weights
    X = X, # The matrix of markers. Note: This can be of dimension n_g x p, where n_g is the number of unique genotypes.
    X_ID = 'Geno', # The column of the data that identifies each genotype. Each level of data$Geno should match a rowname of X
    h2_start = NULL, # An optional vector of h2s to use as starting values for each random effect in the model. If NULL, will be calculated from the error model using GridLMM_ML
    h2_step = 0.01, # step size per random effect for testing alternate values of h2
    max_steps = 100, # maximum number of steps of size h2_step to take from h2_start
    X_map = map, # Optional. The marker positions.
    relmat = list(Plot = K), # A list of Kernel matrices for the random effects. If X_ID (here Geno) is not included in this list, then it is calculated as tcrossprod(Xc)/ncol(Xc) where Xc is the centered (and optionally scaled) X. If any random effects are described in `error_model` but not provided here, the Kernel is assumed to be the identity matrix
    centerX = TRUE, # Should the markers be centered when calculating the GRM (only will be done if needed for calculating the GRM),
    scaleX = FALSE, # Should the markers be scaled to have constant variance when calculating the GRM?
    fillNAX = FALSE, # Should missing marker data be filled in with the mean allele frequency?
    method = 'REML', # REML = Wald test, ML = LRT, BF = calculate Bayes factors
    mc.cores = my_detectCores(), # How many cores should be used for parallel processing. Unless X is large, tends to actually be faster with mc.cores = 1
    verbose = FALSE # Should progress be printed to the screen?
  )
  return(gwas)
}

