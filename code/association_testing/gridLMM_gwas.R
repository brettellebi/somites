
my.invnorm = function(x) {
    res = rank(x)
    res = qnorm(res/(length(res)+0.5))
return(res)
}

gwas_object <- function(geno_file="", pheno_file="", position_file="") {
	d = read.table(geno_file)
	m = read.table(position_file)
	p = read.table(pheno_file, comment.char="n")
return(list("d"=d, "m"=m, "p"=p))
}


run_gwas <- function(d,m,p,invers_norm=F) {
	ids = paste("snp", 1:nrow(m), sep="")
  mm = data.frame(ids, m)
  colnames(d) = ids
  colnames(mm) = c("snp", "Chr", "pos")
  mm = mm[,1:3]
  map = mm
  n_geno=nrow(d)
  Plots = matrix(1:1200,nc = 10)
  Plot_Row = row(Plots)
  Plot_Col = col(Plots)
  data = data.frame(Geno = paste0('Geno',1: n_geno), Plot = sample(Plots)[1:nrow(d)])
  data$Row = Plot_Row[data$Plot]
	data$Col = Plot_Col[data$Plot]
	if(invers_norm) {
		data$y = my.invnorm(p[,2])
	} else {
		data$y = p[,2]
	}
	X = as.matrix(d)
	X[is.na(X)]=0
	row.names(X) = data[,1]
	X_centered = sweep(X,2,colMeans(X),'-') # center marker genotypes
	K = tcrossprod(X_centered) / ncol(X_centered)
	rownames(K) = colnames(K) = data$Plot
	field = data[,c('Row','Col')]
	dists = as.matrix(dist(field))
	h = median(dists)
	K_plot = gausskernel(field,h^2/2); diag(K_plot)=1 # 
	rownames(K_plot) = colnames(K_plot) = data$Plot
	gwas = GridLMM_GWAS(
	                        formula = y~1 + (1|Geno) + (1|Plot), # the same error model is used for each marker. It is specified similarly to lmer
	                        test_formula = ~1, # this is the model for each marker. ~1 means an intercept for the marker. ~1 + cov means an intercept plus a slope on `cov` for each marker
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

plot_gwas <- function(gwas, title="", sig_cut=0, sugest=0) {
	manhattan(gwas$results[complete.cases(gwas$results),],'Chr','pos','p_value_REML','snp', ylim=c(0,8), col=c("light blue", "dark blue"), suggestiveline = F, genomewideline = F, cex = 0.6, main=title)
	abline(h=sig_cut, col="red", lty="dashed")
		abline(h=sugest, col="blue", lty="dashed")
}

permute <- function(d,m,p) {
	p = p[order(runif(nrow(p), 1, nrow(p))),]
	result = run_gwas(d,m,p)
return(max(-log10(result$results$p_value_REML[complete.cases(result$results)])))
}

library(GridLMM)
library(KRLS)
library(qqman)

g_all = gwas_object(geno_file="all_genotypes.txt", pheno_file="phenotypes.txt", position_file="all_positions.txt")
gwas_all = run_gwas(g_all$d, g_all$m, g_all$p)
plot_gwas(gwas_all, "HdrR x Ho5", 3.6, 2.9)



