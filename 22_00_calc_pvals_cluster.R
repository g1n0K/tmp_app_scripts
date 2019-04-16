##### Calc Pvalue and FDR
### calc permutation Pvalue #86
### IV id perm dist, DV id permutation P-value
# 150712

args <- commandArgs()
i0=as.numeric(args[9])
i1=as.numeric(args[10])
fileID <- paste(i0,i1,sep='_')

log <- paste(fileID,' start ',date(),sep='')
write(log,file='log.txt',append=T)

### file path
Eucfile='/scratch/gino/fdr/input/Euc_dist.Rdata'
perm='/scratch/gino/fdr/input/Perm_dist.Rdata'

setwd('/scratch/gino/fdr/log/')

### def functions
# f_fdr: calc False discovery rate
# f_Natural: natural number
# bumfix: beta uniform mixture model from the package, 'dnet'. used 'Broyden–Fletcher–Goldfarb–Shanno'(quasi-Newton methods)

f_fdr <- function(x)
{
	fp <- (lambda + (1 - lambda) * a) * x
	tp <- (1 - lambda) * (x^a - a * x)
	fp/(fp + tp)
}

f_Natural <- function (n) 
{
    stopifnot(is.numeric(n))
    floor(n) == ceiling(n) & n >= 1 & n <= 2^53 - 1
}

# try five times for avoiding fitting error
bumfix <- function(x) 
{
    ntry = 5
    a <- runif(ntry, 0.1, 0.9)
    lambda <- runif(ntry, 0.1, 0.9)
    fbum <- function(x, lambda, a) { lambda + (1 - lambda) * a * x^(a - 1) }
    fn <- function(parms, x) { -1 * sum(log(fbum(x, parms[1], parms[2]))) }
    gr <- function(parms, x) {
        lambda <- parms[1]
        a <- parms[2]
        deno <- lambda + (1 - lambda) * a * x^(a - 1)
        d_lambda <- -1 * sum((1 - a * x^(a - 1))/deno)
        d_a <- -1 * sum(((1 - lambda) * x^(a - 1) + a * (1 - 
            lambda) * x^(a - 1) * log(x))/deno)
        return(c(d_lambda, d_a))
    }
    value <- Inf
    best <- list()
    for (i in 1:ntry) {
        test.optim <- try(opt <- optim(c(lambda[i], a[i]), fn = fn, 
            gr = gr, x = x, lower = rep(1e-05, 3), method = "L-BFGS-B", 
            upper = rep(1 - 1e-05, 3)))
        if ((!class(test.optim) == "try-error") && all(opt$par >= 
            1e-05) && all(opt$par <= 1 - 1e-05)) {
            if (opt$value < value) {
                value <- opt$value
                best <- opt
            }
            if (0) {
                message(sprintf("Try: %d\t with Log-likelihood: %.1f, Mixture parameter (lambda): %1.3f, Shape parameter (a): %1.3f", 
                  i, -1 * opt$value, opt$par[1], opt$par[2]))
            }
        }
    }
    ## fix here
    if (length(best) == 0) {
    	F <- as.list(c(-9,-9))
    	names(F) <- c('a','lambda')
        return(F)
        stop('fail')
    }
    else {
        if (any(opt$par == 1e-05) || any(opt$par == 1 - 1e-05)) {
        }
        fit <- list(lambda = best$par[1], a = best$par[2], NLL = best$value, 
            pvalues = x, call = match.call())
        return(fit)
    }
}

# make result set
# from A to B # from column to row
# ig. Gene A's distance to F # column A, row F
#permp <- perm P-value set #par <- BUM parameter set #fdr <- false rate

# import distance matrix
load(file=Eucfile)
load(file=perm)

# set up result table
headers <- rownames(Euc)
vec <- length(headers)

perm_p <- matrix(NA,nrow=nrow(Euc),ncol=nrow(Euc))
rownames(perm_p) <- headers
colnames(perm_p) <- headers

par <- as.data.frame(headers)
colnames(par) <- c('symbol_column')
par$alpha <- NA
par$lambda <- NA
par$tau <- NA

perm_size<-100000
fdr <- perm_p

log <- paste(fileID,' began ',date(),sep='')
write(log,file='log.txt',append=T)

perm <- as.data.frame(perm_dist)
colnames(perm) <- append('dist',headers)

for(j in i0:i1)
{
	if( (j%%10) == 0) {
		log <- paste(fileID,j,'out of',i0,date(),sep=' ')
		write(log,file='log.txt',append=T)
	}
	
	# m is for euc, k is for perm
	m <- j
	k <- m+1
	pvals <- unlist(lapply(1:nrow(Euc), function(i)	sum( perm[1:min(which(Euc[i,m] <= perm$dist)),k] ) / perm_size ) )
	perm_p[,m] <- pvals
	# Calc FDR by Beta uniform mixture model
	# need replace 0 & 1 to inbound cuz cannor fit for FDR
	pvals[pvals == 0] <- 0.000001
	pvals[pvals == 1] <- 0.999999
	fit <- bumfix(pvals)
	a <- fit$a
	lambda <- fit$lambda
	tau <- lambda+((1-lambda)*a)
	
	## save parameters
	par$alpha[m] <- a
	par$lambda[m] <- lambda
	par$tau[m] <- tau
	
	if(a == -9)
	{
			log <- paste('fitting FAIL',fileID,par$symbol_column[m],sep=' ')
			write(log,file='rec_parameter.txt',append=T)
		}
	
	if(a == 1)
	{
			log <- paste('alpha=1',fileID,par$symbol_column[m],sep=' ')
			write(log,file='rec_parameter.txt',append=T)
		}
	if(lambda == 1)
	{
			log <- paste('lambda=1',fileID,par$symbol_column[m],sep=' ')
			write(log,file='rec_parameter.txt',append=T)
		}
	if(lambda == 0)
	{
			log <- paste('lambda=0',fileID,par$symbol_column[m],sep=' ')
			write(log,file='rec_parameter.txt',append=T)
		}
	
	for(i in 1:nrow(fdr))
	{
		tmp <- perm_p[i,m]
		if(tmp == 0) {fdr[i,m] <- f_fdr(0.000001)} else if(tmp == 1) {
			fdr[i,m] <- f_fdr(0.999999)} else {
			fdr[i,m] <- f_fdr(tmp)
			}
	}
	
}

log <- paste(fileID,'P-value FDR end: ',date(),sep=' ')
write(log,file='log.txt',append=T)

rownames(perm_p) <- headers
rownames(fdr) <- headers
colnames(perm_p) <- headers
colnames(fdr) <- headers

setwd('/scratch/gino/fdr/out/')
Pvals <- as.data.frame(perm_p)
FDR <- as.data.frame(fdr)

save(Pvals,file=paste(fileID,'Pvals.Rdata',sep='_'))
save(FDR,file=paste(fileID,'FDR.Rdata',sep='_'))

#Pvals_max <- Pvals
#Pvals_min <- Pvals

#k <- nrow(Pvals)
#for(i in 1:(nrow(Pvals)-1))
#{
#	j <- i+1
#	tmp1 <- as.numeric(Pvals[i,j:k])
#	tmp2 <- as.numeric(Pvals[j:k,i])
#	maxs <- unlist(lapply(1:length(tmp1), function(x) max(tmp1[x],tmp2[x]) ))
#	mins <- unlist(lapply(1:length(tmp1), function(x) min(tmp1[x],tmp2[x]) ))
#	Pvals_max[i,j:k] <- maxs
#	Pvals_max[j:k,i] <- maxs
#	Pvals_min[i,j:k] <- mins
#	Pvals_min[j:k,i] <- mins
#}

#save(Pvals_max,file='Pvals_max.Rdata')
#save(Pvals_min,file='Pvals_min.Rdata')

