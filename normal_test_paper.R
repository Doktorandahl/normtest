library(tidyverse)
library(foreach)
library(doParallel)
library(ranger)
library(extraDistr)

sim_tp_reg <- function(ydist,
											y_params = list(),
											n = 10,
											xdists = NULL,
											x_params = list(),
											nreps = 1000){

		if(length(x_params) == 0){
			x_params <- rep(list(list()),length(xdists))
		}

	if(str_detect(ydist,'rlnorm')){
		y_correction <- exp(mean(y_params[[1]]) + 0.5*mean(y_params[[2]])^2)
	}else if(str_detect(ydist, 'rbeta')){
		y_correction <- y_params[[1]]/(y_params[[1]] + y_params[[2]])
	}else if(ydist == 'runif'){
		y_correction <- 0.5
	}else{
		y_correction <- 0
	}

	t <- p <- matrix(NA,nreps,1+length(xdists))
	kurt <- skew <- rep(NA,nreps)
	for(i in 1:nreps){
		y <- do.call(ydist,c(n=n,y_params)) - y_correction
		if(!is.null(xdists)){
			x <- map2(xdists,x_params[1:length(xdists)],~do.call(.x,c(n=n,.y))) %>% unlist() %>% matrix(n,length(xdists))
			m <- lm(y~x)
		}else{
			m <- lm(y~1)
		}
		t[i,] <- summary(m)$coef[,3]
		kurt[i] <- moments::kurtosis(resid(m))
		skew[i] <- moments::skewness(resid(m))
	}
	return(list(t = t,kurt = kurt,skew = skew,df = m$df.residual))
}

three_dist_tests <- function(x,df, alpha){

	return(c(ad = goftest::ad.test(x,null = 'pt',df = df,estimated = F)$p.value < alpha,
					 cvm = goftest::cvm.test(x,null = 'pt',df = df,estimated = F)$p.value < alpha,
					 ks = ks.test(x,'pt',df = df)$p.value < alpha))
}

significance_checker <- function(x,df,alpha = 0.05,pars){x

	return(pt(abs(x),df = df,lower.tail = F)*2 < alpha)
}

sim_dist <- function(ydist,
										 y_params = list(),
										 n,
										 xdists = NULL,
										 x_params = list(),
										 nreps_inner,
										 nreps_outer,
										 early_stop = TRUE,
										 early_stop_interval = 200,
										 alpha_dist = 0.05,
										 alpha_t1e = 0.05,
										 correction = 'bonferroni'){

	p_vals_dist <- matrix(NA,nreps_outer,3*(length(xdists)+1))
	t1e_rate <- matrix(NA,nreps_outer,(length(xdists)+1))

	if(is.null(xdists)){
		colnames(p_vals_dist) <- paste('b0',c('ad','cvm','ks'),sep = '_')
		colnames(t1e_rate) <- paste('y','t1e_rate',sep = '_')
	}else{
		names(xdists) <- paste0('x',1:length(xdists))
		colnames(p_vals_dist) <- paste(rep(c('b0',paste0('b',1:length(xdists))),each=3),c('ad','cvm','ks'),sep = '_')
		colnames(t1e_rate) <- paste(c('y',paste0('x',1:length(xdists))),'t1e_rate',sep = '_')

	}

	skews <- kurtoses <- rep(NA,nreps_outer)

	for(i in 1:nreps_outer){
	simulated_reg <- sim_tp_reg(ydist = ydist, y_params = y_params, n = n, xdists = xdists, x_params = x_params, nreps = nreps_inner)
		skews[i] <- mean(simulated_reg$skew)
		kurtoses[i] <- mean(simulated_reg$kurt)
		p_vals_dist[i,] <- foreach(j = 1:(length(xdists)+1),.final=unlist) %do%
				three_dist_tests(simulated_reg$t[,j],df = simulated_reg$df,alpha = alpha_dist)

		t1e_rate[i,] <- foreach(j = 1:(length(xdists)+1),.final=unlist) %do%
				mean(significance_checker(simulated_reg$t[,j],df = simulated_reg$df,alpha = alpha_t1e))

		if(early_stop & i %% early_stop_interval == 0){
			err_rate_min <- p_vals_dist[1:i,] %>% apply(2,mean) %>% max()
			err_rate_max <- p_vals_dist[1:i,] %>% apply(2,mean) %>% min()
			if(err_rate_min < 0.01 | err_rate_max > 0.99){
				break
			}
		}




	}
	prop_reject <- p_vals_dist %>% apply(2,mean,na.rm = T)
	t1e_rate <- t1e_rate %>% apply(2,mean,na.rm = T)
	out1 <- as_tibble_row(c(prop_reject,t1e_rate,r_skew = mean(skews,na.rm=T), r_kurtosis = mean(kurtoses,na.rm=T),n=n,i=i))
	out2 <- as_tibble_row(c(y = ydist, xdists))
	return(bind_cols(out1,out2))
}



runner <- function(ydist,x1dist,x2dist,x3dist,n,i){
	t <- Sys.time()
	if(ydist == 'rlnorm0.5' | ydist == 'rlnorm1' | ydist == 'rlnorm1.5' | ydist == 'rlnorm2'){
		y_params <- list(meanlog = 0,sdlog = as.numeric(str_remove(ydist,'rlnorm')))
		ydist <- 'rlnorm'
	}else if(ydist == 'rbeta0.1'){
		y_params <- list(shape1 = 0.1, shape2 = 0.1)
		ydist <- 'rbeta'
	}else if(ydist == 'rbeta_lowskew'){
		y_params <- list(shape1 = 5, shape2 = 2)
		ydist <- 'rbeta'
	}else if(ydist == 'rbeta_medskew'){
		y_params <- list(shape1 = 5, shape2 = 1)
		ydist <- 'rbeta'
	}else if(ydist == 'rbeta_highskew'){
		y_params <- list(shape1 = 5, shape2 = 0.5)
		ydist <- 'rbeta'
	}else{
		y_params <- list()
	}

	if(x1dist == 'rlnorm0.5' | x1dist == 'rlnorm1' | x1dist == 'rlnorm1.5' | x1dist == 'rlnorm2'){
		x1_params <- list(meanlog = 0,sdlog = as.numeric(str_remove(x1dist,'rlnorm')))
		x1dist <- 'rlnorm'
	}else if(x1dist == 'rbeta0.1'){
		x1_params <- list(shape1 = 0.1, shape2 = 0.1)
		x1dist <- 'rbeta'
	}else if(x1dist == 'rbeta_lowskew'){
		x1_params <- list(shape1 = 5, shape2 = 2)
		x1dist <- 'rbeta'
	}else if(x1dist == 'rbeta_medskew'){
		x1_params <- list(shape1 = 5, shape2 = 1)
		x1dist <- 'rbeta'
	}else if(x1dist == 'rbeta_highskew'){
		x1_params <- list(shape1 = 5, shape2 = 0.5)
		x1dist <- 'rbeta'
	}else if(x1dist == 'NULL'){
		x1_params <- NULL
		x1dist <- NULL
	}else{
		x1_params <- list()
	}

	if(x2dist == 'rlnorm0.5' | x2dist == 'rlnorm1' | x2dist == 'rlnorm1.5' | x2dist == 'rlnorm2'){
		x2_params <- list(meanlog = 0,sdlog = as.numeric(str_remove(x2dist,'rlnorm')))
		x2dist <- 'rlnorm'
	}else if(x2dist == 'rbeta0.1'){
		x2_params <- list(shape1 = 0.1, shape2 = 0.1)
		x2dist <- 'rbeta'
	}else if(x2dist == 'rbeta_lowskew'){
		x2_params <- list(shape1 = 5, shape2 = 2)
		x2dist <- 'rbeta'
	}else if(x2dist == 'rbeta_medskew'){
		x_params <- list(shape1 = 5, shape2 = 1)
		x2dist <- 'rbeta'
	}else if(x2dist == 'rbeta_highskew'){
		x2_params <- list(shape1 = 5, shape2 = 0.5)
		x2dist <- 'rbeta'
	}else if(x2dist == 'NULL'){
		x2_params <- NULL
		x2dist <- NULL
	}else{
		x2_params <- list()
	}

	if(x3dist == 'rlnorm0.5' | x3dist == 'rlnorm1' | x3dist == 'rlnorm1.5' | x3dist == 'rlnorm2'){
		x3_params <- list(meanlog = 0,sdlog = as.numeric(str_remove(x3dist,'rlnorm')))
		x3dist <- 'rlnorm'
	}else if(x3dist == 'rbeta0.1'){
		x3_params <- list(shape1 = 0.1, shape2 = 0.1)
		x3dist <- 'rbeta'
	}else if(x3dist == 'rbeta_lowskew'){
		x3_params <- list(shape1 = 5, shape2 = 2)
		x3dist <- 'rbeta'
	}else if(x3dist == 'rbeta_medskew'){
		x3_params <- list(shape1 = 5, shape2 = 1)
		x3dist <- 'rbeta'
	}else if(x3dist == 'rbeta_highskew'){
		x3_params <- list(shape1 = 5, shape2 = 0.5)
		x3dist <- 'rbeta'
	}else if(x3dist == 'NULL'){
		x3_params <- NULL
		x3dist <- NULL
	}else{
		x3_params <- list()
	}

	tmp <- sim_dist(ydist = ydist, y_params = y_params,
					 x_params = list(x1_params,x2_params,x3_params), xdist = c(x1dist,x2dist,x3dist),
					 n = n,
					 nreps = 1000, nreps_inner = 1000, nreps_outer = 500, early_stop = T, early_stop_interval = 100, alpha_dist = 0.05, alpha_t1e = 0.05)
	cat('Finished sim ',i,' in', Sys.time()-t,'\n')
	saveRDS(tmp,paste0('normsim_res/sim_',i,'.rds'))
}

grid_of_runs <- expand_grid(
	ydist = c('rnorm','runif','rlnorm0.5','rlnorm1','rlnorm1.5','rlnorm2',
						'rlaplace','rbeta0.1','rbeta_lowskew','rbeta_medskew','rbeta_highskew'),
	x1dist = c('NULL','rnorm','runif','rlnorm0.5','rlnorm1','rlnorm1.5','rlnorm2',
						 'rlaplace','rbeta0.1','rbeta_lowskew','rbeta_medskew','rbeta_highskew'),
	x2dist = c('NULL','rnorm','rlnorm1'),
	x3dist = c('NULL','rnorm','rlnorm1'),
	n = c(6:30)
)

grid_of_runs <- grid_of_runs %>% filter(!(x2dist == 'NULL' & x3dist != 'NULL'),
																				!(x1dist == 'NULL' & x2dist != 'NULL'),
																				!(x1dist == 'NULL' & x3dist != 'NULL'),
																				!(x1dist != 'rnorm' & x2dist == 'rnorm'),
																				!(x1dist != 'rnorm' & x3dist == 'rnorm'),
																				!(x2dist != 'rnorm' & x3dist == 'rnorm'))

grid_of_runs <- grid_of_runs %>% arrange(x3dist,x2dist,x1dist,ydist,n)

already_run <- list.files('normsim_res') %>% str_remove('sim_') %>% str_remove('.rds') %>% as.numeric()

if(length(already_run) == 0){
	already_run <- 0
}
remaining <- setdiff(1:nrow(grid_of_runs),already_run) %>% sort()

registerDoParallel(cores = 10)

foreach(k = 1:nrow(grid_of_runs)) %dopar%{
	if(k %in% remaining){
		cat('starting sim ',k,' of', nrow(grid_of_runs), '\n')
		runner(grid_of_runs$ydist[k], grid_of_runs$x1dist[k], grid_of_runs$x2dist[k], grid_of_runs$x3dist[k], grid_of_runs$n[k],k)
	}
}






