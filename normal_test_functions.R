library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(stringr)
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

	t <- p <- rep(NA,nreps)
	kurt <- skew <- rep(NA,nreps)
	for(i in 1:nreps){
		y <- do.call(ydist,c(n=n,y_params)) - y_correction
		x <- map2(xdists,x_params[1:length(xdists)],~do.call(.x,c(n=n,.y))) %>% unlist() %>% matrix(n,length(xdists))
		m <- lm(y~x)

		t[i] <- summary(m)$coef[2,3]
	}
	return(list(t = t,df = m$df.residual))
}

three_dist_tests <- function(x,df, alpha){

	return(c(ad = goftest::ad.test(x,null = 'pt',df = df,estimated = F)$p.value < alpha,
					 cvm = goftest::cvm.test(x,null = 'pt',df = df,estimated = F)$p.value < alpha,
					 ks = ks.test(x,'pt',df = df)$p.value < alpha))
}

significance_checker <- function(x,df,alpha = 0.05,pars){

	return(pt(abs(x),df = df,lower.tail = F)*2 < alpha)
}

sim_dist <- function(ydist,
										 y_params = list(),
										 n,
										 i,
										 xdists = NULL,
										 x_params = list(),
										 nreps_inner,
										 nreps_outer,
										 early_stop = TRUE,
										 early_stop_interval = 200,
										 alpha_dist = 0.05,
										 alpha_t1e = 0.05,
										 correction = 'bonferroni'){

	p_vals_dist <- matrix(NA,nreps_outer,3)
	t1e_rate <- rep(NA,nreps_outer)

	early_stopped <- FALSE
	for(j in 1:nreps_outer){

		simulated_reg <- sim_tp_reg(ydist = ydist, y_params = y_params, n = n, xdists = xdists, x_params = x_params, nreps = nreps_inner)

		p_vals_dist[j,] <- three_dist_tests(simulated_reg$t,df = simulated_reg$df,alpha = alpha_dist)

		t1e_rate[j] <- mean(significance_checker(simulated_reg$t,df = simulated_reg$df,alpha = alpha_t1e))


	}
	colnames(p_vals_dist) <- paste("b1",c('ad','cvm','ks'),sep = '_')
	prop_reject <- p_vals_dist %>% apply(2,mean,na.rm = T)
	t1e_rate <- t1e_rate %>% mean(na.rm = T)
	out1 <- as_tibble_row(c(prop_reject,x1_t1e_rate = t1e_rate,n=n,i=i))
	out2 <- as_tibble_row(c(y = ydist, x1 = xdists[1], x2 = xdists[2], early_stop = early_stopped))
	return(bind_cols(out1,out2))
}



runner3 <- function(ydist,x1dist,x2dist,n,i,id_row){
	filename <- paste0('res4/res4_',i %% 256,'.rds')
	if(file.exists(filename)){
		runs <- readRDS(filename)
		cat('File exists. Reading in\n')
		if(i %in% runs$i){
			cat('Already run. Skipping',i,'\n')
			return(NULL)
		}
	}else{
		runs <- NULL
		cat('File does not exist. Creating\n')
	}

	t <- Sys.time()
	if(ydist == 'rlnorm0.5' | ydist == 'rlnorm1' | ydist == 'rlnorm1.5' | ydist == 'rlnorm2'){
		y_params <- list(meanlog = 0,sdlog = as.numeric(str_remove(ydist,'rlnorm')))
		ydist <- 'rlnorm'
	}else if(ydist == 'rbeta0.1'){
		y_params <- list(shape1 = 0.2, shape2 = 0.2)
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
		x1_params <- list(shape1 = 0.2, shape2 = 0.2)
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
		x2_params <- list(shape1 = 0.2, shape2 = 0.2)
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

	cat('Starting sim ',i,' at', Sys.time()-t,"\n")

	tmp <- sim_dist(ydist = ydist, y_params = y_params,i=i,
									x_params = list(x1_params,x2_params), xdist = c(x1dist,x2dist),
									n = n, nreps_inner = 1000, nreps_outer = 500, early_stop = F, early_stop_interval = 200, alpha_dist = 0.05, alpha_t1e = 0.05)
	cat('Finished sim ',i,' in', Sys.time()-t,"\n")
	if(file.exists(filename)){
		runs <- readRDS(filename)
	}
	saveRDS(bind_rows(runs,bind_cols(tmp,id_row)),filename)
}
