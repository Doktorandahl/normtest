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
	for(i in 1:nreps_outer){

		simulated_reg <- sim_tp_reg(ydist = ydist, y_params = y_params, n = n, xdists = xdists, x_params = x_params, nreps = nreps_inner)

		p_vals_dist[i,] <- three_dist_tests(simulated_reg$t,df = simulated_reg$df,alpha = alpha_dist)

		t1e_rate[i] <- mean(significance_checker(simulated_reg$t,df = simulated_reg$df,alpha = alpha_t1e))

		if(early_stop & i %% early_stop_interval == 0){
			if(min(apply(p_vals_dist[1:i,],2,mean)) < 0.06 | max(apply(p_vals_dist[1:i,],2,mean)) > 0.95){
				early_stopped <- TRUE
				break
			}
		}




	}
	colnames(p_vals_dist) <- paste("b1",c('ad','cvm','ks'),sep = '_')
	prop_reject <- p_vals_dist %>% apply(2,mean,na.rm = T)
	t1e_rate <- t1e_rate %>% mean(na.rm = T)
	out1 <- as_tibble_row(c(prop_reject,x1_t1e_rate = t1e_rate,n=n,i=i))
	out2 <- as_tibble_row(c(y = ydist, x1 = xdists[1], x2 = xdists[2], early_stop = early_stopped))
	return(bind_cols(out1,out2))
}



runner3 <- function(ydist,x1dist,x2dist,n,i,id_row){
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

	tmp <- sim_dist(ydist = ydist, y_params = y_params,
									x_params = list(x1_params,x2_params), xdist = c(x1dist,x2dist),
									n = n, nreps_inner = 1000, nreps_outer = 500, early_stop = T, early_stop_interval = 200, alpha_dist = 0.05, alpha_t1e = 0.05)
	cat('Finished sim ',i,' in', Sys.time()-t,". Early stop:", tmp$early_stop,'\n')
	saveRDS(bind_cols(tmp,id_row),paste0('res3/res3_',i,'.rds'))
}

grid_of_runs3 <- expand_grid(
	ydist = c('rnorm','runif','rlnorm0.5','rlnorm1','rlnorm2',
						'rlaplace','rbeta0.1','rbeta_lowskew','rbeta_medskew','rbeta_highskew'),
	x1dist = c('rnorm','runif','rlnorm0.5','rlnorm1','rlnorm2',
						 'rlaplace','rbeta0.1','rbeta_lowskew','rbeta_medskew','rbeta_highskew'),
	x2dist = c('rnorm','rlaplace','rlnorm1','rlnorm0.5','rlnorm2','rbeta0.1','runif'),

	n = c(seq(from = 6, to = 28, by = 2),
				seq(from = 30, to = 50, by = 5),
				seq(from = 60, to = 100, by = 10),
				seq(from = 120, to = 200, by = 20),
				seq(from = 250, to = 500, by = 50),
				seq(from = 600, to = 1000, by = 100),
				seq(from = 1250, to = 2000, by = 250),
				seq(from = 2500, to = 5000, by = 500),
				seq(from = 6000, to = 10000, by = 1000)))

grid_of_runs3 <- grid_of_runs3 %>% arrange(n)

already_run <- list.files('res3') %>% str_remove('res3_') %>% str_remove('.rds') %>% as.numeric()

if(length(already_run) == 0){
	already_run <- 0
}
remaining <- setdiff(1:nrow(grid_of_runs3),already_run) %>% sort()

registerDoParallel(cores = 6)

foreach(k = 1:nrow(grid_of_runs3)) %dopar%{
	if(k %in% remaining){
		cat('starting sim ',k,' of', nrow(grid_of_runs3),'n=',grid_of_runs3$n[k], '\n')
		runner3(grid_of_runs3$ydist[k], grid_of_runs3$x1dist[k], grid_of_runs3$x2dist[k], grid_of_runs3$n[k],k,grid_of_runs3[k,1:2])
	}
}





