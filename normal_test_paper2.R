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
			x <- do.call(xdists,c(n=n,x_params))
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



runner2 <- function(ydist,x1dist,x2dist,x3dist,n,i, id_row){
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


	tmp <- sim_dist(ydist = ydist, y_params = y_params,
									x_params = x1_params, xdist = x1dist,
									n = n, nreps_inner = 1000, nreps_outer = 500, early_stop = T, early_stop_interval = 100, alpha_dist = 0.05, alpha_t1e = 0.05)
	cat('Finished sim ',i,' in', Sys.time()-t,'\n')
	saveRDS(bind_cols(tmp,id_row),paste0('res2/res2_',i,'.rds'))
}

grid_of_runs2 <- expand_grid(
	ydist = c('rnorm','runif','rlnorm0.5','rlnorm1','rlnorm1.5','rlnorm2',
						'rlaplace','rbeta0.1','rbeta_lowskew','rbeta_medskew','rbeta_highskew'),
	x1dist = c('NULL','rnorm','runif','rlnorm0.5','rlnorm1','rlnorm1.5','rlnorm2',
						 'rlaplace','rbeta0.1','rbeta_lowskew','rbeta_medskew','rbeta_highskew'),
	n = c(seq(from = 4, to = 28, by = 2),
				seq(from = 30, to = 50, by = 5),
				seq(from = 60, to = 100, by = 10),
				seq(from = 120, to = 200, by = 20),
				seq(from = 250, to = 500, by = 50),
				seq(from = 600, to = 1000, by = 100),
				seq(from = 1250, to = 2000, by = 250),
				seq(from = 2500, to = 5000, by = 500),
				seq(from = 6000, to = 10000, by = 1000)))


grid_of_runs2 <- grid_of_runs2 %>% arrange(n)

already_run <- list.files('res2') %>% str_remove('res2_') %>% str_remove('.rds') %>% as.numeric()

if(length(already_run) == 0){
	already_run <- 0
}
remaining <- setdiff(1:nrow(grid_of_runs2),already_run) %>% sort()

registerDoParallel(cores = 10)

foreach(k = 1:nrow(grid_of_runs2)) %dopar%{
	if(k %in% remaining){
		cat('starting sim ',k,' of', nrow(grid_of_runs2),'n=',grid_of_runs2$n[k], '\n')
		runner2(grid_of_runs2$ydist[k], grid_of_runs2$x1dist[k], grid_of_runs2$x2dist[k], grid_of_runs2$x3dist[k], grid_of_runs2$n[k],k,grid_of_runs2[k,1:2])
	}
}





