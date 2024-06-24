library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(foreach)
library(doParallel)
library(ranger)
library(extraDistr)
source('normal_test_functions.R')


grid_of_runs3 <- expand_grid(
	ydist = c('rnorm','runif','rlnorm0.5','rlnorm1','rlnorm2',
						'rlaplace','rbeta0.1','rbeta_lowskew','rbeta_medskew','rbeta_highskew'),
	x1dist = c('rnorm','runif','rlnorm0.5','rlnorm1','rlnorm2',
						 'rlaplace','rbeta0.1','rbeta_lowskew','rbeta_medskew','rbeta_highskew'),
	x2dist = c('rnorm','rlaplace','rlnorm1','rlnorm0.5','rlnorm2','rbeta0.1','runif'),

	n = c(seq(from = 4, to = 28, by = 2),
				seq(from = 30, to = 50, by = 5),
				seq(from = 60, to = 100, by = 10),
				seq(from = 120, to = 200, by = 20),
				seq(from = 250, to = 500, by = 50),
				seq(from = 600, to = 1000, by = 100),
				seq(from = 1250, to = 2000, by = 250),
				seq(from = 2500, to = 5000, by = 500),
				seq(from = 6000, to = 10000, by = 1000)))

grid_of_runs3 <- grid_of_runs3 %>% arrange(n)


# if(length(already_run) == 0){
# 	already_run <- 0
# }
# remaining <- setdiff(1:nrow(grid_of_runs3),already_run) %>% sort()

registerDoParallel(cores = 256)

foreach(k = nrow(grid_of_runs3):1) %dopar%{
	cat('starting sim ',k,' of', nrow(grid_of_runs3),'n=',grid_of_runs3$n[k], '\n')
	runner3(grid_of_runs3$ydist[k], grid_of_runs3$x1dist[k], grid_of_runs3$x2dist[k], grid_of_runs3$n[k],k,grid_of_runs3[k,1:2])
}





