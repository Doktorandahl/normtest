library(tidyverse)

res <- list.files('res2') %>%
	map(~readRDS(paste0('res2/',.x))) %>%
	bind_rows()

res <- res %>% mutate(symmetric_x1 = case_when(x1dist %in% c('rnorm','runif','rlaplace','rbeta0.1') ~ T,
																							 T~F),
											symmetric_y = case_when(ydist %in% c('rnorm','runif','rlaplace','rbeta0.1') ~ T,
																							T~F))




## Figure 1
res %>%
	filter(x1dist %in% c('runif','rlaplace','rbeta0.1','rbeta_lowskew','rnorm','rlnorm0.5','rlnorm1','rlnorm2'),
				 ydist %in% c('runif','rlaplace','rbeta0.1','rbeta_lowskew','rnorm','rlnorm1','rlnorm2','rlnorm0.5')) %>%
	mutate(ydist = factor(ydist,levels = c('rnorm','runif','rlaplace','rbeta0.1','rbeta_lowskew','rlnorm0.5','rlnorm1','rlnorm2'),
												labels = c('Normal','Uniform','Laplace','Beta(0.1. 0.1)','Beta (5,2)','Lognormal (0, 0.5)','Lognormal (0, 1)','Lognormal (0, 2)')),
				 x1dist = factor(x1dist,levels = c('rnorm','runif','rlaplace','rbeta0.1','rbeta_lowskew','rlnorm0.5','rlnorm1','rlnorm2'),
				 								labels = c('Normal','Uniform','Laplace','Beta(0.1. 0.1)','Beta (5,2)','Lognormal (0, 0.5)','Lognormal (0, 1)','Lognormal (0, 2)'))) %>%
	pivot_longer(b1_ad) %>%
	ggplot(aes(x=n,y=value,color=x1dist))+
	geom_line(lwd = 1)+
	geom_vline(xintercept = 30,linetype = 2,lwd=0.25)+
	scale_x_continuous(trans='log', breaks = c(4,8,15,30,60,120,250,500,1000,2000,4000,10000))+
	scale_color_manual(name = 'X Distribution', values = c(wesanderson::wes_palette('Darjeeling1')[1:6],wesanderson::wes_palette('Darjeeling2')[1:2]))+
	facet_wrap(~ydist,nrow=2)+
	theme_minimal()+
	xlab('Sample Size')+
	ylab('Proportion of rejected Anderson-Darling test')+
	theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
				axis.title = element_text(size = 14, face = 'bold'),
				strip.text = element_text(size = 11, face = 'bold'),
				legend.title = element_text(size = 11, face = 'bold'),
				legend.text = element_text(size = 11))

## Figure 1 full
res %>% filter(x1dist != 'NULL') %>%
	mutate(ydist = factor(ydist,levels = c('rnorm','runif','rlaplace','rbeta0.1','rbeta_lowskew','rbeta_medskew','rbeta_highskew','rlnorm0.5','rlnorm1','rlnorm1.5','rlnorm2'),
												labels = c('Normal','Uniform','Laplace','Beta(0.1. 0.1)','Beta (5,2)','Beta (5,1)','Beta (5,0.5)','Lognormal (0, 0.5)','Lognormal (0, 1)','Lognormal (0, 1.5)','Lognormal (0, 2)')),
				 x1dist = factor(x1dist,levels = c('rnorm','runif','rlaplace','rbeta0.1','rbeta_lowskew','rbeta_medskew','rbeta_highskew','rlnorm0.5','rlnorm1','rlnorm1.5','rlnorm2'),
				 								labels = c('Normal','Uniform','Laplace','Beta(0.1. 0.1)','Beta (5,2)','Beta (5,1)','Beta (5,0.5)','Lognormal (0, 0.5)','Lognormal (0, 1)','Lognormal (0, 1.5)','Lognormal (0, 2)'))) %>%
	pivot_longer(b1_ad) %>%
	ggplot(aes(x=n,y=value,color=x1dist))+
	geom_line(lwd = 1)+
	geom_vline(xintercept = 30,linetype = 2,lwd=0.25)+
	scale_x_continuous(trans='log', breaks = c(4,8,15,30,60,120,250,500,1000,2000,4000,10000))+
	scale_color_manual(name = 'X Distribution', values = c(wesanderson::wes_palette('Darjeeling1')[1:6],wesanderson::wes_palette('Darjeeling2')[1:5]))+
	facet_wrap(~ydist,nrow=2)+
	theme_minimal()+
	xlab('Sample Size')+
	ylab('Proportion of rejected Anderson-Darling test')+
	theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
				axis.title = element_text(size = 14, face = 'bold'),
				strip.text = element_text(size = 11, face = 'bold'),
				legend.title = element_text(size = 11, face = 'bold'),
				legend.text = element_text(size = 11))



## Figure 2
res %>%
	filter(x1dist %in% c('runif','rlaplace','rbeta0.1','rbeta_lowskew','rnorm','rlnorm0.5','rlnorm1','rlnorm2'),
				 ydist %in% c('runif','rlaplace','rbeta0.1','rbeta_lowskew','rnorm','rlnorm1','rlnorm2','rlnorm0.5')) %>%
	mutate(ydist = factor(ydist,levels = c('rnorm','runif','rlaplace','rbeta0.1','rbeta_lowskew','rlnorm0.5','rlnorm1','rlnorm2'),
												labels = c('Normal','Uniform','Laplace','Beta(0.1. 0.1)','Beta (5,2)','Lognormal (0, 0.5)','Lognormal (0, 1)','Lognormal (0, 2)')),
				 x1dist = factor(x1dist,levels = c('rnorm','runif','rlaplace','rbeta0.1','rbeta_lowskew','rlnorm0.5','rlnorm1','rlnorm2'),
				 								labels = c('Normal','Uniform','Laplace','Beta(0.1. 0.1)','Beta (5,2)','Lognormal (0, 0.5)','Lognormal (0, 1)','Lognormal (0, 2)'))) %>%
	pivot_longer(x1_t1e_rate) %>%
	ggplot(aes(x=n,y=value,color=x1dist))+
	geom_line(lwd = 1)+
	geom_vline(xintercept = 30,linetype = 2,lwd=0.25)+
	scale_x_continuous(trans='log', breaks = c(4,8,15,30,60,120,250,500,1000,2000,4000,10000))+
	scale_color_manual(name = 'X Distribution', values = c(wesanderson::wes_palette('Darjeeling1')[1:6],wesanderson::wes_palette('Darjeeling2')[1:2]))+
	facet_wrap(~ydist,nrow=2)+
	theme_minimal()+
	xlab('Sample Size')+
	ylab('Type I Error Rate')+
	theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
				axis.title = element_text(size = 14, face = 'bold'),
				strip.text = element_text(size = 11, face = 'bold'),
				legend.title = element_text(size = 11, face = 'bold'),
				legend.text = element_text(size = 11))



## Figure 2 full
res %>% filter(x1dist != 'NULL') %>%
	mutate(ydist = factor(ydist,levels = c('rnorm','runif','rlaplace','rbeta0.1','rbeta_lowskew','rbeta_medskew','rbeta_highskew','rlnorm0.5','rlnorm1','rlnorm1.5','rlnorm2'),
												labels = c('Normal','Uniform','Laplace','Beta(0.1. 0.1)','Beta (5,2)','Beta (5,1)','Beta (5,0.5)','Lognormal (0, 0.5)','Lognormal (0, 1)','Lognormal (0, 1.5)','Lognormal (0, 2)')),
				 x1dist = factor(x1dist,levels = c('rnorm','runif','rlaplace','rbeta0.1','rbeta_lowskew','rbeta_medskew','rbeta_highskew','rlnorm0.5','rlnorm1','rlnorm1.5','rlnorm2'),
				 								labels = c('Normal','Uniform','Laplace','Beta(0.1. 0.1)','Beta (5,2)','Beta (5,1)','Beta (5,0.5)','Lognormal (0, 0.5)','Lognormal (0, 1)','Lognormal (0, 1.5)','Lognormal (0, 2)'))) %>%
	pivot_longer(x1_t1e_rate) %>%
	ggplot(aes(x=n,y=value,color=x1dist))+
	geom_line(lwd = 1)+
	geom_vline(xintercept = 30,linetype = 2,lwd=0.25)+
	scale_x_continuous(trans='log', breaks = c(4,8,15,30,60,120,250,500,1000,2000,4000,10000))+
	scale_color_manual(name = 'X Distribution', values = c(wesanderson::wes_palette('Darjeeling1')[1:6],wesanderson::wes_palette('Darjeeling2')[1:5]))+
	facet_wrap(~ydist,nrow=2)+
	theme_minimal()+
	xlab('Sample Size')+
	ylab('Type I Error Rate')+
	theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
				axis.title = element_text(size = 14, face = 'bold'),
				strip.text = element_text(size = 11, face = 'bold'),
				legend.title = element_text(size = 11, face = 'bold'),
				legend.text = element_text(size = 11))





## Figure 3 full
res %>% filter(x1dist != 'NULL') %>%
	mutate(ydist = factor(ydist,levels = c('rnorm','runif','rlaplace','rbeta0.1','rbeta_lowskew','rbeta_medskew','rbeta_highskew','rlnorm0.5','rlnorm1','rlnorm1.5','rlnorm2'),
												labels = c('Normal','Uniform','Laplace','Beta(0.1. 0.1)','Beta (5,2)','Beta (5,1)','Beta (5,0.5)','Lognormal (0, 0.5)','Lognormal (0, 1)','Lognormal (0, 1.5)','Lognormal (0, 2)')),
				 x1dist = factor(x1dist,levels = c('rnorm','runif','rlaplace','rbeta0.1','rbeta_lowskew','rbeta_medskew','rbeta_highskew','rlnorm0.5','rlnorm1','rlnorm1.5','rlnorm2'),
				 								labels = c('Normal','Uniform','Laplace','Beta(0.1. 0.1)','Beta (5,2)','Beta (5,1)','Beta (5,0.5)','Lognormal (0, 0.5)','Lognormal (0, 1)','Lognormal (0, 1.5)','Lognormal (0, 2)'))) %>%
	pivot_longer(b1_cvm) %>%
	ggplot(aes(x=n,y=value,color=x1dist))+
	geom_line(lwd = 1)+
	geom_vline(xintercept = 30,linetype = 2,lwd=0.25)+
	scale_x_continuous(trans='log', breaks = c(4,8,15,30,60,120,250,500,1000,2000,4000,10000))+
	scale_color_manual(name = 'X Distribution', values = c(wesanderson::wes_palette('Darjeeling1')[1:6],wesanderson::wes_palette('Darjeeling2')[1:5]))+
	facet_wrap(~ydist,nrow=2)+
	theme_minimal()+
	xlab('Sample Size')+
	ylab('Proportion of rejected Cramer von Mies test')+
	theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
				axis.title = element_text(size = 14, face = 'bold'),
				strip.text = element_text(size = 11, face = 'bold'),
				legend.title = element_text(size = 11, face = 'bold'),
				legend.text = element_text(size = 11))

## Figure 4 full
res %>% filter(x1dist != 'NULL') %>%
	mutate(ydist = factor(ydist,levels = c('rnorm','runif','rlaplace','rbeta0.1','rbeta_lowskew','rbeta_medskew','rbeta_highskew','rlnorm0.5','rlnorm1','rlnorm1.5','rlnorm2'),
												labels = c('Normal','Uniform','Laplace','Beta(0.1. 0.1)','Beta (5,2)','Beta (5,1)','Beta (5,0.5)','Lognormal (0, 0.5)','Lognormal (0, 1)','Lognormal (0, 1.5)','Lognormal (0, 2)')),
				 x1dist = factor(x1dist,levels = c('rnorm','runif','rlaplace','rbeta0.1','rbeta_lowskew','rbeta_medskew','rbeta_highskew','rlnorm0.5','rlnorm1','rlnorm1.5','rlnorm2'),
				 								labels = c('Normal','Uniform','Laplace','Beta(0.1. 0.1)','Beta (5,2)','Beta (5,1)','Beta (5,0.5)','Lognormal (0, 0.5)','Lognormal (0, 1)','Lognormal (0, 1.5)','Lognormal (0, 2)'))) %>%
	pivot_longer(b1_ks) %>%
	ggplot(aes(x=n,y=value,color=x1dist))+
	geom_line(lwd = 1)+
	geom_vline(xintercept = 30,linetype = 2,lwd=0.25)+
	scale_x_continuous(trans='log', breaks = c(4,8,15,30,60,120,250,500,1000,2000,4000,10000))+
	scale_color_manual(name = 'X Distribution', values = c(wesanderson::wes_palette('Darjeeling1')[1:6],wesanderson::wes_palette('Darjeeling2')[1:5]))+
	facet_wrap(~ydist,nrow=2)+
	theme_minimal()+
	xlab('Sample Size')+
	ylab('Proportion of rejected Kolmogorov-Smirnoff test')+
	theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
				axis.title = element_text(size = 14, face = 'bold'),
				strip.text = element_text(size = 11, face = 'bold'),
				legend.title = element_text(size = 11, face = 'bold'),
				legend.text = element_text(size = 11))




