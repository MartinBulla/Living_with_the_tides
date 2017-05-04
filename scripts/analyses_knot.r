##################################### info #################################################
## R Script generates statistical outputs for the red knot case study from
## Bulla, M. et al. (2017). Marine biorhythms: bridging chronobiology and ecology. Phylo Trans B. 

# --------------------------------------------------------------------------------------------------------
# Questions can be directed to: Martin Bulla (bulla.mar@gmail.com) and Thomas Oudman (thomas.oudman@gmail.com)
# --------------------------------------------------------------------------------------------------------
# The script runs with 	- R version 3.3.0 (2016-05-03); for packages version see the sessionInfo() at the end of the script
#                      	- following files available from https://osf.io/xby9t/
	#						- low tides.txt
	#						- high tides.txt
	#						- median positions.txt
	#                       - subShab.tif
#						- following folders in your directory (defined as outdir)
#							- supplementary_actograms
#							- models_per_individual
#							- ACF
#							- simulations
# To save time, analyses can be performed without 
#							- initial data preparation using data.txt and 95%distance.txt files from 'data' folder at https://osf.io/xby9t/
#							- generating the posteriory distrubutions and model predictions using files from 'simulartions' folder at https://osf.io/xby9t/

{############################# load/define tools #################################
### for plotting
rm(list=ls()) # remove any objects from R
#palette("default")
PNG = TRUE # output is png file if FALSE prints plot within opened device in R

### set working directory
wd = "C:/Users/mbulla/Documents/Dropbox/Science/Projects/MS/Living_with_the_tides/Living_with_the_tides_shorebirds_analyses/"
outdir = "C:/Users/mbulla/Documents/Dropbox/Science/Projects/MS/Living_with_the_tides/Output/"
setwd("C:/Users/mbulla/Documents/Dropbox/Science/Projects/MS/Living_with_the_tides/Living_with_the_tides_shorebirds_analyses/")


### skip data generation and plots for each individual (TURE or FALSE), as generating data takes some minutes
skip = TRUE

### skip simulations (TRUE or FALSE), as generating simulations takes ~1h
skip2 = TRUE

### packages
library(raster)
library(fields)
library(RColorBrewer)
library(lme4)
library(arm)
library(effects)
library(multcomp)
library(plyr)
library(ggplot2)

### define colors for plotting
mypalette<-brewer.pal(7,"Blues")
mypalette <- mypalette[2:7]
col_f = '#FCB42C'
col_m = '#535F7C'

### get map
map <- raster(paste(wd,"data/bda.tif",sep=""))

### define sunrises and sunsets for plotting (sunrise/sunset 9 January 7:38/18:45, 13 feb 7:32/19:04 LOCAL TIME)
rx <- c(7+38/60,7+32/60)  #sunrise
sx <- c(18+45/60,19+4/60) #sunset
y <- c(36,1) #the day number for each tide

### sampling function (x = vector to sample from, size = how many point to sample, lag = how far apart shall be the closest points
	ssample = function(x, size, lag) { 
 
	r = data.frame(o = sort(sample(x, size)) )
	r$l = c(r$o[-1] - r$o[-nrow(r)], lag+1)
	r = r[r$l> lag, ]
 
	check = 0
 
	while(nrow(r) < size) {
		check = check + 1
 
		if(check > 10000) stop('lag  and/or size are probably too large')
 
		rj = data.frame(o = sort(sample(x, size)) )
		rj$l = c(rj$o[-1] - rj$o[-nrow(rj)], lag+1)
		r = rbind(r, rj)
		r = r[order(r$o), ]
		r$l = c(r$o[-1] - r$o[-nrow(r)], lag+1)
		r = r[r$l> lag, ]
 
		r
 
 
	}
 
	return(r$o[1:size])
 
 
 }

 }

############################# prepare data for analyses #########################################
 if(skip == FALSE{	
### high and low tide time
low <- read.table("data/low tides.txt",header=T,sep="\t")
high <- read.table("data/high tides.txt",header=T,sep="\t")
low$Timestamp.loc <- as.POSIXct(strptime(low$Timestamp.loc, format =  "%d-%m-%Y %H:%M"))
high$Timestamp.loc <- as.POSIXct(strptime(high$Timestamp.loc, format =  "%d-%m-%Y %H:%M"))

### define main high tide roosts
roosts <- data.frame(X=c(364306.2,362906.2,364631.2,366131.2,367331.2,367506.2,362306.2),
                     Y=c(2197575,2199150,2199125,2201550,2200425,2196075,2196275),ID=1:7)



### get bird locations
tabo <- read.table("data/median positions.txt",header=T,sep="\t")
tabo$Timestamp.loc <- as.POSIXct(strptime(tabo$Timestamp.loc, format =  "%d-%m-%Y %H:%M"))

### calculate distance to the closest roost for each bird
tags <- unique(tabo$Tag)
l = list() # empty list to save the 95% largest distance for each bird
for (t in tags){
  print(t)
  
  #### load positions
  tab1 <- tabo[tabo$Tag==t,]

  # calculate distance to nearest roost, and which roost that is
  tab1$Dist  <- NA
  tab1$Roost <- NA
  
  
  for(q in 1:dim(tab1)[1]){
    tab1$Dist[q] <- min(sqrt((tab1$x[q] - roosts$X)^2 + (tab1$y[q] - roosts$Y)^2))
    tab1$Roost[q] <- which(sqrt((tab1$x[q] - roosts$X)^2 + (tab1$y[q] - roosts$Y)^2)==tab1$Dist[q])
  }
  
  # standardize distance
  temp <- tab1$Dist[order(tab1$Dist)]
  maxd <- temp[round(length(temp)*0.95)]
  tab1$RDist <- tab1$Dist/maxd                     # distance divided by 95% largest distance
  tab1$RDist2 <- ifelse(tab1$RDist>1,1,tab1$RDist) # same but floored at 1, just for plotting purposes
  tab1$lRDist <- log10(tab1$Dist)/log10(maxd)      # log-transformed distance divided by log-transformed 95% largest distance
  
  
  # make new table including all calculated values
  if(t!=1){tabx <- rbind(tabx,tab1)}else{tabx <- tab1}
  
  # add the 95% largest distance for each bird
  l[[t]] = data.frame(tag = tabo$Tag[tabo$Tag==t][1], maxd = maxd, stringsAsFactors = FALSE) 
}
md = do.call(rbind, l)

### define explanatory variables 
# time from high tide
tabx$htide    <- NA
tabx$TideTime <- NA
tabx$nTideTime   <- NA
for(i in 1:dim(tabx)[1]){
  #i=10
  
  tabx$TideTime[i] <- min(abs(difftime(high$Timestamp.loc,tabx$Timestamp.loc[i],units="hours")))
  tabx$htide[i] <- which(abs(difftime(high$Timestamp.loc,tabx$Timestamp.loc[i],units="hours"))==tabx$TideTime[i])
  tabx$nTideTime[i] <- difftime(tabx$Timestamp.loc[i],high$Timestamp.loc[tabx$htide[i]],units="hours")
  
}
d=tabx
### time in radians
	d$rad_dn=2*pi*d$rhour / 24  
	d$rad_t=2*pi*d$nTideTime / 12.4

### keep only standardized distances smaller than 2
	d = d[d$RDist<2,]
	
### save the data-sets for later	
write.table(d,paste(wd,"data/data.txt",sep=''),row.names=F,sep="\t")
write.table(md,paste(wd,"data/95%distance.txt",sep=''),row.names=F,sep="\t")
#save(md, d, file = 'data/data.Rdata')	
}

############################# ANALYSES ###################################
if(skip2 == FALSE){
{# RUN FIRST - generate model output and model predictions for each bird			
	d = read.table(paste(wd,"data/data.txt",sep=''),header=T,sep="\t", stringsAsFactors = FALSE)
	d$Timestamp.loc <- as.POSIXct(strptime(d$Timestamp.loc, format =  "%d-%m-%Y %H:%M"))
	md = read.table(paste(wd,"data/95%distance.txt",sep=''),header=T,sep="\t", stringsAsFactors = FALSE)
	
	load(file = 'data/data.Rdata')		
	l = list()
	d$rad_t=2*pi*d$nTideTime / 12.4
	for(i in 1:length(unique(d$Tag))){
			print(i)
			di = d[d$Tag==unique(d$Tag)[i],]
			if(nrow(di)>100){ # include only nests with more than 50 h of date
					
				{# run the model and print rough model output and temporal auto-correlations
						di$t_=1:nrow(di)
						for(k in 1:100){
							lp = list()
							foo=ssample(di$t_, round(nrow(di)*0.15),2)	
							di_=di[di$t_%in%foo,]		
							m <- lm(RDist ~ sin(rad_dn) + cos(rad_dn) + sin(rad_t)+cos(rad_t), data=di_)
							#x = summary(glht(m))
							nsim <- 10000
							bsim <- sim(m, n.sim=nsim)
							lp[[k]]=data.frame(bsim@coef)
						}
						xx = do.call(rbind, lp)
						v = apply(xx, 2, quantile, prob=c(0.5))
						lwr = apply(xx, 2, quantile, prob=c(0.025))
						upr = apply(xx, 2, quantile, prob=c(0.975))
										
						l[[i]] = data.frame(tag = i, n = nrow(di), est = v,lwr = lwr,upr = upr, period = c(NA, 24,24,12.4,12.4), vari = c('int','sin','cos', 'sin','cos'))
						
						#png(paste(outdir,"models_per_individual/Thight_no_log", di_$Tag[1], "sampled.png", sep=""), width=8,height=8,units="in",res=600) 
						#plot(allEffects(m))
						#dev.off()
						
						#png(paste(outdir,"ACF/Thight_no_log", di_$Tag[1], "sampled.png", sep=""), width=8,height=8,units="in",res=600) 
						#acf(resid(m), type="p", main=list("Temporal autocorrelation:\npartial series residual",cex=0.8))
										
						#dev.off()
					}		
				{# predicted values for day				
				# values to predict for		
					newD=data.frame(rhour = seq(0,23.9, length.out=300),
									rad_t=mean(di_$rad_t)
									)
					newD$rad_dn = 2*pi*newD$rhour / 24  
					#newD$rad_t = 2*pi*newD$rhour / 12.4  
				# exactly the model which was used has to be specified here 
					X <- model.matrix(~ sin(rad_dn) + cos(rad_dn) + sin(rad_t)+cos(rad_t),data=newD)	
								
				# calculate predicted values and creditability intervals
					newD$pred <-(X%*%v) 
							predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
							for(j in 1:nsim) predmatrix[,j] <- (X%*%bsim@coef[j,])
							newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
							newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
					pp=newD	
					
					pp$tag=di$Tag[1]
					if(i==1){ write.table(pp, "simulations/fits_day_100s.txt", row.names = FALSE, col.names=TRUE)}else{ write.table(pp, "simulations/fits_day_100s.txt", append = TRUE, row.names = FALSE,col.names=FALSE)	}	
					}
				{# predicted values for tide			
				# values to predict for		
					newD=data.frame(nTideTime = seq(-6.2,6.2, length.out=300),
									rad_dn=mean(di_$rad_dn)
									)
					#newD$rad_dn = 2*pi*newD$rhour / 24  
					newD$rad_t = 2*pi*newD$nTideTime / 12.4  
				# exactly the model which was used has to be specified here 
					X <- model.matrix(~ sin(rad_dn) + cos(rad_dn) + sin(rad_t)+cos(rad_t),data=newD)	
								
				# calculate predicted values and creditability intervals
					newD$pred <-(X%*%v) 
							predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
							for(j in 1:nsim) predmatrix[,j] <- (X%*%bsim@coef[j,])
							newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
							newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
					pp=newD	
					
					pp$tag=di$Tag[1]
					if(i==1){ write.table(pp, "simulations/fits_tide_100s.txt", row.names = FALSE, col.names=TRUE)}else{ write.table(pp, "simulations/fits_tide_100s.txt", append = TRUE, row.names = FALSE,col.names=FALSE)	}
					}
				}
			}
	x = do.call(rbind,l) # modelu summaries
	#save(x, file = paste(wd, 'simulations/model_output_100.Rdata', sep=''))	
	write.table(x,paste(wd,"simulations/model_output_100.txt",sep=''),row.names=F,sep="\t")	
}								
}
				
{# proportion of individuals with given rhythms	
				x = read.table(paste(wd,"simulations/model_output_100.txt",sep=''),header=T,sep="\t", stringsAsFactors = FALSE)	#load(file = 'simulations/model_output_100.Rdata')		
				x=x[!x$vari=='int',]
				y=x[-which(x$lwr<0 & x$upr>0),] # only those with intervals not overlapping zero
				
				length(unique(y$tag[y$period==12.4]))/length(unique(x$tag)) # proportion of birds with clear tidal pattern
				length(unique(y$tag[y$period==24]))/length(unique(x$tag)) # proportion of birds with clear circadian pattern
				yy = y[y$period==12.4,]
				yy=unique(yy$tag[yy$tag%in%unique(y$tag[y$period==24])])
				length(yy)/length(unique(x$tag)) # proportion of birds with both patterns
				length(unique(x$tag)) # number of individuals
				
				# days of observation per tag
					d = read.table(paste(wd,"data/data.txt",sep=''),header=T,sep="\t", stringsAsFactors = FALSE)	
					n = ddply(d[d$Tag%in%unique(x$tag),],.(Tag), summarize, n = length(Tag)/48) 
					summary(n)
}				
{# check birds with reversed tidal rhythms
			pe = read.table('simulations/fits_tide_100s.txt', header = TRUE, stringsAsFactors = FALSE 
			pe_ = ddply(pe,.(tag), summarise, time_ = nTideTime[pred == max(pred)]) 
			pe_[pe_$time_>(-3.1) & pe_$time_<3.1,]
			table(d$Tag) 
			
			px = pe[pe$tag%in%pe_$tag[pe_$time_>(-3.1) & pe_$time_<3.1],]
			px=px[order(px$tag, px$nTideTime),]
			ggplot(px, aes(x = nTideTime, y = pred, col = as.factor(tag)))+geom_point()		
}	
{# Figure 2	- km scale
	{# run first
				md = read.table(paste(wd,"data/95%distance.txt",sep=''),header=T,sep="\t", stringsAsFactors = FALSE)
				x = read.table(paste(wd,"simulations/model_output_100.txt",sep=''),header=T,sep="\t", stringsAsFactors = FALSE)	#load(file = 'simulations/model_output_100.Rdata')		
				pd = read.table('simulations/fits_day_100s.txt', header = TRUE, stringsAsFactors = FALSE )
				pe = read.table('simulations/fits_tide_100s.txt', header = TRUE, stringsAsFactors = FALSE )
	}	
	{# tide - plot on log or standardized scale
				if(PNG == TRUE) {
					png(paste(outdir,"Figure_1A_orig_100.png", sep=""), width=1.85+0.3,height=1.5,units="in",res=600) 
					}else{
					dev.new(width=1.85+0.3,height=1.5)
					}	
				par(mar=c(0.8,0.2,0.2,2.5),oma = c(1, 2, 0, 0),ps=12, mgp=c(1.2,0.35,0), las=1, cex=1, col.axis="black",font.main = 1, col.lab="black", col.main="black", fg="black", cex.lab=0.5,cex.main=0.7, cex.axis=0.5, tcl=-0.1,bty="n",xpd=TRUE) #
			
						
				plot(pred ~ nTideTime , data = pe, 
										#ylab =NULL, 
										xaxt='n',
										#yaxt='n',
										#xaxs = 'i',
										#yaxs = 'i',
										ylim=c(0,2.5),
										#xlim=c(0,12.4),
										xlim=c(-6.2,6.2),
										type='n'
										) # col=z_g$cols, border=z_g$cols
										
				axis(1, at=c(-6.2,-3.1,0,3.1,6.2), label=c(-6.2,3.1,0,3.1,6.2), mgp=c(0,-0.20,0), lwd=0.5)
				#axis(1, at=c(0,6.2,12,18,24), label=c(0,'',12,'',24), mgp=c(0,-0.20,0))
				axis(2, at=seq(0,2.5, by=0.5), label=c('0.0','0.5','1.0','1.5','2','2.5'), lwd = 0.5)
				
				mtext("Time to/from high tide [h]",side=1,line=0.3, cex=0.5, las=1, col='black')
				mtext("Distance from the roost [km]",side=2,line=1, cex=0.5, las=3, col='black')
				
				text(x=6.2, y=2.5 ,"A", font = 2, cex =0.6 )
				
				for(i in 1:length(unique(pe$tag))){ # check tag 4, 7, 8, 11, 16
							pdi = pe[pe$tag==unique(pe$tag)[i],]
							pdi=pdi[order(pdi$nTideTime),]
							mdi = md$maxd[md$tag==unique(md$tag)[i]]/1000
							lines(pdi$nTideTime,pdi$pred*mdi, col = 'black')
						}
				
				
				
				 if(PNG == TRUE) {dev.off()}
			}
	{# 24h - plot on log or standardized scale
				if(PNG == TRUE) {
					png(paste(outdir,"Figure_1B_ori_100.png", sep=""), width=1.85+0.3,height=1.5,units="in",res=600) 
					}else{
					dev.new(width=1.85+0.3,height=1.5)
					}	
				par(mar=c(0.8,0.2,0.2,2.5),oma = c(1, 2, 0, 0),ps=12, mgp=c(1.2,0.35,0), las=1, cex=1, col.axis="black",font.main = 1, col.lab="black", col.main="black", fg="black", cex.lab=0.5,cex.main=0.7, cex.axis=0.5, tcl=-0.1,bty="n",xpd=TRUE) #
			
						
				plot(pred ~ rhour , data = pd, 
										#ylab =NULL, 
										xaxt='n',
										yaxt='n',
										#xaxs = 'i',
										#yaxs = 'i',
										ylim=c(0,2.5),
										#xlim=c(0,12.4),
										xlim=c(0,24),
										type='n'
										) # col=z_g$cols, border=z_g$cols
										
				#axis(1, at=c(0,6.2,12.4), label=c(0,6.2,12.4), mgp=c(0,-0.20,0))
				axis(1, at=c(0,6,12,18,24), label=c(0,'',12,'',24), mgp=c(0,-0.20,0), lwd = 0.5)
				#axis(2, at=c(0,1,2,3), label=c(0,1,2,3))
				
				mtext("Time of day [h]",side=1,line=0.3, cex=0.5, las=1, col='black')
				#mtext("Distance from the roost\n[standardized]",side=2,line=1, cex=0.6, las=3, col='grey30')
				
				text(x=24, y=2.5, "B", font = 2, cex =0.6 )
				text(x=12, y=2.5*1.1/1.2, "Farthest at",  cex =0.5 )
				text(x=10, y=2.5*1/1.2, "night", col = col_m, cex =0.5 )
				text(x=15, y=2.5*1/1.2, "day", col = col_f, cex =0.5 )
				
				for(i in 1:length(unique(pd$tag))){
							pdi = pd[pd$tag==unique(pd$tag)[i],]
							pdi=pdi[order(pdi$rhour),]
							if(pdi$rhour[pdi$pred == min(pdi$pred)]>19.5 | pdi$rhour[pdi$pred == min(pdi$pred)]<7){col_=col_f}else{col_=col_m}
							mdi = md$maxd[md$tag==unique(md$tag)[i]]/1000
							lines(pdi$rhour,pdi$pred*mdi, col = col_)
						}
						
				
				
				 if(PNG == TRUE) {dev.off()}
			}
}					
	{# not in paper - Figure 2 - on standardized sclae
		{# run first
					md = read.table(paste(wd,"data/95%distance.txt",sep=''),header=T,sep="\t", stringsAsFactors = FALSE)
					x = read.table(paste(wd,"simulations/model_output_100.txt",sep=''),header=T,sep="\t", stringsAsFactors = FALSE)	#load(file = 'simulations/model_output_100.Rdata')		
					pd = read.table('simulations/fits_day_100s.txt', header = TRUE, stringsAsFactors = FALSE )
					pe = read.table('simulations/fits_tide_100s.txt', header = TRUE, stringsAsFactors = FALSE)
		}	
		{# tide - plot on log or standardized scale
					if(PNG == TRUE) {
						png(paste(outdir,"Figure_1A_100.png", sep=""), width=1.85+0.3,height=1.5,units="in",res=600) 
						}else{
						dev.new(width=1.85+0.3,height=1.5)
						}	
					par(mar=c(0.8,0.2,0.2,2.5),oma = c(1, 2, 0, 0),ps=12, mgp=c(1.2,0.35,0), las=1, cex=1, col.axis="grey30",font.main = 1, col.lab="grey30", col.main="grey30", fg="grey40", cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.1,bty="n",xpd=TRUE) #
				
							
					plot(pred ~ nTideTime , data = pe, 
											#ylab =NULL, 
											xaxt='n',
											#yaxt='n',
											#xaxs = 'i',
											#yaxs = 'i',
											ylim=c(0,1.2),
											#xlim=c(0,12.4),
											xlim=c(-6.2,6.2),
											type='n'
											) # col=z_g$cols, border=z_g$cols
											
					axis(1, at=c(-6.2,-3.1,0,3.1,6.2), label=c(-6.2,3.1,0,3.1,6.2), mgp=c(0,-0.20,0))
					#axis(1, at=c(0,6.2,12,18,24), label=c(0,'',12,'',24), mgp=c(0,-0.20,0))
					#axis(2, at=seq(-3,0.5, by=0.5), label=T)
					
					mtext("Time from high tide [h]",side=1,line=0.3, cex=0.6, las=1, col='grey30')
					mtext("Distance from the roost\n[standardized]",side=2,line=1, cex=0.6, las=3, col='grey30')
					
					text(x=6.2, y=1.2, "A", font = 2, cex =0.6 )
					
					for(i in 1:length(unique(pe$tag))){ # check tag 4, 7, 8, 11, 16
								pdi = pe[pe$tag==unique(pe$tag)[i],]
								pdi=pdi[order(pdi$nTideTime),]
								lines(pdi$nTideTime,pdi$pred)
							}
					
					
					
					 if(PNG == TRUE) {dev.off()}
				}
		{# 24h - plot on log or standardized scale
					if(PNG == TRUE) {
						png(paste(outdir,"Figure_1B_100.png", sep=""), width=1.85+0.3,height=1.5,units="in",res=600) 
						}else{
						dev.new(width=1.85+0.3,height=1.5)
						}	
					par(mar=c(0.8,0.2,0.2,2.5),oma = c(1, 2, 0, 0),ps=12, mgp=c(1.2,0.35,0), las=1, cex=1, col.axis="grey30",font.main = 1, col.lab="grey30", col.main="grey30", fg="grey40", cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.1,bty="n",xpd=TRUE) #
				
							
					plot(pred ~ rhour , data = pd, 
											#ylab =NULL, 
											xaxt='n',
											yaxt='n',
											#xaxs = 'i',
											#yaxs = 'i',
											ylim=c(0,1.2),
											#xlim=c(0,12.4),
											xlim=c(0,24),
											type='n'
											) # col=z_g$cols, border=z_g$cols
											
					#axis(1, at=c(0,6.2,12.4), label=c(0,6.2,12.4), mgp=c(0,-0.20,0))
					axis(1, at=c(0,6,12,18,24), label=c(0,'',12,'',24), mgp=c(0,-0.20,0))
					#axis(2, at=c(0,1,2,3), label=c(0,1,2,3))
					
					mtext("Time of day [h]",side=1,line=0.3, cex=0.6, las=1, col='grey30')
					#mtext("Distance from the roost\n[standardized]",side=2,line=1, cex=0.6, las=3, col='grey30')
					text(x=24, y=1.2, "B", font = 2, cex =0.6 )
					text(x=12, y=1.1, "Farthest at",  cex =0.5 )
					text(x=10, y=1, "night", col = col_m, cex =0.5 )
					text(x=15, y=1, "day", col = col_f, cex =0.5 )
					
					for(i in 1:length(unique(pd$tag))){
								pdi = pd[pd$tag==unique(pd$tag)[i],]
								pdi=pdi[order(pdi$rhour),]
								if(pdi$rhour[pdi$pred == min(pdi$pred)]>19.5 | pdi$rhour[pdi$pred == min(pdi$pred)]<7){col_=col_f}else{col_=col_m}
								lines(pdi$rhour,pdi$pred, col = col_)
							}
							
					
					
					 if(PNG == TRUE) {dev.off()}
				}
		}										
{# Figure 3
 ### load distance data
	 d = read.table(paste(wd,"data/data.txt",sep=''),header=T,sep="\t", stringsAsFactors = FALSE)
	 md = read.table(paste(wd,"data/95%distance.txt",sep=''),header=T,sep="\t", stringsAsFactors = FALSE)
	 t = 1 # tag 1
	 tab1 = d[d$Tag==t,]
	 maxd = md$maxd[md$tag==t]
 ### prepare double plot
  tab2 <- tab1
  #head(tab2)
  tab2$rhour <- tab1$rhour+24
  tab2$day <- tab1$day+1
  tab3 <- rbind(tab1,tab2)

  
 ### prepare low tide and day-night data for B panel
	low <- read.table("data/low tides.txt",header=T,sep="\t")
	low$Timestamp.loc <- as.POSIXct(strptime(low$Timestamp.loc, format =  "%d-%m-%Y %H:%M"))
	low <- low[low$Date!="02-14",]
	low$id <- 1:dim(low)[1]
	low$sun <- ifelse(low$hour>7.5 & low$hour<19.75,0,1) #night = 1, day = 0

  ### define main high tide roosts
	roosts <- data.frame(X=c(364306.2,362906.2,364631.2,366131.2,367331.2,367506.2,362306.2),
                     Y=c(2197575,2199150,2199125,2201550,2200425,2196075,2196275),ID=1:7)

  #### make plot
  if(PNG == TRUE) {
					png(paste(outdir,"Figure 3.png", sep=""),width=14, height=8, units="cm", pointsize = 8, res = 1000)
					}else{
					dev.new(width=5.51181, height=3.14961, pointsize = 8)
					}	
  split.screen(rbind(c(0,0.6,0,1), c(0.6,1,0,1)))
  screen(1)
  
  par(mar=c(4,5,4,0))
  quilt.plot(y=tab3$day,x=tab3$rhour,z=tab3$RDist2,
             nlevel=6,
             nx=96,
             ny=max(tab1$day)-min(tab1$day)+1,
             xlim=c(0,48),
             ylim=c(0.5,36.5),
             col=mypalette,
             axes=F,
             xlab="Time",
             add.legend=F
  )
  axis(2,at=c(4,15,25,35),labels=c("10 Feb","30 Jan","20 Jan","10 Jan"),las=1)
  abline(lm(y~rx),lwd=2) # sunrise
  abline(lm(y~sx),lwd=2)  # sunset
  rx2 <- rx+24
  sx2 <- sx+24
  abline(lm(y~rx2),lwd=2) # sunrise
  abline(lm(y~sx2),lwd=2)  # sunset
  points(day~hour,data=low[low$Number==1,],type="l",lwd=2,lty=2)
  points(day~hour,data=low[low$Number==2,],type="l",lwd=2,lty=2)
  points(day~hour,data=low[low$Number==3,],type="l",lwd=2,lty=2)
  points(day~hour,data=low[low$Number==4,],type="l",lwd=2,lty=2)
  low$day2 <- low$day+1
  low$hour2 <- low$hour+24
  points(day2~hour2,data=low[low$Number==1,],type="l",lwd=2,lty=2)
  points(day2~hour2,data=low[low$Number==2,],type="l",lwd=2,lty=2)
  points(day2~hour2,data=low[low$Number==3,],type="l",lwd=2,lty=2)
  points(day2~hour2,data=low[low$Number==4,],type="l",lwd=2,lty=2)
  abline(h=seq(1.5,40.5,1))
  box()
  points(44,34.1,pch=21,cex=4,bg=0)
  text(x=44,y=34.1,labels=c("A"),cex=1.2)
  par(xpd=T)
  axis(1,at=seq(0,50,12),labels=c("0:00","12:00","24:00","36:00","48:00"))
  polygon(x=c(0,48,48,0),y=c(37.5,37.5,38.5,38.5),col=0)
  polygon(x=c(0,7.6,7.6,0),y=c(37.5,37.5,38.5,38.5),col=1)
  polygon(x=c(18.7,31.6,31.6,18.7),y=c(37.5,37.5,38.5,38.5),col=1)
  polygon(x=c(42.8,48,48,42.8),y=c(37.5,37.5,38.5,38.5),col=1)
  text(x=-1,y=38,labels="dark/light",pos=2)
  x <- round(seq(0.1,0.9,0.2)*maxd/100)*100
  par(xpd=F)
  
  screen(2)
  par(mar=c(4,0,4,9))
  plot(rev(id)~Height,data=low,ylim=c(2.5,68.5),pch=21,bg=sun,
       xlab="",
       ylab="",
       axes=F
       )
  points(0.65,66,pch=21,cex=4,bg=0)
  text(x=0.65,y=66,labels=c("B"),cex=1.2)
  axis(1,at=c(0.17,0.4,0.6),labels=c("48:00",0.4,0.6))
  #axis(2,at=c(7,29,49,68),labels=c("","","",""),las=1)
  box()
  par(xpd=T)
  mtext(side=1,"Water height (m)",line=3)
 legend("topright",inset=c(-1.13,0),
         fill=c(0,mypalette),
         legend=c("no data",paste("<",x[1]," m",sep=""),
                  paste(x[1]," - ",x[2]," m",sep=""),
                  paste(x[2]," - ",x[3]," m",sep=""),
                  paste(x[3]," - ",x[4]," m",sep=""),
                  paste(x[4]," - ",x[5]," m",sep=""),
                  paste(">",x[5]," m",sep="")),
         bty="n")
  legend("topright",inset=c(-1.13,0.45),
         lwd=c(2,2,NA,NA),
         lty=c(1,2,NA,NA),
         pt.bg=c(NA,NA,"white","black"),
         pch=c(NA,NA,21,21),
         legend=c("sunrise/sunset","low tide","day low tide","night low tide"),bty="n")
  par(xpd=F)
  close.screen(all.screens = TRUE)
  if(PNG == TRUE) {dev.off()}
}
{# Supplementary actograms
 ### load distance data
	 d = read.table(paste(wd,"data/data.txt",sep=''),header=T,sep="\t", stringsAsFactors = FALSE)
		# use only birds with more than 50h of data
		dx = ddply(d,.(Tag), summarise, n = length(Tag))
		d = d[d$Tag%in%dx$Tag[dx$n>100],] 
	 md = read.table(paste(wd,"data/95%distance.txt",sep=''),header=T,sep="\t", stringsAsFactors = FALSE)
	 tags <- unique(d$Tag)
 ### load low tide data	
	low <- read.table("data/low tides.txt",header=T,sep="\t")
	low$Timestamp.loc <- as.POSIXct(strptime(low$Timestamp.loc, format =  "%d-%m-%Y %H:%M"))
	
  ### define main high tide roosts
	roosts <- data.frame(X=c(364306.2,362906.2,364631.2,366131.2,367331.2,367506.2,362306.2),
                     Y=c(2197575,2199150,2199125,2201550,2200425,2196075,2196275),ID=1:7)
	
 #### loop per bird for plotting purposes
 for (t in tags){
  #t=1
  print(t)
  
  #### load positions + prepare double plot
  tab1 <- d[d$Tag==t,]
  tab2 <- tab1
  #head(tab2)
  tab2$rhour <- tab1$rhour+24
  tab2$day <- tab1$day+1
  tab3 <- rbind(tab1,tab2)
  
  ### 95% distance 
	maxd = md$maxd[md$tag==t]
  
  ### plot
	tit = if(nchar(t) == 1){paste('0',t,sep='')}else{t}
	png(paste(outdir,"supplementary_actograms/tag_",tit,".png", sep=""), width=12, height=8, units="cm", pointsize = 8, res = 1000)
	  par(mar=c(4,4,4,9))
	  quilt.plot(y=tab3$day,x=tab3$rhour,z=tab3$RDist2,
				 nlevel=6,
				 nx=96,
				 ny=max(tab1$day)-min(tab1$day)+1,
				 xlim=c(0,48),
				 ylim=c(0.5,36.5),
				 col=mypalette,
				 axes=F,
				 add.legend=F
	  )
	  title(main=paste("Distance to nearest roost tag ",t,sep="")) #for some strange reason axis doesn't work without this
	  axis(1,at=seq(0,50,12),labels=c("0:00","12:00","24:00","36:00","48:00"))
	  axis(2,at=c(4,15,25,35),labels=c("10 Feb","30 Jan","20 Jan","10 Jan"),las=1)
	  mtext(side=1,"Time [h]",line=2)
	  abline(lm(y~rx),lwd=2) # sunrise
	  abline(lm(y~sx),lwd=2)  # sunset
	  rx2 <- rx+24
	  sx2 <- sx+24
	  abline(lm(y~rx2),lwd=2) # sunrise
	  abline(lm(y~sx2),lwd=2)  # sunset
	  points(day~hour,data=low[low$Number==1,],type="l",lwd=2,lty=2)
	  points(day~hour,data=low[low$Number==2,],type="l",lwd=2,lty=2)
	  points(day~hour,data=low[low$Number==3,],type="l",lwd=2,lty=2)
	  points(day~hour,data=low[low$Number==4,],type="l",lwd=2,lty=2)
	  low$day2 <- low$day+1
	  low$hour2 <- low$hour+24
	  points(day2~hour2,data=low[low$Number==1,],type="l",lwd=2,lty=2)
	  points(day2~hour2,data=low[low$Number==2,],type="l",lwd=2,lty=2)
	  points(day2~hour2,data=low[low$Number==3,],type="l",lwd=2,lty=2)
	  points(day2~hour2,data=low[low$Number==4,],type="l",lwd=2,lty=2)
	  abline(h=seq(1.5,40.5,1))
	  box()
	  par(xpd=T)
	  polygon(x=c(0,48,48,0),y=c(37.5,37.5,38.5,38.5),col=0)
	  polygon(x=c(0,7.6,7.6,0),y=c(37.5,37.5,38.5,38.5),col=1)
	  polygon(x=c(18.7,31.6,31.6,18.7),y=c(37.5,37.5,38.5,38.5),col=1)
	  polygon(x=c(42.8,48,48,42.8),y=c(37.5,37.5,38.5,38.5),col=1)
	  text(x=50,y=38,labels="dark/light",pos=4)
	  x <- round(seq(0.1,0.9,0.2)*maxd/100)*100
	  legend("topright",inset=c(-0.38,0),
			 fill=c(0,mypalette),
			 legend=c("no data",paste("<",x[1]," m",sep=""),
					  paste(x[1]," - ",x[2]," m",sep=""),
					  paste(x[2]," - ",x[3]," m",sep=""),
					  paste(x[3]," - ",x[4]," m",sep=""),
					  paste(x[4]," - ",x[5]," m",sep=""),
					  paste(">",x[5]," m",sep="")),
			 bty="n")
	  legend("topright",inset=c(-0.38,0.5),lwd=1,lty=c(1,2),legend=c("sunrise/sunset","low tide"),bty="n")

	  par(xpd=F)
  dev.off()  

 
}

}

############################## sessionInfo() ##############################		
R version 3.3.0 (2016-05-03)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 7 x64 (build 7601) Service Pack 1

locale:
[1] LC_COLLATE=English_United States.1252 
[2] LC_CTYPE=English_United States.1252   
[3] LC_MONETARY=English_United States.1252
[4] LC_NUMERIC=C                          
[5] LC_TIME=English_United States.1252    

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] ggplot2_2.1.0      plyr_1.8.3         multcomp_1.4-5     TH.data_1.0-7     
 [5] survival_2.39-4    mvtnorm_1.0-5      effects_3.1-1      arm_1.8-6         
 [9] MASS_7.3-45        lme4_1.1-12        Matrix_1.2-6       RColorBrewer_1.1-2
[13] fields_8.4-1       maps_3.1.0         spam_1.3-0         raster_2.5-8      
[17] sp_1.2-3          

loaded via a namespace (and not attached):
 [1] Rcpp_0.12.5      splines_3.3.0    munsell_0.4.3    colorspace_1.2-6
 [5] lattice_0.20-33  minqa_1.2.4      rgdal_1.1-10     nnet_7.3-12     
 [9] gtable_0.2.0     nlme_3.1-128     coda_0.18-1      abind_1.4-3     
[13] nloptr_1.0.4     codetools_0.2-14 sandwich_2.3-4   scales_0.4.0    
[17] zoo_1.7-13
		