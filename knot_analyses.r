##################################### info #################################################
## R Script generates statistical outputs for the red knot case study from
## Bulla, M. et al. (2017). Marine biorhythms: bridging chronobiology and ecology. Phylo Trans B. 

# --------------------------------------------------------------------------------------------------------
# Questions can be directed to: Martin Bulla (bulla.mar@gmail.com) and Thomas Oudman (thomas.oudman@gmail.com)
# --------------------------------------------------------------------------------------------------------
# The script runs with 	- R version 3.3.0 (2014-07-10); for packages version see the sessionInfo() at the end of the script
#                      	- following files available from https://osf.io/xby9t/
#						- low tides.txt
#						- high tides.txt
#						- median positions.txt
#                       - subShab.tif

##################################### load tools ###########################################

### set working directory
setwd("C:/Users/mbulla/Documents/Dropbox/Science/Projects/MS/Living_with_the_tides/For Martin")
outdir = "C:/Users/mbulla/Documents/Dropbox/Science/Projects/MS/Living_with_the_tides/For Martin/"

### packages
library(raster)
library(fields)
library(RColorBrewer)
library(lme4)
library(AICcmodavg)
library(arm)
library(effects)
library(multcomp)
library(plyr)
library(lattice)
library(ggplot2)

### for plotting
rm(list=ls())
palette("default")
PNG = TRUE # output is png file if FALSE prints plot within opened device in R

### define colors for plotting
mypalette<-brewer.pal(7,"Blues")
mypalette <- mypalette[2:7]
col_f = '#FCB42C'
col_m = '#535F7C'

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
	

##################################### get data ###########################################

### high and low tide time
low <- read.table("low tides.txt",header=T,sep="\t")
high <- read.table("high tides.txt",header=T,sep="\t")
low$Timestamp.loc <- as.POSIXct(strptime(low$Timestamp.loc, format =  "%d-%m-%Y %H:%M"))
high$Timestamp.loc <- as.POSIXct(strptime(high$Timestamp.loc, format =  "%d-%m-%Y %H:%M"))

### define main high tide roosts
roosts <- data.frame(X=c(364306.2,362906.2,364631.2,366131.2,367331.2,367506.2,362306.2),
                     Y=c(2197575,2199150,2199125,2201550,2200425,2196075,2196275),ID=1:7)

### get map
map <- raster("sup5hab.tif")

### get bird locations
tabo <- read.table("median positions.txt",header=T,sep="\t")
tabo$Timestamp.loc <- as.POSIXct(strptime(tabo$Timestamp.loc, format =  "%d-%m-%Y %H:%M"))


############################# PREPARE DATA ###############################
### CALCULATE DISTANCE TO CLOSEST ROOST and PLOT
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
  
  #### make plot
  png(paste("Plots/Plots per bird/tag_",t,".png", sep=""), width=10.5, height=8, units="cm", pointsize = 8, res = 1000)
  par(mar=c(4,4,4,9))
  quilt.plot(y=tab1$day,x=tab1$rhour,z=tab1$RDist2,
             nlevel=6,
             nx=48,
             ny=max(tab1$day)-min(tab1$day)+1,
             ylim=c(0.5,36.5),
             col=mypalette,
             axes=F,
             add.legend=F
  )
  title(main=paste("Standardized distance from nearest roost tag ",t,sep="")) #for some strange reason axis doesn't work without this
  axis(1,at=seq(0,24,6),labels=c("0:00","6:00","12:00","18:00","24:00"))
  axis(2,at=c(4,15,25,35),labels=c("10 Feb","30 Jan","20 Jan","10 Jan"),las=1)
  abline(lm(y~rx),lwd=3) # sunrise
  abline(lm(y~sx),lwd=3)  # sunset
  points(day~hour,data=low[low$Number==1,],type="l",col=1,lwd=3,lty=2)
  points(day~hour,data=low[low$Number==2,],type="l",col=1,lwd=3,lty=2)
  points(day~hour,data=low[low$Number==3,],type="l",col=1,lwd=3,lty=2)
  points(day~hour,data=low[low$Number==4,],type="l",col=1,lwd=3,lty=2)
  abline(h=seq(1.5,40.5,1))
  box()
  par(xpd=T)
  legend("topright",inset=c(-0.47,0),
         fill=c(0,mypalette),
         legend=c("no data","< 0.1","0.1 - 0.3","0.3 - 0.5","0.5 - 0.7","0.7 - 0.9",">0.9"),bty="n", cex=0.8)
  legend("topright",inset=c(-0.47,0.45),
         lwd=2,
         lty=c(1,2),
         legend=c("sunrise/sunset","low tide"),bty="n", cex=0.8)
 legend("topright",inset=c(-0.47,0.65),
         legend=paste("95% furthest distance:\n", round(maxd/1000,1)," km", sep=""),bty="n", cex=0.8) 
  par(xpd=F)
  dev.off()
  
  # make new table including all calculated values
  if(t!=1){tabx <- rbind(tabx,tab1)}else{tabx <- tab1}
  
  # add the 95% largest distance for each bird
  l[[t]] = data.frame(tag = tabo$Tag[tabo$Tag==t][1], maxd = maxd, stringsAsFactors = FALSE) 
}
md = do.call(rbind, l)
save(md, tabx, file = 'data.Rdata')
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

### time in radians
d=tabx
d$rad_dn=2*pi*d$rhour / 24  
d$rad_t=2*pi*d$nTideTime / 12.4

### keep only standardized distances smaller than 2
	d = d[d$RDist<2,]

########################## ANALYSES #####################################

{# generated model output and model predicttions for each bird				
	l = list()
	d$rad_t=2*pi*d$nTideTime / 12.4
	for(i in 1:length(unique(d$Tag))){
			print(i)
			di = d[d$Tag==unique(d$Tag)[i],]
			if(nrow(di)>100){ # include only nests with more than 50 h of date
					
				{# run the model and print rough model output and temporal auto-correlations
						di$t_=1:nrow(di)
						foo=ssample(di$t_, round(nrow(di)*0.15),2)	
						di_=di[di$t_%in%foo,]		
						
						m <- lm(RDist ~ sin(rad_dn) + cos(rad_dn) + sin(rad_t)+cos(rad_t), data=di_)
						x = summary(glht(m))
						nsim <- 10000
						bsim <- sim(m, n.sim=nsim)
						v = apply(bsim@coef, 2, quantile, prob=c(0.5))
						lwr = apply(bsim@coef, 2, quantile, prob=c(0.025))
						upr = apply(bsim@coef, 2, quantile, prob=c(0.975))
										
						l[[i]] = data.frame(tag = i, n = nrow(di), est = v[2:5],lwr = lwr[2:5],upr = upr[2:5], period = c(24,24,12.4,12.4), vari = c('sin','cos', 'sin','cos'), p = x$test$pvalues[2:5])
						
						png(paste(outdir,"models_per_individual/Thight_no_log", di_$Tag[1], "sampled.png", sep=""), width=8,height=8,units="in",res=600) 
						plot(allEffects(m))
						dev.off()
						
						png(paste(outdir,"ACF/Thight_no_log", di_$Tag[1], "sampled.png", sep=""), width=8,height=8,units="in",res=600) 
						acf(resid(m), type="p", main=list("Temporal autocorrelation:\npartial series residual",cex=0.8))
										
						dev.off()
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
					if(i==1){ write.table(pp, "fits_day_htide_nolog.txt", row.names = FALSE, col.names=TRUE)}else{ write.table(pp, "fits_day_htide_nolog.txt", append = TRUE, row.names = FALSE,col.names=FALSE)	}	
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
					if(i==1){ write.table(pp, "fits_tide_hTide_nolog.txt", row.names = FALSE, col.names=TRUE)}else{ write.table(pp, "fits_tide_hTide_nolog.txt", append = TRUE, row.names = FALSE,col.names=FALSE)	}
					}
				}
			}
	x = do.call(rbind,l) # modelu summaries
	save(x, file = 'rhythms_htide_extremes_out_no_log.Rdata')			
				
				
}								
				
{# proportion of individuals with given rhythms						
				
				load(file = 'rhythms_htide_extremes_out.Rdata')
				pd = read.table('fits_day_htide_nolog.txt', header = TRUE, stringsAsFactors = FALSE )
				pe = read.table('fits_tide_hTide_nolog.txt', header = TRUE, stringsAsFactors = FALSE )
							
				y=x[-which(x$lwr<0 & x$upr>0),] # only those with intervals not overlapping zero
				length(unique(y$tag[y$period==12.4]))/length(unique(x$tag)) # proportion of birds with clear tidal pattern
				length(unique(y$tag[y$period==24]))/length(unique(x$tag)) # proportion of birds with clear circadian pattern
				yy = y[y$period==12.4,]
				yy=unique(yy$tag[yy$tag%in%unique(y$tag[y$period==24])])
				length(yy)/length(unique(x$tag)) # proportion of birds with both patterns
				length(unique(x$tag)) # number of individuals
				n = ddply(d[d$Tag%in%unique(x$tag),],.(Tag), summarize, n = length(Tag)/48)
				summary(n)
}				
{# check birds with reversed tidal rhythms
			pe_ = ddply(pe,.(tag), summarise, time_ = nTideTime[pred == max(pred)]) 
			pe_[pe_$time_>(-3.1) & pe_$time_<3.1,]
			table(d$Tag) 
			
			px = pe[pe$tag%in%pe_$tag[pe_$time_>(-3.1) & pe_$time_<3.1],]
			px=px[order(px$tag, px$nTideTime),]
			ggplot(px, aes(x = nTideTime, y = pred, col = as.factor(tag)))+geom_point()		
}	
{# Figure 2 - standardized sclae
		
	{# tide - plot on log or standardized scale
				if(PNG == TRUE) {
					png(paste(outdir,"Figure_1A_no_log.png", sep=""), width=1.85+0.3,height=1.5,units="in",res=600) 
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
					png(paste(outdir,"Figure_1B_col_inv.png", sep=""), width=1.85+0.3,height=1.5,units="in",res=600) 
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
{# Figure 2	- km scale
	{# tide - plot on log or standardized scale
				if(PNG == TRUE) {
					png(paste(outdir,"Figure_1A_orig.png", sep=""), width=1.85+0.3,height=1.5,units="in",res=600) 
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
										ylim=c(0,2.5),
										#xlim=c(0,12.4),
										xlim=c(-6.2,6.2),
										type='n'
										) # col=z_g$cols, border=z_g$cols
										
				axis(1, at=c(-6.2,-3.1,0,3.1,6.2), label=c(-6.2,3.1,0,3.1,6.2), mgp=c(0,-0.20,0))
				#axis(1, at=c(0,6.2,12,18,24), label=c(0,'',12,'',24), mgp=c(0,-0.20,0))
				#axis(2, at=seq(-3,0.5, by=0.5), label=T)
				
				mtext("Time from high tide [h]",side=1,line=0.3, cex=0.6, las=1, col='grey30')
				mtext("Distance from the roost\n[km]",side=2,line=1, cex=0.6, las=3, col='grey30')
				
				text(x=6.2, y=2.5 ,"A", font = 2, cex =0.6 )
				
				for(i in 1:length(unique(pe$tag))){ # check tag 4, 7, 8, 11, 16
							pdi = pe[pe$tag==unique(pe$tag)[i],]
							pdi=pdi[order(pdi$nTideTime),]
							mdi = md$maxd[md$tag==unique(md$tag)[i]][i]/1000
							lines(pdi$nTideTime,pdi$pred*mdi)
						}
				
				
				
				 if(PNG == TRUE) {dev.off()}
			}
	{# 24h - plot on log or standardized scale
				if(PNG == TRUE) {
					png(paste(outdir,"Figure_1B_ori.png", sep=""), width=1.85+0.3,height=1.5,units="in",res=600) 
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
										ylim=c(0,2.5),
										#xlim=c(0,12.4),
										xlim=c(0,24),
										type='n'
										) # col=z_g$cols, border=z_g$cols
										
				#axis(1, at=c(0,6.2,12.4), label=c(0,6.2,12.4), mgp=c(0,-0.20,0))
				axis(1, at=c(0,6,12,18,24), label=c(0,'',12,'',24), mgp=c(0,-0.20,0))
				#axis(2, at=c(0,1,2,3), label=c(0,1,2,3))
				
				mtext("Time of day [h]",side=1,line=0.3, cex=0.6, las=1, col='grey30')
				#mtext("Distance from the roost\n[standardized]",side=2,line=1, cex=0.6, las=3, col='grey30')
				
				text(x=24, y=2.5, "B", font = 2, cex =0.6 )
				text(x=12, y=2.5*1.1/1.2, "Farthest at",  cex =0.5 )
				text(x=10, y=2.5*1/1.2, "night", col = col_m, cex =0.5 )
				text(x=15, y=2.5*1/1.2, "day", col = col_f, cex =0.5 )
				
				for(i in 1:length(unique(pd$tag))){
							pdi = pd[pd$tag==unique(pd$tag)[i],]
							pdi=pdi[order(pdi$rhour),]
							if(pdi$rhour[pdi$pred == min(pdi$pred)]>19.5 | pdi$rhour[pdi$pred == min(pdi$pred)]<7){col_=col_f}else{col_=col_m}
							mdi = md$maxd[md$tag==unique(md$tag)[i]][i]/1000
							lines(pdi$rhour,pdi$pred*mdi, col = col_)
						}
						
				
				
				 if(PNG == TRUE) {dev.off()}
			}
}					

########################## sessionInfo() #####################################		
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
 [1] ggplot2_2.1.0      lattice_0.20-33    plyr_1.8.3         multcomp_1.4-5    
 [5] TH.data_1.0-7      survival_2.39-4    mvtnorm_1.0-5      effects_3.1-1     
 [9] arm_1.8-6          MASS_7.3-45        AICcmodavg_2.0-4   lme4_1.1-12       
[13] Matrix_1.2-6       RColorBrewer_1.1-2 fields_8.4-1       maps_3.1.0        
[17] spam_1.3-0         raster_2.5-8       sp_1.2-3          

loaded via a namespace (and not attached):
 [1] Rcpp_0.12.5      nloptr_1.0.4     digest_0.6.9     nlme_3.1-128    
 [5] gtable_0.2.0     rgdal_1.1-10     coda_0.18-1      stats4_3.3.0    
 [9] nnet_7.3-12      reshape_0.8.5    VGAM_1.0-2       minqa_1.2.4     
[13] scales_0.4.0     codetools_0.2-14 splines_3.3.0    abind_1.4-3     
[17] unmarked_0.11-0  xtable_1.8-2     colorspace_1.2-6 labeling_0.3    
[21] sandwich_2.3-4   munsell_0.4.3    zoo_1.7-13
		