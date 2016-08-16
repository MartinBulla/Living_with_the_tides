	# define working directory
	wd="C:/Users/mbulla/Documents/Dropbox/Science/Projects/MS/Comparative/Tidal_rhythms//Data/"
	
	# load data - used in Bulla et al. 2016 and available here https://osf.io/eydx5/ 
	load(paste(wd,"comparative_all.Rdata", sep=""))
		zdn=z[!is.na(z$period),]	 # period dataset	
		
	# number of species
		length(unique(zdn$sp))
		
	# proporiton of populations foraging with tide in their breeding grounds
		length(unique(zdn$pop[which(zdn$tidal_pop=='y')]))/length(unique(zdn$pop[which(zdn$tidal_pop%in%c('y','n'))])) 
	
	# proporiton of nests with period entrainable by tide	
		zdn$tid_s=ifelse(is.na(zdn$period), NA, ifelse(zdn$period%in%c(24.75,12.5,6.25,3),1,0))
		summary(factor(zdn$tid_s))
		nrow(zdn[zdn$tid_s==1,])/nrow(zdn) 
	
	
	# proportion of nests with period entrainable by tide from populations known to used tide for foraging
		zdnt=zdn[which(zdn$tidal_pop=='y'),]
		nrow(zdnt[zdnt$tid_s==1,])/nrow(zdnt) 
		nrow(zdnt)
		length(unique(zdnt$pop))
		length(unique(zdnt$sp))
		
		length(zdnt$pop[zdnt$tid_s==1]) 
		length(zdnt$sp[zdnt$tid_s==1]) 
		zdnt$nest[zdnt$tid_s==1]
		