# R Script generates statistical outputs for the biparental shorebird incubation case study from
# Bulla, M. et al. (2017). Marine biorhythms: bridging chronobiology and ecology. Phylo Trans B. 

	# define working directory
	wd="C:/Users/mbulla/Documents/Dropbox/Science/Projects/MS/Living_with_the_tides/Living_with_the_tides_shorebirds_analyses/"
	
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
		
	# proportion of nests with period entrainable by tide from populations known not to used tide for foraging
		zdnt=zdn[which(zdn$tidal_pop=='n'),]
		nrow(zdnt[zdnt$tid_s==1,])/nrow(zdnt) 
		nrow(zdnt)
		length(unique(zdnt$pop))
		length(unique(zdnt$sp))
		
		nrow(zdnt[zdnt$tid_s==1,])
		length(zdnt$pop[zdnt$tid_s==1]) 
		length(zdnt$sp[zdnt$tid_s==1]) 
		zdnt$nest[zdnt$tid_s==1]	

# sessionInfo()		
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
[1] stats     graphics  grDevices utils     datasets  methods   base 