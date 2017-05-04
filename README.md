# Living_with_the_tides
				
Description of the Supporting information, including scripts and data to generate the results and visualisations, from 				
  'Marine biorhythms: bridging chronobiology and ecology' 				
   by Bulla et al 2017 (Phyl Trans B)				
				
--------------------------------------------------------------------------------------------------------				
				
WHEN USING any of the Supporting information, PLEASE CITE both the original paper and the Supporting information				
	Bulla, M. at al. (2017).  'Marine biorhythms: bridging chronobiology and ecology'. Phyl Trans. 			
	Bulla, M., Oudman, T., & Bijleveld, A.I. 2017. Supporting information for “Marine biorhythms: bridging chronobiology and ecology” Open Science Framework. http://doi.org/10.17605/OSF.IO/XBY9T Retrieved: ADD DATETIME.			
				
For any publication making substantial use of the data, please contact Martin Bulla (bulla.mar@gmail.com), as the authors welcome the opportunity for collaboration and wish to comment prior to publication.				
				
--------------------------------------------------------------------------------------------------------				
				
CONTENT				
1.  Supplementary Actograms				
2.  Supplementary Methods				
3.  Scripts 
4.  Data
5.  Simulations				
				
--------------------------------------------------------------------------------------------------------				
				
1. Supplementary Actograms - contain distances to the roost for each studied individual:				
				
--------------------------------------------------------------------------------------------------------				
				
2. Supplementary Methods - describe data collection, preparation, analyses and visualisation of the red knot data.				
				
--------------------------------------------------------------------------------------------------------				
				
3. Scripts - contain R-code used to generate the output presented in the paper	
			
	- analyses_knot.R generates the resutls, figures and supporting actograms for the red knot case study			
	- analyses_shorebird_incubation_tide.R contains analyses for shorebird incubation case study 			
				
--------------------------------------------------------------------------------------------------------				
				
4. Data contain				
	   
	   median positions.txt - median half hourly positions to the closest roost for each bird"			
	 	1. Tag			unique ID of the bird	
 		2. timestamp		time in UTC	
 		3. Timestamp.loc	local time	
 		4. date			Date taken from timestamp	
 		5. day			reversed count day (for plotting only)	
 		6. rhour		hour relative to low tide	
 		7. x			latitude in UTM zone 28	
 		8. y			longitude in UTM zone 28	
 		9. ClosestRoost		The closest of the seven defined roost positions	
	   
	   low tides.txt - contains times of low tides at the study site obtained from 			
	  	1. Date			Date	
		2. DayNr		Day number starting at 1 January	
 		3. Time			Time of low tide	
		4. Height		Astronomical water height at low tide	
 		5. Timestamp.loc	Local time of low tide	
 		6. hour			Decimal hours (for plotting)	
		7. min			Minutes after the whole hour	
		8. day			Inverse day number for plotting	
		9. Number		For plotting low tide lines (one number per line)	
	   
	   high tides.txt - contains times of high tides at the study site obtained from 			
	 	1. Date			Date	
		2. DayNr		Day number starting at 1 January	
 		3. Time			Time of high tide	
		4. Height		Astronomical water height at high tide	
 		5. Timestamp.loc	Local time of high tide	
 		7. hour			Decimal hours (for plotting)	
		9. min			Minutes after the whole hour	
		10. day			Inverse day number (for plotting)	
		11. Number		For plotting high tide lines (one number per line)	
	   
	   bda.tif - an image of the Banc d'Arguin, based on band 4 of a Landsat 5 image. Coordinates in UTM zone 28				
	   
	   data.txt - knot-data prepared for analyses (i.e. can be used without need to generate) obtained by 'prepare data for analyses' part in analyses_knot.R script 
		     - contains the off-roost distances and relevant variables) object md contains tag_ID and 95% farthest distance used for standardization of all distance for the given bird"			
	   
	   95%distance.txt - obtained by 'prepare data' part in analyses_knot.R; 
					- contains tag_ID and 95% farthest distance used for standardization of all distance for the given bird			
				
--------------------------------------------------------------------------------------------------------

5. Simulations contain output of simulations runs used to generate resutls and Figure 2 in the paper 			
	   
	   model_output_100.RData - contains object 'x' with summary model estimates for each individual			
		1. tag		: unique number of the bird 
		2. n		: number of half hourly observations
		3. est		: model estimate (median of a joint posterior distribution of 100 model runs with 10 000 simulated values per run )
		4. lwr		: 0.025 percentile of a joint posterior distribution of 100 model runs with 10 000 simulated values per run
		5. upr		: 0.975 percentile of a joint posterior distribution of 100 model runs with 10 000 simulated values per run
		6. period	: period 24h or 12.4h to wich the estimate refers (if vari = int, period = NA)	
		7. vari		: int (intercept), sin or cos
				
	   fits_tide_100s.txt - contains predicted values with 95% credible intervals for relationship of distance with tidal rhythm while keeping the relationship with time of day constant			
		1. rhour	: time of day in hours 	
		2. rad_t	: mean tidal time in radians	
		3. rad_dn	: time of day in radians	
		4. pred		: prediction for each individual (that is standardized distance to the closest roost) based on joint posterior distribution of 10 000 simulated values for each of 100 model runs (each run on a sample of 15% data points for the given individual with two observations at least 1.5h apart)
		5. lwr		: 0.025 percentile of joint posterior distribution of 10 000 predicte values for each of 100 model runs (each run on a sample of 15% data points for the given individual with two observations at least 1.5h apart)
		6. upr		: 0.975 percentile of joint posterior distribution of 10 000 predicte values for each of 100 model runs (each run on a sample of 15% data points for the given individual with two observations at least 1.5h apart)
		7. tag		: unique number of the bird 
				
	   fits_day_100s.txt - contains predicted values with 95% credible intervals for relationship of distance with time of day while keeping the relationship with tide constant			
		1. nTideTime	: time to (negative) and from (positive) high tide 	
		2. rad_t	: tidal time in radians	
		3. rad_dn	: mean day time in radians	
		4. pred		: prediction for each individual (that is standardized distance to the closest roost) based on joint posterior distribution of 10 000 simulated values for each of 100 model runs (each run on a sample of 15% data points for the given individual with two observations at least 1.5h apart)
		5. lwr		: 0.025 percentile of joint posterior distribution of 10 000 predicte values for each of 100 model runs (each run on a sample of 15% data points for the given individual with two observations at least 1.5h apart)
		6. upr		: 0.975 percentile of joint posterior distribution of 10 000 predicte values for each of 100 model runs (each run on a sample of 15% data points for the given individual with two observations at least 1.5h apart)
		7. tag		: unique number of the bird 
				
				
--------------------------------------------------------------------------------------------------------				
	 			

