--------------------------------------------------------------------------------------------------------

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
3.  Scripts and Data (including Analyses_ready_data)

--------------------------------------------------------------------------------------------------------

1. Supplementary Actograms contain distances to the roost for each studied individual:

--------------------------------------------------------------------------------------------------------

2. Supplementary Methods describe data collections, preparation, analyses and visualisation of the red knot data.

--------------------------------------------------------------------------------------------------------

3a. Scripts used to generate the output presented in the paper
	- knot_analyses generates the resutls, figures and supporting actograms for the red knot case study
	- Analyses_shorebird_incubation_tide.R contains analyses for shorebird incubation case study 
	
--------------------------------------------------------------------------------------------------------

3b. Data contain
	a) "median positions.txt" - median half hourly positions to the closest roost for each bird
	   - column definitions:
 		1. Tag		: unique ID of the bird
 		2. timestamp	: 
 		3. Timestamp.loc: 
 		4. date		:
 		5. day		:
 		6. rhour	: 
 		7. x		:
 		8. y		:
 		9. ClosestRoost	:

	b) low tides.txt - contains times of low tides at the study site obtained from 
	  	1. Date	
		2. DayNr	
 		3. Time	
		4. Height		
 		5. HighLow	
 		6. Timestamp.loc	
 		7. Timestamp.loc2	
 		8. hour	
		9. min	
		10. day	
		11. Number
   
	c) high tides.txt - contains times of high tides at the study site obtained from 
	  	1. Date	
		2. DayNr	
 		3. Time	
		4. Height		
 		5. HighLow	
 		6. Timestamp.loc	
 		7. Timestamp.loc2	
 		8. hour	
		9. min	
		10. day	
		11. Number

	d) subShab.tif

	e) data.Rdata contains ready to be analysed (i.e. can be used without need to generate); object 'd' contains the off-roost distances and relevant variables, object md contains tag_ID and 95% farthest distance used for standardization of all distance for the given bird
	
	f) simulations folder contains simulation runs used to generate resutls in the paper 
	   - model_outputs.R contains object 'x' with utput of single models for each nest
		1. tag		: unique number of the bird 
		2. n		: number of half hourly observations
		3. est		: model estimate (median of posterior distribution of 10 000 simulated values)
		4. lwr		: 0.025 percentile posterior distribution of 10 000 simulated values
		5. upr		: 0.975 percentile of posterior distribution of 10 000 simulated values 
		6. period	: period 24h or 12.4h to wich the estimate refers
		7. vari		: sin or cos
	
	   - fits_day_htide_nolog.txt contains predicted values with 95% credible intervals for relationship of distance with tidal rhythm
		1. rhour	: time of day in hours 
		2. rad_t	: mean tidal time in radians
		3. rad_dn	: time of day in radians
		4. pred		: predictions (that is standardized sistance to the closest roost)
		5. lwr		: 0.025 percentile posterior distribution of 10 000 predicted values
		6. upr		: 0.975 percentile of posterior distribution of 10 000 predicted values 
		7. tag		: unique number of the bird 

	   - fits_day_htide_nolog.txt contains predicted values with 95% credible intervals for relationship of distance with tidal rhythm
		1. nTideTime	: time to (negative) and from (positive) high tide 
		2. rad_t	: tidal time in radians
		3. rad_dn	: mean day time in radians
		4. pred		: predictions (that is standardized sistance to the closest roost)
		5. lwr		: 0.025 percentile posterior distribution of 10 000 predicted values
		6. upr		: 0.975 percentile of posterior distribution of 10 000 predicted values 
		7. tag		: unique number of the bird 


--------------------------------------------------------------------------------------------------------
	 
	
