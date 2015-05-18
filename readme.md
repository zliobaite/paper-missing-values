Data and code in this folder is to support the article:

**Zliobaite, I., Hollmen, J. and Junninen, H. (2014). Regression models tolerant to massively missing data: a case study in solar radiation nowcasting. Atmospheric Measurement Techniques 7, 4387-4399.**

The code is in MATLAB. 

The data and the code can be used for research purposes, provided that the above article is cited.

Mailto: zliobaite at gmail.com

Last updated: 2014 06 18


## Contents ##

DATA

	data_smearii.csv
	data_theoretical_radiation.csv
	
MAIN CODE (EXPERIMENTS)

	run_statistics.m
	run_case_study.m
	run_case_study_sensitivity.m
	
SUPPORTING CODE (FUNCTIONS)

	error_reg.m
	nipals_train_batch_nomean.m
	reg_regression_train.m
	remove_missing_values.m
	standardize_data_nan_train.m
	standardize_data_nan.m
	standardize_back.m
	pca_reg.m


## Reproducing the experiments ##

	Figure 2 and Figure 4
	run_statistics.m

	Table 3 and Table 4
	run_case_study.m

	Figure 5
	run_case_study_sensitivity.m

## Data ##

data_smearii.csv 
contains the dataset from SMEAR II station (http://www.atm.helsinki.fi/SMEAR/index.php/smear-ii). 

Number of observations: 140 576
Number of variables: 6 (date) + 37 (sensor readings)

Variables:

1. year
2. month
3. day
4. hour
5. minute
6. second

7. Rain 18.0m
8. SWS 18.0m 
9. Dew point 18.0m 
10. P 0.0m 
11. T 4.2m 
12. T 8.4m 
13. T 16.8m 
14. T 33.6m 
15. T 50.4m 
16. T 67.2m 
17. WS 33.6m 
19. WS 8.4m 
20. WS 16.8m 
21. WS 33.6m 
22. WS 74.0m 
23. WD avr 
24. WD ultrasonic 8.4m 
25. WD ultrasonic 16.8m 
26. WD ultrasonic 33.6m 
27. WD ultrasonic 74.0m 
28. RH 4.2m 
29. RH 8.4m 
30. RH 16.8m 
31. RH 33.6m 
32. RH 50.4m 
33. RH 67.2m 
34. RH Td 18.0m 
35. PTG 
36. Visibility 18.0m 
37. Vis-min 18.0m 
38. Vis-max 18.0m 
39. Precipitation intensity 18.0m 
40. Preci-min 18.0m 
41. Preci-max 18.0m 
42. Precipitation 18.0m 
43. Snowfall 18.0m 
44. Global radiation 18.0m 


data_theoretical_radiation.csv 	
contains the dataset using MIDC SOLPOS
Calculator (http://www.nrel.gov/midc/solpos/solpos.html). 

Parameters used: Surface pressure 990 mbar, Ambient dry-bulb
temperature 3oC$, Azimuth of panel surface 180o$, Degrees tilt
from horizontal of panel 0, Solar constant 1367 W/m2, Shadow-band
width 7.6 cm, Shadow-band radius 31.7 cm, Shadow-band sky factor
0.04, Interval of a measurement period 0 sec.


Number of observations: 140 576
Nuber of variables: 6 (date) + 1 (theoretical radiation)

Variables:

1. year
2. month
3. day
4. hour
5. minute
6. second

7. theoretical radiation

## Code ##

* run_statistics.m - reproducing results in Figure 2 and Figure 4
* run_case_study.m - reproducing results in Table 3 and Table 4
* run_case_study_sensitivity.m - reproducing results in Figure 5


* error_reg.m - computes prediction error
* nipals_train_batch_nomean.m - implements PLS regression (batch training)
* reg_regression_train.m - implements Ridge Regression (training)
* remove_missing_values.m - removes instances with missing labels and replaces missing 		values in the input data by 0
* standardize_data_nan_train.m - transforms data to 0 mean and 1 variance and saves the 		original mean and variance
* standardize_data_nan.m - transforms data subtracting a given mean and dividing by a given 	variance 
* standardize_back.m - reverts standardized data back to the original mean and variance
* pca_reg.m - finds PCA rotation matrix