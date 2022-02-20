
              EX2_Project: Louis LEMAIR et Cathy GELLENONCOURT

In this file, we provide the codes related to the EX2 project on the STELLA experiment. You will find:
- a file named Data;
- libPerso.h;
- Calibration_Eu.C;
- calib_each_Time.C;
- fit_selfRadio.C and fit_selfRadio.h;
- GetAllMeans.C and GetAllMeans_21.C;
- hist_superposition.c;
- meanEvolvTime.C and meanEvolvTime_21.C;
- energy_drift.C;
- selectTime.C.


In the Data file, you will find the various data needed by our code and the data provided by our code. 
The data needed is in the following root file : 
- eu152_20220203.root => data for calibration with an europium source made the 03.02.2022;
- eu152_20220204.root => data for calibration with an europium source made the 04.02.2022;
- fixData.root => data without source;
- la138_array_fullShell.root => data related to the simulation;

The data we provide is :  
- Eu152_calibrated_20220203.root => data calibrated with the data stored the 03.02.2022 with the file eu152_20220203.root;
- Eu152_calibrated_20220204.root => data calibrated with the data stored the 04.02.2022 with the file eu152_20220204.root;
- Eu152_wrong_cal.root => data calibrated with the coefficient of the eu152_20220204.root calibration data file but applied on the data stored the 03.04.2022 with the file eu152_20220203.root;
- Simulation_factors.root => results (ratio and offset) of the simulation;
- Simulation_results.root => data extracted from la138_array_fullShell.root;
- Superposition_calibrated.root => storage of the calibrated histogram superposition.

The codes provided are : 
- libPerso.h => contains libraries needed by the code;
- Calibration_Eu.C => make the calibration with the europium source. You have to specify as an option if you want a calibration with data of the 03.02.2022 (option 3) or with data of the 04.02.2022 (option 4);
- calib_each_Time.C => allows to get the parameters of the calibration over the 21h of data taking. This is an interpolation of the ratio and offset coefficients (represents the way we go from QDC to Energy);
- GetAllMeans.C and GetAllMeans_21.C => this code allows to get the data file Simulation_factors.root. We get all ratio and offset coefficient form the simulation. GetAllMeans.C gives us the coefficients from the third to the 16th hour while GetAllMeans_21.C gives us the coefficients over the 21h of data taking;
- hist_superposition.c => makes the superposition of the cumulated histogram from : calibration_20220203, interpolation and simulation;
- meanEvolvTime.C and meanEvolvTime_21.C =>plot the evolution of the mean in function of time. meanEvolvTime.C gives us the evolution from the third to the 16th hour while meanEvolvTime_21.C gives us the evolution over the 21h of data taking. meanEvolvTime.C has to be run AFTER GetAllMeans.C and meanEvolvTime_21.C has to be run AFTER GetAllMeans_21.C. The user can specify different options : 
	- 0 each peek is fitted but nothing is calibrated;
	- 1 peeks are fitted and only one calibration (20220203) is used;
    	- 2 peeks are fitted and calibration is done via Intrepolation of 2-times calibration;
    	- 3 peeks are fitted on simulation calibration;
    	- 4 : 1 and 2 together;
    	- 5 : 2 and 3 together.
- fit_selfRadio.C and fit_selfRadio.h => the original code provided to fit the data to the simulation;
- energy_drift.C => the original code provided to display the energy drift of the background from 2022-02-03;
- selectTime.C => the original code provided to select a background spectrum by specifying the time window.

