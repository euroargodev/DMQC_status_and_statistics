# DMQC_status_and_statistics
       A. Scripts
1.	DMQC statistics (get_DMQC_stats.m)<br />
Description of the script:<br />
Gets DMQC statistics for a given float list.<br />
Notes:<br />
(1)	Reading index file may take 2 or 3 minutes<br />
(2)	It is better to use only ascent profiles (Dprof = 0) for coherence with figures from “get_DMQC_adjustment.m” script<br />
(3) The input list of floats must be divided in 3 columns: WMO, RT (dac in charge of real time
processing) and DM (institution in charge of delayed mode processing), separated by ';'.<br />
(4) For now, some errors might occur if the user input a list of floats that are not all on the GDAC. This problem should be fixed soon.<br />
(5) Also, the supposed "txt" output is for now commented because it was source for error. Will be fixed in a future update.<br />
(6) The bio-argo-index-detailled is, since mid 2021, contaning the quality flags of the different BGC variables, improving the outputs of this script.<br />

Different Outputs produced:
- State of the DMQC per DACs, in number of floats and/or observations<br />
![alt text](https://github.com/euroargodev/DMQC_status_and_statistics/blob/main/Images/DMQC_stats_FSD_float_DMQCstatus_20210414.png?raw=true)
- State of the DMQC per years, in number of observations
![alt text](https://github.com/euroargodev/DMQC_status_and_statistics/blob/main/Images/DMQC_stats_FSD_obs_DMQCstatus_byyear_20210414.png?raw=true)
- State of the DMQC per variables, in number of observations
![alt text](https://github.com/euroargodev/DMQC_status_and_statistics/blob/main/Images/DMQC_stats_DOXY_params_data_mode_flags_20210407.png?raw=true)
- Number of floats greylisted in the input list and for which variable
![alt text](https://github.com/euroargodev/DMQC_status_and_statistics/blob/main/Images/DMQC_stats_FSD_greylisted_floats_per_variables_20210414.png?raw=true)

   
2.	Variable adjustments for a list of floats (DMQC_adjustment.m)<br />
Description of the script:<br />
Shows variable adjustment and state of DMQC for a list of floats. Contrary to "DMQC_stats.m", this script gathers data from the profiles files and not an index.<br />
Notes: <br />
(1) Getting data process may take some minutes (depending on number of floats) because the script is checking all profile files<br />
(2) Descent profiles are not included<br />
(3) Format options (number of floats per figure and yaxis ticks size) can be modified<br />

Different Outputs produced:
- State of the DMQC per floats numbers in terms of number of cycle<br />
![alt text](https://github.com/euroargodev/DMQC_status_and_statistics/blob/main/Images/DMQC_adjustments_FSD_CTD_DMQCstatusScatter_lot2_20210412.png?raw=true)
- Most frequent QC flags for the variable defined as input<br />
![alt text](https://github.com/euroargodev/DMQC_status_and_statistics/blob/main/Images/DMQC_adjustments_FSD_CTD_mostfreqPSALQC_lot2_20210412.png?raw=true)
- Variable QC flags per cycles<br />
![alt text](https://github.com/euroargodev/DMQC_status_and_statistics/blob/main/Images/DMQC_adjustments_FSD_CTD_PSALQCscatter_lot2_20210412.png?raw=true)


      B. TOOLBOX
Name & Description of the auxiliary functions:<br />
- read_csv:	read a csv file and generates an struct with file variables	
- get_data_from_index:gets chosen variables from index file for given floats
- export_fig: Export figures as png or jpg
- Flex_legends: Permits to create a more customizable legend in a plot
- Split_figures : Function splitting the figures in different outputs if the number of floats exceed a defined value.
- plotBarStackGroups: Permit to make a bar plot with stacked bars for one graph tick.
