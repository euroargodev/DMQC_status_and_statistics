# DMQC_status_and_statistics
    A. DMQC related scripts
1.	DMQC statistics (get_DMQC_stats.m)<br />
Description of the script:<br />
Gets DMQC statistics for given floats.<br />
Notes:<br />
(1)	Reading index file may take 2 or 3 minutes<br />
(2)	It is better to use only ascent profiles (Dprof = 0) for coherence with figures from “get_DMQC_adjustment.m” script<br />
    
   
 
2.	Salinity adjustments for a list of floats (get_DMQC_adjustment.m)<br />
Description of the script:<br />
Shows psal adjustment for a list of floats<br />
Notes: <br />
(1) Getting data process may take some minutes (depending on number of floats) because we are checking all profile files<br />
(2) Descent profiles are not included<br />
(3) Format options (number of floats per figure and yaxis ticks size) can be modified<br />
 
 
      B. TOOLBOX
Name & Description of the auxiliary functions:<br />
- map_tech_param: plot a scatter of a tech, traj or config parameter in a map (lon, lat) using a threshold	read_csv
- get_floats_data_gdac
- M_MAP: matlab package
- export_fig: package
- get_csv_QGIS
- get_phases_times: calculates duration and speed of floats in each cycle phase (descent to park, parking, descent to profile, profile drift and ascent to surface)	suptitle
- get_floats_data_gdac
- get_floats_filespath
- export_fig
- get_last_config_values: get list of last configuration parameters values	
- read_csv
- get_floats_data_gdac<br />

Important functions:
- read_csv:	read a csv file and generates an struct with file variables	
- get_floats_files_paths: gets files paths from ar_index_global_meta.txt file for given floats and optionally creates a .txt file	
- get_data_from_index:gets chosen variables from index file for given floats	
- get_floats_data_gdac: gets data from tech, traj, aux, meta or index files	get_traj_param
- get_tech_param
- get_configparam_meta
- get_data_from_index
- get_floats_data_gadc_format: gets data from tech, traj, aux or meta file and formats it to have the same number of cycles for the same float	get_traj_param
- get_tech_param
- get_configparam_meta
- get_data_from_index
- get_configparam_meta:	gets all cycles values of a list of given configuration parameters from a meta file (one float)	
- get_tech_param:	gets all cycles values of a list of given technical parameters from a tech file (one float)	
- get_traj_param:	gets all values of a list of given parameters from a trajectory file	
- get_traj_param_AUX		
