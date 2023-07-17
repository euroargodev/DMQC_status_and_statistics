# DMQC_status_and_statistics

## A. Script **get_DMQC_stats.m** <br />

This script computes DMQC statistics for a given list of floats <br />
<br />
**INPUTS**
 - **argo_profile_detailled_index.txt** and <br />
   **argo_synthetic-profile_detailled_index.txt** if i_BGC=1
 - **ar_greylist.txt**
 - **floats list**: csv file, separated by ";", with 4 fields:<br />
    WMO; COUNTRY; LAUNCH_DATE; PROGRAM<br />
        COUNTRY is the country in charge of the delayed mode processing (from OceanOPS)<br />
        LAUNCH_DATE must be in the format "YYYY/MM/DD".<br />
   e.g.:<br />
   WMO;COUNTRY;LAUNCH_DATE;PROGRAM<br />
   3901496;United Kingdom;2014/10/20;Argo UK Bio<br />
   3901497;United Kingdom;2014/10/24;Argo UK Bio<br />
   These 4 fields can be extracted from the OceanOPS website or directly downloaded from <br />
   https://www.ocean-ops.org/share/Argo/Status/<br />
   In the argo_all.csv from OceanOPS, the corresponding columns to account for are :<br />
         REF (-> WMO), COUNTRY , DEPL_DATE (-> LAUNCH_DATE) and PROGRAM 
 - **country_code.csv**: csv file, separated by ";", with 2 fields:<br />
    COUNTRY; COUNTRY_CODE<br />
       COUNTRY should follow the OceanOPS conventions<br />
       COUNTRY_CODE is a 3-digit CODE that will be used on graphs and outputs.<br />
    The default template country_code_template.csv from the script directory can be used.<br />
    This default contains all the countries associated to Argo in OceanOPS with the<br />
    corresponding OTAN (https://en.wikipedia.org/wiki/List_of_NATO_country_codes) official <br />
    3-digits codes (list of relevant countries extracted the 2023/07/07 from OceanOPS database)<br />
<br />

 **CONFIGURATION PARAMETERS**
 - **i_descending_profile**: 1 means descending profiles are considered, <br />
                         0 means descending profiles are not considered.
 - **sage**: threshold (in days) for floats and observations statistics.
 - **i_bgc**: 1 means bgc profiles/parameters are considered (the detailed argo synthetic index will be <br />
            read for BGC parameters, the detailed argo index will be read for CTD)<br />
          0 means core information (CTD) is analysed from the detailed argo index.
 - **input_list_of_parameters_to_treat** is the list of parameters to analyse
 - **print_svg**: 1 means figures will be saved in .svg format as well<br />
 (interesting for high quality, but a little longer to save).
 - **output_graphs_per_float**: flag to indicate if graphs per float should be<br />
   recorded. (Graphs will be recorded by bunch of 40 floats max. For treatment with a large<br />
   number of floats, this may not be relevant)
 - **n_max_float_per_graph**: associated to output_graphs_per_float.
<br />

**OUTPUTS**
 - **Figures**   saved in folder outputs_yyyy-mm-dd/Plots 
 - **Analyses**  saved in folder outputs_yyyy-mm-dd/Syntheses 
 - **Copy of input files** saved in folder outputs_yyyy-mm-dd

 **Auxiliary functions needed**
  - read_csv
  - get_data_from_index 
  - plotBarStackGroups 

 **WARNING** : Profile_QC for PRES information is not yet available in the Argo detailed index. <br />
 It is filled with qc="X" in the scriptfor the moment, and plots related to pres profile<br />
 qc are skipped.<br />

 **Author**: Euro-Argo ERIC (contact@euro-argo.eu)<br />

 **Version**: 3.1 (2023/07/12)<br />

 **Historic**:<br />
 - V1.0 : This script originally created by Andrea Garcia Juan and Romain<br />
        Cancouët, and updated by Luca Arduini Plaisant.
 - V2.0 (2023/06/19): 
   - The script architecture was reviewed on 2023/06/19 by Delphine Dobler <br />
     to include the processing of BGC floats and to enhance performances.
 - V2.01 (2023/07/07): 
   - adding one test for the existence of found indices
   - only keep BGC parameters that were found in the synthetic index file
   - output an additional file with this information
 - V2.1  (2023/07/07): 
   - fetch DAC information from the index
   - record D-profile last update date in output
   - associate country_codes from a config file
   - add -f option to final zip to allow overwrite if script is<br />
        launched several times on the same day
 - V3.0  (2023/07/10):
   - merge CTD and BGC processings (no need for the user to know<br />
          which WMO is BGC, which one is not)
 - V3.1  (2023/07/12):
   - add plot with the number of R/A/D profiles per variable
   - add plots with information by cycle and by WMO:
     - R/A/D status
     - profile QC status
     - PSAL_adjustement<br />
     => These plots replace the old get_DMQC_adjustment.m script
   - stop duplicating graphs for CTD mode
   - skip plots with PRES profile_QC (information not yet available in index)


## B. Graphical outputs for **get_DMQC_stats.m** 
Different outputs are produced: graphical and textual. Here after are examples of graphical outputs obtained for floats from the MOCCA project.

- __Plot 01 and 02__: State of the DMQC per country, in number of floats and/or profiles (one plot per parameter) <br />

<p float="left">
<img 
src="https://github.com/delphinedobler/DMQC_status_and_statistics/blob/main/OUTPUT_examples/MOCCA_case/Plots/01_MOCCA_Fleet_CTD_DMQC_status_nb_floats_by_country_20230713.png" width="400" /> 
<img src="https://github.com/delphinedobler/DMQC_status_and_statistics/blob/main/OUTPUT_examples/MOCCA_case/Plots/02_MOCCA_Fleet_CTD_DMQC_status_nb_profiles_by_country_20230713.png" width="400" />
</p>

- __Plot 03__: Profile quality (all profiles and only Delayed Mode - i.e. consolidated - profiles) in number of floats profiles (one plot per parameter) <br />
<p float="center">
<img 
src="https://github.com/delphinedobler/DMQC_status_and_statistics/blob/main/OUTPUT_examples/MOCCA_case/Plots/03_MOCCA_Fleet_PSAL_profile_QC_20230713.png" width="400" /> 
 <img 
src="https://github.com/delphinedobler/DMQC_status_and_statistics/blob/main/OUTPUT_examples/MOCCA_case/Plots/03_MOCCA_Fleet_TEMP_profile_QC_20230713.png" width="400" /> 
</p>

- __Plot 04 and 05__ Profile quality evolution (all profiles and only Delayed Mode - i.e. consolidated - profiles) in percentage of floats profiles (one plot per parameter) <br />
<p float="center">
<img 
src="https://github.com/delphinedobler/DMQC_status_and_statistics/blob/main/OUTPUT_examples/MOCCA_case/Plots/04_MOCCA_Fleet_PSAL_profile_QC_evolution_20230713.png" width="400" /> 
 <img 
src="https://github.com/delphinedobler/DMQC_status_and_statistics/blob/main/OUTPUT_examples/MOCCA_case/Plots/05_MOCCA_Fleet_PSAL_Dprofile_QC_evolution_20230713.png" width="400" /> 
</p>

- __Plot 06__ Global DMQC status and grey list information (one plot per parameter) <br />
<p float="center">
<img 
src="https://github.com/delphinedobler/DMQC_status_and_statistics/blob/main/OUTPUT_examples/MOCCA_case/Plots/06_MOCCA_Fleet_PSAL_DMQC_status_and_grey_list_20230713.png" width="400" /> 
</p>

- __Plot 07 and 08__ DMQC status per profile year, and age histogram of non-DMQC profiles (one plot per parameter) <br />
<p float="center">
<img 
src="https://github.com/delphinedobler/DMQC_status_and_statistics/blob/main/OUTPUT_examples/MOCCA_case/Plots/07_MOCCA_Fleet_CTD_prof_DMQCstatus_byyear_20230713.png" width="400" /> 
 <img 
src="https://github.com/delphinedobler/DMQC_status_and_statistics/blob/main/OUTPUT_examples/MOCCA_case/Plots/08_MOCCA_Fleet_CTD_prof_DMQCstatus_agehist_20230713.png" width="400" /> 
</p>

- __Plot 09__ R-A-D status for all parameters <br />
<p float="center">
<img 
src="https://github.com/delphinedobler/DMQC_status_and_statistics/blob/main/OUTPUT_examples/MOCCA_case/Plots/09_MOCCA_Fleet_prof_RAD_mode_per_param_20230713.png" width="400" /> 
</p>

- __Plot 10 and 11__ DMQC and quality profile status by batch of WMOs per cycle (one plot per parameter) <br />
<p float="center">
 <img 
src="https://github.com/delphinedobler/DMQC_status_and_statistics/blob/main/OUTPUT_examples/MOCCA_case/Plots/10_MOCCA_Fleet_CTD_RAD_mode_per_wmo_per_cycle_001_20230713.png" width="400" /> 
 <img 
src="https://github.com/delphinedobler/DMQC_status_and_statistics/blob/main/OUTPUT_examples/MOCCA_case/Plots/11_MOCCA_Fleet_PSAL_profile_QC_per_wmo_per_cycle_001_20230713.png" width="400" /> 
</p>

- __Plot 12__ PSAL adjustment by batch of WMOs per cycle  <br />
 <img 
src="https://github.com/delphinedobler/DMQC_status_and_statistics/blob/main/OUTPUT_examples/MOCCA_case/Plots/12_MOCCA_Fleet_PSAL_PSAL_adj_per_wmo_per_cycle_001_20230713.png" width="400" /> 
</p>

## C. Syntheses outputs for **get_DMQC_stats.m** 

There are 4 kinds of syntheses produced by the script:

-  Additional_Info_{yyyymmdd}.txt
 This file indicates which WMO, if any, were not found in the Argo detailed index, and in the Argo detailed synthetic index if i_bgc = 1.
If i_bgc=1, it also indicates which BGC parameters, if any, were not found in the Argo detailed synthetic index for the input list of WMOs. On the other hand, it indicates which BGC parameters were found and were not requested.

-  DMQC_status_per_country_for_{param}_{yyyymmdd}.txt
 This file is the numbered information for plots 01 and plots 02.

-  DMQC_status_per_wmo_for_{param}_{yyyymmdd}.txt is a synthesis by wmo. This output can be used to define priorities and the list of floats to be treated in DMQC with respect to personnal criteria. There is quite a number of indication. Some may be added if needed. The associated parameter is recalled in the first column to make sure that the information is understood as being relative to this parameter only. The column name are quite straightforward, however a few may deserve some more details:
   - The nb_prof_QC_X and nb_prof_DM_QC_X refer to profiles with no value for profile QC (see Plot 03 comment above).
   - DM_done column refers to the fact that this WMO has been seen in delayed mode at least once. To get the delayed mode completeness, refer to percentage_DM_prof column.
   - greylist means that the float was put on greylist with QC3 or QC4 for the given parameter.

- The DMQC_status_warnings_per_wmo_for_{param}_{yyyymmdd}.txt is an output of what was thought should raise the attention such as:
  - there is at least one profile_QC set to F (prof_QC=F_once)
  - there is at least one profile_QC with no value (prof_QC_not_filled_once)
  - the DMQC was never performed whereas the float is older than 1 year (No_DMQC_older_1yr) (Nota bene: this may be a criterion little bit too restrictive ...)
  - The float is in grey list.   

## D. TOOLBOX
Name & Description of the auxiliary functions:
- read_csv:	read a csv file and generates an struct with file variables	
- get_data_from_index: gets chosen variables from index file for a given list of floats
- plotBarStackGroups: Permit to make a bar plot with stacked bars for one graph tick.
