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
  - grep (Matlab grep equivalent function)

 **WARNING 1** : Profile_QC for PRES information is not yet available in the Argo detailed index. <br />
 It is filled with qc="X" in the scriptfor the moment, and plots related to pres profile<br />
 qc are skipped.<br />

 **WARINING 2** : the detailed index has an issue with some psal adjustments 
 (not filled when it should). A correction was asked to the dev team. 
 In the meanwhile, plots with PSAL adjustment may not be complete.

 **Author**: Euro-Argo ERIC (contact@euro-argo.eu)<br />

 **Version**: 3.2 (2023/07/21)<br />

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
 - V3.2 (2023/07/21) :
   - correct bug at the final zipping step when file does not exist<br />
     (in case i_bgc = 0).
   - add case All Argo Fleet and manage optimization section for large number of floats.
   - dealing with float_wmo type for cases when not all WMOs are coded on 7 characters.
   - dealing with graph layout when n_countries is large
   - special workaround for float 4900566 that used QC 1 instead of QC A for profile QC.
   - add quotes in synthese output for program, in case comma is used.


## B. Graphical outputs for **get_DMQC_stats.m** 
Different outputs are produced: graphical and textual. Here after are examples of graphical outputs obtained for floats from the MOCCA project.

- __Plot 01 and 02__: State of the DMQC per country, in number of floats and/or profiles (one plot per parameter) <br />

These plots present the number of floats (resp. profiles) by country:
total number, number at least DMQCed once, number older than {sage} year  and number at least DMQCed once and older than {sage} year. By default, sage = 1 year.
Here country refers to country associated to the float in the OceanOPS database. It is coded using the OTAN trigramme code. The used correspondence can be found in the file country_code_{yyy-mm-dd}.txt. N.B.: For the European case and the European BGC case, the countries outside Europe but for which floats are decoded at Coriolis where gathered under the 'EXT' trigramme.

<p float="left">
<img 
src="https://github.com/delphinedobler/DMQC_status_and_statistics/blob/main/OUTPUT_examples/MOCCA_case/Plots/01_MOCCA_Fleet_CTD_DMQC_status_nb_floats_by_country_20230713.png" width="400" /> 
<img src="https://github.com/delphinedobler/DMQC_status_and_statistics/blob/main/OUTPUT_examples/MOCCA_case/Plots/02_MOCCA_Fleet_CTD_DMQC_status_nb_profiles_by_country_20230713.png" width="400" />
</p>

- __Plot 03__: Profile quality (all profiles and only Delayed Mode - i.e. consolidated - profiles) in number of floats profiles (one plot per parameter) <br />

This plot presents the number of profiles with respect to the profile_QC code, both for all processed mode and for D-mode only. 
Profile QC codes are defined in the Reference table 2a of the Argo QC manual (https://archimer.ifremer.fr/doc/00228/33951/32470.pdf) and is recalled hereafter:
GOOD data = QC flag values of 1, 2, 5 or 8
BAD data = QC flag values of 3 or 4
- profile QC A: 100% GOOD data (All the measurement points of the profile are GOOD data);
- profile QC B: 75% to 100% GOOD data;
- profile QC C: 50% to 75% GOOD data;
-	profile QC D: 25% to 50% GOOD data;
-	profile QC E:  0% to 25% GOOD data;
-	profile QC F: no good data GOOD data

The exact number and relative percentages are also indicated above the bars. 
The relative percentages for the profiles processed in delayed mode provide a consolidated view.
A few profiles do not have a profile QC in the index file. This observation deserves further analysis.


<p float="center">
<img 
src="https://github.com/delphinedobler/DMQC_status_and_statistics/blob/main/OUTPUT_examples/MOCCA_case/Plots/03_MOCCA_Fleet_PSAL_profile_QC_20230713.png" width="400" /> 
 <img 
src="https://github.com/delphinedobler/DMQC_status_and_statistics/blob/main/OUTPUT_examples/MOCCA_case/Plots/03_MOCCA_Fleet_TEMP_profile_QC_20230713.png" width="400" /> 
</p>

- __Plot 04 and 05__ Profile quality evolution (all profiles and only Delayed Mode - i.e. consolidated - profiles) in percentage of floats profiles (one plot per parameter) <br />

The upper panel of this plot is a time evolution view of plot 03 for all profiles (plot 04) and D-profiles only (plot 05), with respect to the float launch year (rapid proxi for sensor generation). To indicate the significance of the statistics, the number of profiles for the corresponding year is also provided on the lower panel.

<p float="center">
<img 
src="https://github.com/delphinedobler/DMQC_status_and_statistics/blob/main/OUTPUT_examples/MOCCA_case/Plots/04_MOCCA_Fleet_PSAL_profile_QC_evolution_20230713.png" width="400" /> 
 <img 
src="https://github.com/delphinedobler/DMQC_status_and_statistics/blob/main/OUTPUT_examples/MOCCA_case/Plots/05_MOCCA_Fleet_PSAL_Dprofile_QC_evolution_20230713.png" width="400" /> 
</p>

- __Plot 06__ Global DMQC status and grey list information (one plot per parameter) <br />

The first four coloured bars provide the same information as Plot 01 but summed for all countries. 
The three bars on the right side provide information about grey listing (QC3 or QC4):
   - number for active floats (and in the legend, the relative part of all active floats is mentioned)
   - number of inactive floats that still have profiles not processed in delayed mode
   - number of inactive floats that have been fully processed in delayed mode.
Here active floats refers to floats having emitted a profiles within the last 30 days (with respect to the index file update date). 
This limit is arbitrary and can be tuned.


<p float="center">
<img 
src="https://github.com/delphinedobler/DMQC_status_and_statistics/blob/main/OUTPUT_examples/MOCCA_case/Plots/06_MOCCA_Fleet_PSAL_DMQC_status_and_grey_list_20230713.png" width="400" /> 
</p>

- __Plot 07 and 08__ DMQC status per profile year, and age histogram of non-DMQC profiles (one plot per parameter) <br />

Plot 07 presents the time evolution of the percentage of profiles processed in delayed mode with respect to the profile date.
Plot 08 presents the age histogram of profiles with no DMQC performed yet.

<p float="center">
<img 
src="https://github.com/delphinedobler/DMQC_status_and_statistics/blob/main/OUTPUT_examples/MOCCA_case/Plots/07_MOCCA_Fleet_CTD_prof_DMQCstatus_byyear_20230713.png" width="400" /> 
 <img 
src="https://github.com/delphinedobler/DMQC_status_and_statistics/blob/main/OUTPUT_examples/MOCCA_case/Plots/08_MOCCA_Fleet_CTD_prof_DMQCstatus_agehist_20230713.png" width="400" /> 
</p>

- __Plot 09__ R-A-D status for all parameters <br />
Plot 09 presents by parameter (x axis), the number of R-profiles, A-profiles and D-profiles.

<p float="center">
<img 
src="https://github.com/delphinedobler/DMQC_status_and_statistics/blob/main/OUTPUT_examples/MOCCA_case/Plots/09_MOCCA_Fleet_prof_RAD_mode_per_param_20230713.png" width="400" /> 
</p>

- __Plot 10 and 11__ DMQC and quality profile status by batch of WMOs per cycle (one plot per parameter) <br />
These plots show the DMQC (plot 10) and quality profile status (plot 11) by batch of WMOs per cycle (one plot per parameter).
These plots are output only on demand. The number of WMOs shown by graph can be tuned.

<p float="center">
 <img 
src="https://github.com/delphinedobler/DMQC_status_and_statistics/blob/main/OUTPUT_examples/MOCCA_case/Plots/10_MOCCA_Fleet_CTD_RAD_mode_per_wmo_per_cycle_001_20230713.png" width="400" /> 
 <img 
src="https://github.com/delphinedobler/DMQC_status_and_statistics/blob/main/OUTPUT_examples/MOCCA_case/Plots/11_MOCCA_Fleet_PSAL_profile_QC_per_wmo_per_cycle_001_20230713.png" width="400" /> 
</p>

- __Plot 12__ PSAL adjustment by batch of WMOs per cycle  <br />
This plot shows PSAL_adjustment by batch of WMOs per cycle.
This plot is output only on demand. The number of WMOs shown by graph can be tuned.
In grey color: the profiles that are not yet processed in delayed mode and that are not profile QC F.
In black: the profiles (either real-time or delayed mode) that are QC-F.
In jet colorscale, the value of the PSAL adjustment bounded by [-0.07 0.07]. The same bounds are used for all plots for better intercomparison and to ensure that "no adjustment" case will always appear in green.

/!\ WARNING 2: There is an issue with Argo detailed index: for a few floats, PSAL_adjustment is not computed (https://gitlab.ifremer.fr/coriolis/actions/actions-argo/-/issues/63).

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

-  DMQC_status_per_wmo_for_{param}_{yyyymmdd}.txt is a synthesis by wmo. This output can be used to define priorities and the list of floats to be treated in DMQC with respect to some given-criteria. There is quite a number of indication. Some may be added if needed. The associated parameter is recalled in the first column to make sure that the information is understood as being relative to this parameter only. The column name are quite straightforward, however a few may deserve some more details:
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
- grep : Matlab equivalent of the unix grep command (not as performant).
