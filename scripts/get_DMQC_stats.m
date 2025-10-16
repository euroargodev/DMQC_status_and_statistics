% get_DMQC_stats
% This script computes DMQC statistics for a given list of floats
%
%
% INPUTS (paths to update below)
%
% - argo_profile_detailled_index.txt and
%   argo_synthetic-profile_detailled_index.txt if i_BGC=1
%
% - ar_greylist.txt (to replace with the exclusion list/ask CC)
%
% - floats list: csv file, separated by ";", with 4 fields:
%    WMO; COUNTRY; LAUNCH_DATE; PROGRAM
%        COUNTRY is the country in charge of the delayed mode processing (from OceanOPS)
%        LAUNCH_DATE must be in the format "YYYY/MM/DD".
%   e.g.:
%   WMO;COUNTRY;LAUNCH_DATE;PROGRAM
%   3901496;United Kingdom;2014/10/20;Argo UK Bio
%   3901497;United Kingdom;2014/10/24;Argo UK Bio
%   These 4 fields can be extracted from the OceanOPS website or directly downloaded from 
%   https://www.ocean-ops.org/share/Argo/Status/
%   In the argo_all.csv from OceanOPS, the corresponding columns to account for are :
%         REF (-> WMO), COUNTRY , DEPL_DATE (-> LAUNCH_DATE) and PROGRAM 
%
% - country_code.csv: csv file, separated by ";", with 2 fields:
%    COUNTRY; COUNTRY_CODE
%       COUNTRY should follow the OceanOPS conventions
%       COUNTRY_CODE is a 3-digit CODE that will be used on graphs and outputs.
%    The default template country_code_template.csv from the script directory can be used.
%    This default contains all the countries associated to Argo in OceanOPS with the
%    corresponding NATO (https://en.wikipedia.org/wiki/List_of_NATO_country_codes) official 
%    3-digits codes (list of relevant countries extracted the 2023/07/07 from OceanOPS database)
%
%
% CONFIGURATION PARAMETERS:
% - i_descending_profile: 1 means descending profiles are considered, 
%                         0 means descending profiles are not considered.
% - profile_age_method = 'date' or 'days': choose the method to filter "old" profiles
% - profile_age_min_days: "old" threshold (in days) for floats and observations statistics.
% - profile_age_max_date: "old" threshold (in date, format yyyy/mm/dd) for floats and observations statistics.
% - i_bgc: 1 means bgc profiles/parameters are considered (the detailed argo synthetic index will be 
%            read for BGC parameters, the detailed argo index will be read for CTD)
%          0 means core information (CTD) is analysed from the detailed argo index.
% - input_list_of_parameters_to_treat is the list of parameters to analyse
% - print_svg: 1 means figures will be saved in .svg format as well
% (interesting for high quality, but a little longer to save).
% - output_graphs_per_float: flag to indicate if graphs per float should be
%   recorded. (Graphs will be recorded by bunch of 40 floats max. For treatmant with a large
%   number of floats, this may not be relevant)
% - n_max_float_per_graph: associated to output_graphs_per_float. By
% default, set to 30.
%
% Outputs
% - Figures   saved in folder output_files_yyyy-mm-dd_hhmmss/Plots 
% - Analyses  saved in folder output_files_yyyy-mm-dd_hhmmss/Syntheses
% - Copy of input files saved in folder output_files_yyyy-mm-dd_hhmmss.
%
% Auxiliary functions needed:
%    read_csv
%    get_data_from_index
%    plotBarStackGroups
%    grep (Matlab grep equivalent function)
%
% WARNING : Profile_QC for PRES information is not yet available. IndexData
% is filled with qc="X" for the moment, and plots related to pres profile
% qc are skipped.
%
% Author: Euro-Argo ERIC (contact@euro-argo.eu)
%
% Version: 3.5 (2025/10/17)
%
% History:
% V1.0 : This script was originally created by Andrea Garcia Juan and Romain
%        Cancouët, and updated by Luca Arduini Plaisant.
% V2.0 (2023/06/19): 
%        - The script architecture was reviewed on 2023/06/19 by Delphine Dobler 
% to include the processing of BGC floats and to enhance performances.
% V2.01 (2023/07/07): 
%        - adding one test for the existence of found indices
%        - only keep BGC parameters that were found in the synthetic index file
%        - output an additional file with this information
% V2.1  (2023/07/07): 
%        - fetch DAC information from the index.
%        - record D-profile last update date in output.
%        - associate country_codes from a config file
%        - add -f option to final zip to allow overwrite if script is
%        launched several times on the same day.
% V3.0  (2023/07/10):
%        - merge CTD and BGC processings (no need for the user to know
%          which WMO is BGC, which one is not).
% V3.1  (2023/07/13):
%        - add plot with the number of R/A/D profiles per variable
%        - add plots with information by cycle and by WMO:
%               - R/A/D status
%               - profile QC status
%               - PSAL_adjustement
%            => These plots replace the old get_DMQC_adjustment.m script.
%        - stop duplicating graphs for CTD mode, 
%        - skip plots with PRES profile_QC (information not yet available in index)
%        - add work-around treatments when command system are not available
%        - add scalability performance information
%        - correct code warnings 
%        - correct bug for wmos_operational that was overwritten in loop
% V3.2 (2023/07/21) :
%        - correct bug at the final zipping step when file does not exist
%          (in case i_bgc = 0).
%        - add case All Argo Fleet and manage optimization section for large number
%        of floats.
%        - dealing with float_wmo type for cases when not all WMOs are coded on 7 characters.
%        - dealing with graph layout when n_countries is large
%        - special workaround for float 4900566 that used QC 1 instead of QC A for
%        profile QC.
%        - add quotes in synthese output for program, in case comma is
%        used.
% V3.3 (2023/10/02) :
%        - add hhmmss in the output directory name
%        - change search for param name in index for a more robust means
%        - add x grid and minor grid for psal adjustment display by wmo.
%        - add an option to group prof QC A and B
% V3.4 (2024/02/23) :
%        - remove CTD from plot 09 when i_bgc=1.
%        - 1-yr static string was replaced by the config value
%        - add a new graph per profile year with DMQC and profile QC F + 
%          output values in a text file
%        - add a log (diary)
% V3.5 (2025/10/17) :
%        - modify psal adj plot to separate RA mode not QC-F from RA mode QC-F
%        profiles.
%        - correct cycles and descending profiles processing in
%        get_data_from_index.m routine
%        - replace "grey list" by "exclusion list"
%        - correct issue with print('-dpng') for more recent Matlab
%        versions, which do not support a space.
%        - externalise the configuration in a conf file
%        - choose the method to select old profiles/floats: either by profile age
%        in days or by the profile date

%option explicit

clear variables
close all

% inputs and configuration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp([newline '######################'])
disp('Load configuration')

cf=load_configuration('get_DMQC_status_config.txt');

test_case_dir=cf.test_case_dir;

wmo_list_file=[cf.test_case_dir '/' cf.input_file_dir '/' cf.wmo_list_file];
country_code_file=[cf.test_case_dir '/' cf.input_file_dir '/' cf.country_code_file];
index_file=[cf.index_file_dir '/' cf.index_file_core];
index_file_synthetic=[cf.index_file_dir '/' cf.index_file_bgc_synthetic];
exclusionlist_file=[cf.exclusion_file_dir '/' cf.exclusion_file];
output_file_dir_prefix=cf.output_file_dir_prefix;
script_path=cf.script_dir;
log_file=cf.log_file;


project_name=cf.project_name;

input_list_of_parameters_to_treat = string(split(cf.input_list_of_parameters_to_treat,";"));

i_bgc=str2double(cf.i_bgc);
input_list_of_BGC_parameters_to_treat = string(split(cf.input_list_of_BGC_parameters_to_treat,";"));

output_graphs_per_float = str2double(cf.output_graphs_per_float);
n_max_float_per_graph   = str2double(cf.n_max_float_per_graph);

i_descending_profile = str2double(cf.i_descending_profile);
profile_age_method   = cf.profile_age_method; % either 'days' or 'date'
profile_age_min_days = str2double(cf.profile_age_min_days);
profile_age_max_date = datenum(cf.profile_age_max_date,'yyyy/MM/dd');
print_svg            = str2double(cf.print_svg);
i_group_AB_profQC    = str2double(cf.i_group_AB_profQC);

wmo_list_file_short = [cf.test_case_dir '/' cf.input_file_dir '/wmo_list.txt'];
index_file_short=[cf.test_case_dir '/' cf.input_file_dir '/argo_profile_detailed_index_subset.txt'];
index_file_synthetic_short = [cf.test_case_dir '/' cf.input_file_dir '/argo_profile_detailled_index_subset.txt'];
if strcmp(profile_age_method,'days')
    profile_old_a= sprintf('%.1f',profile_age_min_days/365);
    profile_old_b=['older than ' profile_old_a ' yr'];
    profile_old_c=['older_than_' profile_old_a '_yr'];
end
if strcmp(profile_age_method,'date')
    profile_old_a= sprintf('%s',datestr(profile_age_max_date,'yyyy-mm-dd'));
    profile_old_b=['older than ' profile_old_a];
    profile_old_c=['older_than_' profile_old_a];
end


test_date=char(datetime('now','TimeZone','local','Format','yyyy-MM-dd_HHmmSS'));
output_dir=[ cf.test_case_dir cf.output_file_dir_prefix '_' test_date '/'];
if ~exist(output_dir, 'dir')
    disp('creating the output_directory')
    mkdir(output_dir)
end

log_file=[output_dir '/' cf.log_file];
   
if i_bgc == 1
    disp([newline ' Treating CTD and BGC parameters for ' project_name ' case']);
else
    disp([newline ' Treating CTD parameters for ' project_name ' case']);
end

% add paths (packages and auxiliary functions)
% aux_functions_path = [script_path '/aux_functions'];
addpath(genpath(script_path))

disp([newline 'End of configuration'])
disp('######################')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read input files
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp([newline '######################'])
disp('Reading input files')
tStart = tic;

% Copy the instance of the index file that was used for the computation
% This makes it available for later understanding of any bugs / or rerun 
% with data from the same date. For the seek of space: to be zipped once
% finished.

diary(log_file)
diary on

disp('-- creating a local copy of the index file ...')
local_index_file=[output_dir 'argo_profile_detailled_index_' test_date '.txt'];
copyfile(index_file,local_index_file)

if i_bgc == 1
    local_index_synthetic_file=[output_dir 'argo_synthetic-profile_detailled_index_' test_date '.txt'];
    copyfile(index_file_synthetic,local_index_synthetic_file)
end


local_exclusionlist_file=[output_dir 'ar_exclusionlist_' test_date '.txt'];
disp('-- creating a local copy of the exclusionlist file ...')
copyfile(exclusionlist_file,local_exclusionlist_file)


local_wmo_list_file=[output_dir 'wmo_list_processed_' test_date '.txt'];
disp('-- creating a local copy of the wmo_list file ...')
copyfile(wmo_list_file,local_wmo_list_file)


local_country_code_file=[output_dir 'country_code_' test_date '.txt'];
disp('-- creating a local copy of the country_code file ...')
copyfile(country_code_file,local_country_code_file)

% read wmo list of floats to be treated
disp('-- reading wmo list to treat ...')
[Floats_list] = read_csv(wmo_list_file,';');
n_floats = size(Floats_list.WMO,1);
float_list_str=cellstr(Floats_list.WMO);
float_list_str_wo_blank=strrep(float_list_str,' ','');
clearvars Floats_list.WMO 
Floats_list.WMO=float_list_str_wo_blank;
fprintf('-- number of floats : %d\n',n_floats)

% read country_codes
disp('-- reading country_codes mapping ...')
[Country_codes] = read_csv(country_code_file,';');


% map country_codes into floats_list structure:
l_min=min(size(Floats_list.COUNTRY,2),size(Floats_list.COUNTRY,2));
[~,loc]=ismember(string(Floats_list.COUNTRY(:,1:l_min)),string(Country_codes.COUNTRY(:,1:l_min)));
Floats_list.COUNTRYCODE=string(Country_codes.COUNTRYCODE(loc,:));

% read exclusion list
disp('-- reading exclusionlist ...')
[exclusion_list] = read_csv(local_exclusionlist_file,',');
disp('-- end of exclusionlist reading ...')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% optimizing sequence
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('-- optimizing index_file content ...')

data=cellstr(strcat(Floats_list.WMO,"/"))';
file_id = fopen(wmo_list_file_short, 'w');
fprintf(file_id, '%s\n', data{:});
fclose(file_id);

nb_header_lines_core=0;
nb_header_lines_bgc=0;

optim=1;
try

    if n_floats > 5000
        
        optim = 0;
        
    else
        disp('--- using unix commands')
        disp('--- (approx 86 seconds for 4500 floats)')
        disp(['--- ' char(datetime('now','TimeZone','local'))])
        tStart = tic;
        cmd = ['grep ^file ' local_index_file ' > ' index_file_short];
        system(cmd);
        cmd = ['grep -f ' wmo_list_file_short ' ' local_index_file ' >> ' index_file_short];
        system(cmd);

        if i_bgc ==1
            cmd = ['grep ^file ' local_index_synthetic_file ' > ' index_file_synthetic_short];
            system(cmd);
            cmd = ['grep -f ' wmo_list_file_short ' ' local_index_synthetic_file ' >> ' index_file_synthetic_short];
            system(cmd);
            cmd = ['wc -l ' index_file_synthetic_short '| awk ' char(39) '{print $1}' char(39)];
            [~,n_line]=system(cmd);
            n_line=double(string(n_line));
            if n_line == 1
                disp('There is no entry (BGC) in the synthetic index file for this list of float: back to i_bgc=0')
                i_bgc=0;
            end
        end
    end
    
    
catch
    % If unix commands are not available:
    % switching to an equivalent in Matlab
    % This equivalent is far less performant:
    %    Matlab vs  Unix for Load
    %  -   50 s vs  4 s  for  100 floats
    %  -  119 s vs  4 s  for  200 floats
    %  -  360 s vs  8 s  for  726 floats
    %  - 2400 s vs 86 s  for 4443 floats 
    % as the readtable function from get_data_from_index takes 300 second
    % for the whole Argo fleet, the optimization using equivalent
    % unix-function in Matlab will be performed only when n_floats < 200
    % N.B: the Matlab grep command has a different
    % order in outputs compared to the unix command.
    % once sorted, both outputs match.

    
    if n_floats > 200
        
        optim=0;
        
    else    
        disp('--- Optimization through Matlab unix-like commands')
        disp('--- (approx 2 minutes for 200 floats)')
        

        disp('---- header for core')
        [~,PH]=grep('-s','latitude',local_index_file);
        data=cellstr(PH.match);
        file_id = fopen(index_file_short, 'w');
        fprintf(file_id, '%s\n', data{:});
        fclose(file_id);


        disp('---- grep pattern file for core')
        disp(['---- ' char(datetime('now','TimeZone','local'))])
        % warning: option '-f' and pattern file (wmo_list_file_short) must
        % be consecutive in the arguments of the function.
        [~,P]=grep('-s','-f',wmo_list_file_short,local_index_file);

        data=cellstr(P.match);
        file_id = fopen(index_file_short, 'a');
        fprintf(file_id, '%s\n', data{:});
        fclose(file_id);


        if i_bgc ==1


            disp('---- header for synthetic')
            [~,PH]=grep('-s','file,date',local_index_synthetic_file);

            data=cellstr(PH.match);
            file_id = fopen(index_file_synthetic_short, 'w');
            fprintf(file_id, '%s\n', data{:});
            fclose(file_id);

            disp('---- grep pattern file for synthetic')
            disp(['---- ' char(datetime('now','TimeZone','local'))])
            [~,P]=grep('-s','-f',wmo_list_file_short,local_index_synthetic_file);

            if size(P.match) > 0
                data=cellstr(P.match);
                file_id = fopen(index_file_synthetic_short, 'a');
                fprintf(file_id, '%s\n', data{:});
                fclose(file_id);
            else
                disp('There is no entry (BGC) in the synthetic index file for this list of float: back to i_bgc=0')
                i_bgc=0;
            end
        end
        
        clearvars P PH
    end
    
end

% delete temporary file wmo_list_file_short
delete(wmo_list_file_short);

if optim == 0
    disp('--- too large number of floats requested for index file optimization')
    disp('--- the whole index will be read - approx 5 minutes')
    %tStart = tic;
    [~,P]=grep('-s','#',local_index_file); % 3 seconds
    %tEnd = toc(tStart);
    nb_header_lines_core=P.lcount;
    copyfile(local_index_file,index_file_short)

    if i_bgc ==1
        %tStart = tic;
        [~,P]=grep('-s','#',local_index_synthetic_file); % 0.3 seconds
        %tEnd = toc(tStart);
        nb_header_lines_bgc=P.lcount;
        copyfile(local_index_synthetic_file,index_file_synthetic_short)
    end

    clearvars P
end

tEnd = toc(tStart);
disp (['-- End of optimizing index_file content (' char(string(floor(tEnd))) ' s)'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% end of index optimizing sequence
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Read index file
tStart = tic;
% Performances are approximately the following on datarmor (and using 15 GB of RAM):
%       -  23s for   591 000 lines (4443 floats).
%       - 300s for 2 832 660 lines (all argo floats at 2023/07/13)
%       - 391s for 3 500 000 lines (+ 23 % ~ + 4   yr - scalability test)
%       - 485s for 4 000 000 lines (+ 41 % ~ + 6.7 yr - scalability test)
%      unknown for 5 665 311 lines (+100 % ~ +16   yr => failed with mem=15g)
% The time estimates associated to the increase in line numbers assumes:
% + 4700 * 36 = + 169 200 profiles per year (OneArgo target).


disp('-- Entering get_data_from_index for CTD parameters')
disp(['-- ' char(datetime('now','TimeZone','local'))])
[IndexData_CTD] = get_data_from_index(index_file_short, ...
                                      nb_header_lines_core, ...
                                      float_list_str_wo_blank, ...
                                      i_descending_profile, ...
                                      0, ...
                                      input_list_of_parameters_to_treat);
list_of_CTD_parameters_to_treat=IndexData_CTD.ParamList;

if i_bgc ==1
    disp('-- Entering get_data_from_index for BGC parameters')
    disp(['-- ' char(datetime('now','TimeZone','local'))])
    [IndexData_BGC] = get_data_from_index(index_file_synthetic_short, ...
                                          nb_header_lines_bgc, ...
                                          float_list_str_wo_blank, ...
                                          i_descending_profile,...
                                          i_bgc,...
                                          input_list_of_BGC_parameters_to_treat);
    list_of_BGC_parameters_to_treat=IndexData_BGC.ParamList;
    
% delete temporary file short
delete(index_file_short)
if i_bgc ==1
    delete(index_file_synthetic_short)
end
    
    % DEBUG
    % save PSAL QC info.
%     data_to_record=cat(2,IndexData_CTD.profile_WMO, IndexData_CTD.cycle, IndexData_CTD.param.PSAL.qc);
%     file= [ test_case_dir '/dbg_psal_qc_from_index.csv'];
%     data_to_record=array2table(data_to_record);
%     writetable(data_to_record,file);
%     
%     data_to_record=cat(2,IndexData_BGC.profile_WMO, IndexData_BGC.cycle, IndexData_BGC.param.PSAL.qc);
%     file= [ test_case_dir '/dbg_psal_qc_from_synthetic_index.csv'];
%     data_to_record=array2table(data_to_record);
%     writetable(data_to_record,file);
    
    
    disp('-- Merging information')
    disp('--- merging level1 information')
    list_of_parameters_to_treat=cat(1,list_of_CTD_parameters_to_treat,list_of_BGC_parameters_to_treat);
    nb_lines_ctd=size(IndexData_CTD.profile_WMO,1);
    nb_lines_bgc=size(IndexData_BGC.profile_WMO,1);
    nb_lines_ctd_bgc=nb_lines_ctd+nb_lines_bgc;
    IndexData.source=strings(nb_lines_ctd_bgc,1);
    IndexData.source(1:size(IndexData_CTD.profile_WMO,1))="CTD";
    IndexData.source(size(IndexData_CTD.profile_WMO,1)+1:end)="BGC";
    
    IndexData.profile_WMO=cat(1,IndexData_CTD.profile_WMO,IndexData_BGC.profile_WMO);
    IndexData.profile_dac=cat(1,IndexData_CTD.profile_dac,IndexData_BGC.profile_dac);
    IndexData.cycle=cat(1,IndexData_CTD.cycle,IndexData_BGC.cycle);
    IndexData.date=cat(1,IndexData_CTD.date,IndexData_BGC.date);
    IndexData.latitude=cat(1,IndexData_CTD.latitude,IndexData_BGC.latitude);
    IndexData.longitude=cat(1,IndexData_CTD.longitude,IndexData_BGC.longitude);
    IndexData.prof_update_date=cat(1,IndexData_CTD.prof_update_date,IndexData_BGC.prof_update_date);
    IndexData.comment=IndexData_CTD.comment;
    IndexData.ParamList=cat(1,IndexData_CTD.ParamList,IndexData_BGC.ParamList);
    IndexData.ParamNotFound_In_BGC_Index_Requested=IndexData_BGC.ParamNotFoundInIndexList;
    IndexData.ParamFound_In_BGC_Index_notRequested=IndexData_BGC.ParamFoundInIndex_notRequested;
    IndexData.WMOs_unfound_in_IndexFile.CTD=IndexData_CTD.WMOs_unfound_in_IndexFile;
    IndexData.WMOs_unfound_in_IndexFile.BGC=IndexData_BGC.WMOs_unfound_in_IndexFile;
    
    disp('--- merging level2 information')
    for i=1:size(list_of_parameters_to_treat,1)
        i_param=list_of_parameters_to_treat(i);
        
        % initialisation
        IndexData.param.(i_param).qc = strings(nb_lines_ctd_bgc,1);
        IndexData.param.(i_param).mode = strings(nb_lines_ctd_bgc,1);
        IndexData.param.(i_param).presence = zeros(nb_lines_ctd_bgc,1);
        if i_param == "PSAL"
            IndexData.param.(i_param).adj = zeros(nb_lines_ctd_bgc,1);
        end
        
    end
    
    for i=1:size(list_of_CTD_parameters_to_treat,1)
        i_param=list_of_CTD_parameters_to_treat(i);
        
        IndexData.param.(i_param).nb_profile=IndexData_CTD.param.(i_param).nb_profile;
        IndexData.param.(i_param).qc(1:nb_lines_ctd)=IndexData_CTD.param.(i_param).qc;
        IndexData.param.(i_param).mode(1:nb_lines_ctd)=IndexData_CTD.param.(i_param).mode;
        IndexData.param.(i_param).presence(1:nb_lines_ctd)=IndexData_CTD.param.(i_param).presence;
        if i_param == "PSAL"
            IndexData.param.(i_param).adj(1:nb_lines_ctd)=IndexData_CTD.param.(i_param).adj;
        end
    end
    
    
    for i=1:size(list_of_BGC_parameters_to_treat,1)
        i_param=list_of_BGC_parameters_to_treat(i);
        
        IndexData.param.(i_param).nb_profile=IndexData_BGC.param.(i_param).nb_profile;
        IndexData.param.(i_param).qc(nb_lines_ctd+1:end)=IndexData_BGC.param.(i_param).qc;
        IndexData.param.(i_param).mode(nb_lines_ctd+1:end)=IndexData_BGC.param.(i_param).mode;
        IndexData.param.(i_param).presence(nb_lines_ctd+1:end)=IndexData_BGC.param.(i_param).presence;
    end
    
    
else
    list_of_parameters_to_treat=list_of_CTD_parameters_to_treat;
    IndexData=IndexData_CTD;
    IndexData.source=strings(size(IndexData_CTD.profile_WMO));
    IndexData.source(1:end)="CTD";
    %IndexData.WMOs_unfound_in_IndexFile.CTD=IndexData_CTD.WMOs_unfound_in_IndexFile;
end

clearvars IndexData_CTD IndexData_BGC


disp('-- End of get_data_from_index')

n_param=size(list_of_parameters_to_treat,1);

% Get the update date of the index file
try
    cmd = ['grep "Date of update" ' local_index_file ' | awk ' char(39) '{print $NF}' char(39)];
    [status,IndexData.index_update]=system(cmd);
catch
    [~,P]=grep('-s',' Date of update', local_index_file);
    tmp=split(P.match," ");
    IndexData.index_update=char(tmp(end));
end
update_date_str=datestr(datenum(IndexData.index_update,'yyyymmddHHMMSS'),'yyyy-mm-dd');

% Creating output directory:
working_date = IndexData.index_update(1:8);


tEnd = toc(tStart);
disp([newline 'End of reading index files (' char(string(floor(tEnd))) ' s)'])
disp('######################')
diary off
%%
diary on

% Compute values for graphs:
disp('')
disp('######################')
disp('A few computation before making graphs')
tStart = tic;

% Attribute the country_code within the IndexData structure:
[~,loc]=ismember(IndexData.profile_WMO,Floats_list.WMO);
IndexData.country=string(Floats_list.COUNTRYCODE(loc,:));
% Compute the numbers of countries:
n_countries = size(unique(string(Floats_list.COUNTRYCODE)),1);

% Compute profile age when the index file was updated

% Format to date type
index_file_date = datenum(IndexData.index_update,'yyyymmddHHMMSS');
prof_date = datenum(IndexData.date,'yyyymmddHHMMSS');
prof_update_date = datenum(IndexData.prof_update_date,'yyyymmddHHMMSS');
% Compute the delta time between index file date and each profile's date
IndexData.profile_age=index_file_date-prof_date;
% Index the profiles that are considered as old enough
if strcmp(profile_age_method,'days')
    i_old_profiles = (IndexData.profile_age > profile_age_min_days);
end
if strcmp(profile_age_method,'date')
    i_old_profiles = (prof_date < profile_age_max_date);
end



for i=1:n_param
    i_param=list_of_parameters_to_treat(i);
    fprintf('Indexing %s profile, mode and operational wmos information\n',i_param);
    
    % Index the wmos that are still operational (to compute exclusionlist
    % proportion)
    i_recent_profiles=(IndexData.profile_age < 30);
    wmos_operational.(i_param)=unique(IndexData.profile_WMO(i_recent_profiles & IndexData.param.(i_param).presence==1));
    nb_operational_floats.(i_param)=size(wmos_operational.(i_param),1);
    
    % Index the profiles that were DMQCed
    i_DMQCed.(i_param)=(IndexData.param.(i_param).mode == 'D');
    i_R_or_A_profile.(i_param)=(IndexData.param.(i_param).mode == 'R' | IndexData.param.(i_param).mode == 'A');
    wmos_with_R_or_A_profile.(i_param)=unique(IndexData.profile_WMO(i_R_or_A_profile.(i_param)));
    
    % Index the profiles that have profile_QC = 'F'
    i_Fprof.(i_param)=(IndexData.param.(i_param).qc == 'F');
    
    % Index the profiles that have param
    i_pres.(i_param)=(IndexData.param.(i_param).presence == 1);

end


% Figure 1 and 2: data mode in number of floats, in number of profiles (per country)
% ---------------------------------

% data computation
disp('compute DMQC status per country by float and by profile')

% Number of floats per country
disp('- compute nb of floats per country')
% Define one unique x axis:
nb_per_country_x=unique(string(Floats_list.COUNTRYCODE));
for i=1:n_param
    i_param=list_of_parameters_to_treat(i);
    
    %output initialisation
    nb_floats_per_country.(i_param)=zeros(n_countries,1);
    nb_floats_xage_per_country.(i_param)=zeros(n_countries,1);
    nb_floats_DMQCed_per_country.(i_param)=zeros(n_countries,1);
    nb_floats_xage_DMQCed_per_country.(i_param)=zeros(n_countries,1);
    
    
    fprintf('-- compute nb of floats for %s per country\n',i_param)
    wmos.param.(i_param)=unique(IndexData.profile_WMO((IndexData.param.(i_param).presence==1))); %retrieve unique wmos for which there is a profile with the parameter
    [~,loc]=ismember(wmos.param.(i_param),Floats_list.WMO);
    wmos_country=string(Floats_list.COUNTRYCODE(loc,:)); % associate the country_code
    [xx,~,ic]=unique(wmos_country);
    val=accumarray(ic,1);
    [~,loc]=ismember(xx,nb_per_country_x);
    nb_floats_per_country.(i_param)(loc)=val;
    nb_float_tot.(i_param)=sum(val);

    fprintf('-- compute nb of floats with at least 1 profile %s for %s per country\n',profile_old_b,i_param)
    wmos_xage.param.(i_param)=unique(IndexData.profile_WMO(i_old_profiles & (IndexData.param.(i_param).presence==1)));
    [~,loc]=ismember(wmos_xage.param.(i_param),Floats_list.WMO);
    wmos_xage_country=string(Floats_list.COUNTRYCODE(loc,:));
    [xx,~,ic]=unique(wmos_xage_country);
    val=accumarray(ic,1);
    [~,loc]=ismember(xx,nb_per_country_x);
    nb_floats_xage_per_country.(i_param)(loc)=val;
    nb_float_xage_tot.(i_param)=sum(val);
    
    
    fprintf('-- compute nb of floats with at least 1 profile DMQCed for %s per country\n',i_param)
    wmos_DMQCed=unique(IndexData.profile_WMO(i_DMQCed.(i_param)));
    [~,loc]=ismember(wmos_DMQCed,Floats_list.WMO);
    wmos_dmqced_country=string(Floats_list.COUNTRYCODE(loc,:));
    [xx,~,ic]=unique(wmos_dmqced_country);
    val=accumarray(ic,1);
    [~,loc]=ismember(xx,nb_per_country_x);
    nb_floats_DMQCed_per_country.(i_param)(loc)=val;
    nb_float_DMQCed_tot.(i_param)=sum(val);

    fprintf('-- compute nb of floats with at least 1 profile DMQCed for %s and with at least 1 profile %s per country\n',i_param,profile_old_b)
    wmos_xage_DMQCed=unique(IndexData.profile_WMO(i_DMQCed.(i_param) & i_old_profiles));
    [~,loc]=ismember(wmos_xage_DMQCed,Floats_list.WMO);
    wmos_xage_dmqced_country=string(Floats_list.COUNTRYCODE(loc,:));
    [xx,~,ic]=unique(wmos_xage_dmqced_country);
    val=accumarray(ic,1);
    [~,loc]=ismember(xx,nb_per_country_x);
    nb_floats_xage_DMQCed_per_country.(i_param)(loc)=val;
    nb_floats_xage_DMQCed_tot.(i_param)=sum(val);
    
end


% Number of profile per country
disp('- compute nb of profiles per country')
for i=1:n_param
    i_param=list_of_parameters_to_treat(i);
    
    
    %output initialisation
    nb_profiles_per_country.(i_param)=zeros(n_countries,1);
    nb_profiles_xage_per_country.(i_param)=zeros(n_countries,1);
    nb_profiles_DMQCed_per_country.(i_param)=zeros(n_countries,1);
    nb_profiles_xage_DMQCed_per_country.(i_param)=zeros(n_countries,1);
    
    
    fprintf('-- compute nb of profiles for %s per country\n',i_param)
    [xx,~,ic]=unique(IndexData.country(IndexData.param.(i_param).presence==1));
    val=accumarray(ic,1);
    [~,loc]=ismember(xx,nb_per_country_x);
    nb_profiles_per_country.(i_param)(loc)=val;

    fprintf('-- compute nb of profiles %s for %s per country\n',profile_old_b,i_param)
    [xx,~,ic]=unique(IndexData.country(i_old_profiles & IndexData.param.(i_param).presence==1));
    val=accumarray(ic,1);
    [~,loc]=ismember(xx,nb_per_country_x);
    nb_profiles_xage_per_country.(i_param)(loc)=val;

    fprintf('-- compute nb of profiles DMQCed for %s per country\n',i_param)
    [xx,~,ic]=unique(IndexData.country(i_DMQCed.(i_param)));
    val=accumarray(ic,1);
    [~,loc]=ismember(xx,nb_per_country_x);
    nb_profiles_DMQCed_per_country.(i_param)(loc)=val;

    fprintf('-- compute nb of profiles DMQCed and %s for %s per country\n',profile_old_b, i_param)
    [xx,~,ic]=unique(IndexData.country(i_DMQCed.(i_param) & i_old_profiles));
    val=accumarray(ic,1);
    [~,loc]=ismember(xx,nb_per_country_x);
    nb_profiles_xage_DMQCed_per_country.(i_param)(loc)=val;
end

% Figure 3 : Quality control per parameter
% -----------------------------------------
disp('- compute profile quality control stats')
nb_per_qc_x=string(transpose('ABCDEFX'));
for i=1:n_param
    i_param=list_of_parameters_to_treat(i);
    
    nb_profiles_per_qc.(i_param)=zeros(7,1);
    nb_profiles_DMQCed_per_qc.(i_param)=zeros(7,1);
    
    IndexData.param.(i_param).qc(IndexData.param.(i_param).qc==" " & ...
                                 IndexData.param.(i_param).presence==1)='X';
    IndexData.param.(i_param).qc(IndexData.param.(i_param).qc==""  & ...
                                 IndexData.param.(i_param).presence==1)='X';
    
    % workaround for float 4900566 that filled profile QC with "1" instead
    % of "A" (DM file were last updated in 2013/2014).
    IndexData.param.(i_param).qc(IndexData.param.(i_param).qc=="1")='A';
    
    fprintf('-- compute nb of %s profiles per qc_code\n',i_param)
    [xx,~,ic]=unique(IndexData.param.(i_param).qc(IndexData.param.(i_param).presence==1));
    val=accumarray(ic,1);
    [~,loc]=ismember(xx,nb_per_qc_x);
    nb_profiles_per_qc.(i_param)(loc)=val;
    
    fprintf('-- compute nb of %s D-profiles per qc_code\n',i_param)
    [xx,~,ic]=unique(IndexData.param.(i_param).qc(IndexData.param.(i_param).presence==1 & ...
                                                  string(IndexData.param.(i_param).mode)=="D"     ));
    val=accumarray(ic,1);
    [~,loc]=ismember(xx,nb_per_qc_x);
    nb_profiles_DMQCed_per_qc.(i_param)(loc)=val;
    
end

% Figure 4 : Quality control evolution per parameter per launch year
% ------------------------------------------------------------------
disp('- compute profile quality control stats time evolution')
for i=1:n_param
    i_param=list_of_parameters_to_treat(i);
    % First match profile info with launch_years:
    [~,loc]=ismember(IndexData.profile_WMO,Floats_list.WMO);
    wmos_launch_year=string(Floats_list.LAUNCHDATE(loc,1:4));

    % Define one unique x axis:
    nb_per_launch_year_x=unique(wmos_launch_year);
    n_launch_years=size(nb_per_launch_year_x,1);

    %Initatisations:
    percent_prof_QC_A_per_launch_year.(i_param)=zeros(n_launch_years,1);
    percent_prof_QC_B_per_launch_year.(i_param)=zeros(n_launch_years,1);
    percent_prof_QC_C_per_launch_year.(i_param)=zeros(n_launch_years,1);
    percent_prof_QC_D_per_launch_year.(i_param)=zeros(n_launch_years,1);
    percent_prof_QC_E_per_launch_year.(i_param)=zeros(n_launch_years,1);
    percent_prof_QC_F_per_launch_year.(i_param)=zeros(n_launch_years,1);

    percent_prof_DMQCed_QC_A_per_launch_year.(i_param)=zeros(n_launch_years,1);
    percent_prof_DMQCed_QC_B_per_launch_year.(i_param)=zeros(n_launch_years,1);
    percent_prof_DMQCed_QC_C_per_launch_year.(i_param)=zeros(n_launch_years,1);
    percent_prof_DMQCed_QC_D_per_launch_year.(i_param)=zeros(n_launch_years,1);
    percent_prof_DMQCed_QC_E_per_launch_year.(i_param)=zeros(n_launch_years,1);
    percent_prof_DMQCed_QC_F_per_launch_year.(i_param)=zeros(n_launch_years,1);
    
    nb_prof_per_launch_year.(i_param)=zeros(n_launch_years,1);
    nb_prof_DMQCed_per_launch_year.(i_param)=zeros(n_launch_years,1);


    % Then compute percentage of A, B, C, D, E, F by launch years:
    for iyear = 1:n_launch_years
       cur_year= nb_per_launch_year_x(iyear);

       i_prof_launch_year=(wmos_launch_year==cur_year & IndexData.param.(i_param).presence==1);
       nb_prof_per_launch_year.(i_param)(iyear)=sum(i_prof_launch_year);
       percent_prof_QC_A_per_launch_year.(i_param)(iyear)=100*sum(i_prof_launch_year & IndexData.param.(i_param).qc == "A")/sum(i_prof_launch_year);
       percent_prof_QC_B_per_launch_year.(i_param)(iyear)=100*sum(i_prof_launch_year & IndexData.param.(i_param).qc == "B")/sum(i_prof_launch_year);
       percent_prof_QC_C_per_launch_year.(i_param)(iyear)=100*sum(i_prof_launch_year & IndexData.param.(i_param).qc == "C")/sum(i_prof_launch_year);
       percent_prof_QC_D_per_launch_year.(i_param)(iyear)=100*sum(i_prof_launch_year & IndexData.param.(i_param).qc == "D")/sum(i_prof_launch_year);
       percent_prof_QC_E_per_launch_year.(i_param)(iyear)=100*sum(i_prof_launch_year & IndexData.param.(i_param).qc == "E")/sum(i_prof_launch_year);
       percent_prof_QC_F_per_launch_year.(i_param)(iyear)=100*sum(i_prof_launch_year & IndexData.param.(i_param).qc == "F")/sum(i_prof_launch_year);

       i_prof_DMQCed_launch_year=(wmos_launch_year==cur_year & ...
                                  IndexData.param.(i_param).presence==1 & ...
                                  string(IndexData.param.(i_param).mode)=="D");
       nb_prof_DMQCed_per_launch_year.(i_param)(iyear)=sum(i_prof_DMQCed_launch_year);
       percent_prof_DMQCed_QC_A_per_launch_year.(i_param)(iyear)=100*sum(i_prof_DMQCed_launch_year & IndexData.param.(i_param).qc == "A")/sum(i_prof_DMQCed_launch_year);
       percent_prof_DMQCed_QC_B_per_launch_year.(i_param)(iyear)=100*sum(i_prof_DMQCed_launch_year & IndexData.param.(i_param).qc == "B")/sum(i_prof_DMQCed_launch_year);
       percent_prof_DMQCed_QC_C_per_launch_year.(i_param)(iyear)=100*sum(i_prof_DMQCed_launch_year & IndexData.param.(i_param).qc == "C")/sum(i_prof_DMQCed_launch_year);
       percent_prof_DMQCed_QC_D_per_launch_year.(i_param)(iyear)=100*sum(i_prof_DMQCed_launch_year & IndexData.param.(i_param).qc == "D")/sum(i_prof_DMQCed_launch_year);
       percent_prof_DMQCed_QC_E_per_launch_year.(i_param)(iyear)=100*sum(i_prof_DMQCed_launch_year & IndexData.param.(i_param).qc == "E")/sum(i_prof_DMQCed_launch_year);
       percent_prof_DMQCed_QC_F_per_launch_year.(i_param)(iyear)=100*sum(i_prof_DMQCed_launch_year & IndexData.param.(i_param).qc == "F")/sum(i_prof_DMQCed_launch_year);


    end
end


% Figure 5 : Operational floats in exclusion list for PSAL
% ---------------------------------------------------
disp('- compute exclusion list stats')
for i=1:n_param
    i_param=list_of_parameters_to_treat(i);
    fprintf('-- for %s \n',i_param)

    i_wmo_ope_in_exclusionlist_with_qc_3_or_4.(i_param)  = contains(string(exclusion_list.PARAMETERNAME),i_param) & ...
                                                 (exclusion_list.QUALITYCODE == '3' | exclusion_list.QUALITYCODE == '4') & ...
                                                 ismember(string(exclusion_list.PLATFORMCODE),string(wmos_operational.(i_param)));
    nb_wmo_ope_in_exclusionlist_with_qc_3_or_4.(i_param) = sum(i_wmo_ope_in_exclusionlist_with_qc_3_or_4.(i_param));

    i_wmo_inactiveR_in_exclusionlist_with_qc_3_or_4.(i_param)  = contains(string(exclusion_list.PARAMETERNAME),i_param) & ...
                                                 (exclusion_list.QUALITYCODE == '3' | exclusion_list.QUALITYCODE == '4') & ...
                                                 ismember(string(exclusion_list.PLATFORMCODE),string(wmos_with_R_or_A_profile.(i_param))) & ...
                                                 ~ismember(string(exclusion_list.PLATFORMCODE),string(wmos_operational.(i_param)));
    nb_wmo_inactiveR_in_exclusionlist_with_qc_3_or_4.(i_param) = sum(i_wmo_inactiveR_in_exclusionlist_with_qc_3_or_4.(i_param));


    i_wmo_inactiveD_in_exclusionlist_with_qc_3_or_4.(i_param)  = contains(string(exclusion_list.PARAMETERNAME),i_param) & ...
                                                 (exclusion_list.QUALITYCODE == '3' | exclusion_list.QUALITYCODE == '4') & ...
                                                 ismember(string(exclusion_list.PLATFORMCODE),string(Floats_list.WMO)) & ...
                                                 ~ismember(string(exclusion_list.PLATFORMCODE),string(wmos_with_R_or_A_profile.(i_param))) & ...
                                                 ~ismember(string(exclusion_list.PLATFORMCODE),string(wmos_operational.(i_param)));
    nb_wmo_inactiveD_in_exclusionlist_with_qc_3_or_4.(i_param) = sum(i_wmo_inactiveD_in_exclusionlist_with_qc_3_or_4.(i_param));

end

% Figure 6 : DMQC status w.r.t profile year
% ---------------------------------------------------
disp('- compute DMQC status w.r.t profile year')
% Define one unique x axis:
IndexData.date(ismissing(IndexData.date))='19800106000000';
tmp=char(IndexData.date);
nb_per_year_x=unique(string(tmp(:,1:4)));
n_years=size(nb_per_year_x,1);

for i=1:n_param
    i_param=list_of_parameters_to_treat(i);
    fprintf('-- for %s \n',i_param)
    
    %output initialisation
    nb_profiles_per_year.(i_param)=zeros(n_years,1);
    nb_profiles_DMQCed_per_year.(i_param)=zeros(n_years,1);
    nb_profiles_QCF_per_year.(i_param)=zeros(n_years,1);
    nb_profiles_DMQCed_QCF_per_year.(i_param)=zeros(n_years,1);
    
    % Extract year from the profile date
    %i_missing_profile_date=ismissing(IndexData.date);
    %nb_missing_profile_date=sum(i_missing_profile_date);
    
    IndexData.profile_year.(i_param)=string(1980*ones(size(IndexData.date,1),1));
    
    i_not_missing_profile_date=IndexData.param.(i_param).presence & ~ismissing(IndexData.date);

    if sum(i_not_missing_profile_date) > 0
        tmp=char(IndexData.date(i_not_missing_profile_date));
        IndexData.profile_year.(i_param)(i_not_missing_profile_date)=string(tmp(:,1:4));
    end

    % compute nb of profiles per profile year
    [xx,~,ic]=unique(IndexData.profile_year.(i_param)(i_pres.(i_param)));
    val=accumarray(ic,1);
    [~,loc]=ismember(xx,nb_per_year_x);
    nb_profiles_per_year.(i_param)(loc)=val;
    
    % compute nb of D-profile per profile year
    [xx,~,ic]=unique(IndexData.profile_year.(i_param)(i_DMQCed.(i_param) & i_pres.(i_param)));
    val=accumarray(ic,1);
    [~,loc]=ismember(xx,nb_per_year_x);
    nb_profiles_DMQCed_per_year.(i_param)(loc)=val;
    
    % compute nb of F-profile per profile year
    [xx,~,ic]=unique(IndexData.profile_year.(i_param)(i_Fprof.(i_param) & i_pres.(i_param)));
    val=accumarray(ic,1);
    [~,loc]=ismember(xx,nb_per_year_x);
    nb_profiles_QCF_per_year.(i_param)(loc)=val;
    
    % compute nb of DF-profile per profile year
    [xx,~,ic]=unique(IndexData.profile_year.(i_param)(i_DMQCed.(i_param) & i_Fprof.(i_param) & i_pres.(i_param)));
    val=accumarray(ic,1);
    [~,loc]=ismember(xx,nb_per_year_x);
    nb_profiles_DMQCed_QCF_per_year.(i_param)(loc)=val;
    
    
end

%nb_profiles_DMQCed_per_year.(i_param)
% nb_profiles_per_year.(i_param)

% Figure 6 : DMQC status w.r.t profile year
% -------------------------------------------
% No need to to special computations for this one.

tEnd = toc(tStart);
disp(['End of computations for graphics (' char(string(floor(tEnd))) ' s)'])
disp('######################')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

diary off
%% Make graphics
diary on

close all

disp(' ')
disp('Making graphics...')

if i_descending_profile== 0
    prof_included = '  [only ascending profiles]';
else
    prof_included = '  [including descending profiles]';
end

% create output folder
output_plots_dir = [output_dir '/Plots' ];

if ~exist(output_plots_dir, 'dir')
    disp('creating the output_directory')
    mkdir(output_plots_dir)
end
disp(['  plots will be saved in '  output_plots_dir])

% colors (https://www.rapidtables.org/fr/web/color/html-color-codes.html)
bars_colors = [0.5273    0.8047    0.9180; ... % 1  - total obs/floats (dark blue) old color: 0.2422    0.1504    0.6603
               0.9769    0.9839    0.0805; ... % 2  - >1year obs/floats (yellow)
               0.1953    0.8008    0.1953; ... % 3  - DMQC done (green)
               0.5938    0.9833    0.5938; ... % 4  - >1year obs/floats + DMQC done (pale green)
               0.5430    0         0.5430; ... % 5  - qc A "violet"
               0.9375    0.5000    0.5000; ... % 6  - qc B "saumon"
               1.0000    0.6470         0; ... % 7  - qc C "orange"
               0.8516    0.6445    0.1250; ... % 8  - qc D "verge d'or"
               0.1250    0.6953    0.6641; ... % 9  - qc E "vert clair"
               0.5273    0.8047    0.9180; ... % 10 - qc F "bleu ciel"
               0.4102    0.4102    0.4102; ... % 11 - exclusion list for operational floats AND qc X (no profile qc in index file)
               0.3000    0.3000    0.3000; ... % 12 - exclusion list for inactive floats with R or A profiles
               0.2000    0.2000    0.2000; ... % 13 - exclusion list for inactive floats with only D profiles
               0.2940    0         0.5090; ... % 14 - qc A for D profiles "indigo"
               0.8039    0.3608    0.3608; ... % 15 - qc B for D profiles "indianred"
               1.0000    0.5469         0; ... % 16 - qc C for D profiles "orange sombre"
               0.6863    0.5059    0.0157; ... % 17 - qc D for D profiles "verge d'or" with darker coefficients (not standard)
                    0    0.5020    0.5020; ... % 18 - qc E for D profiles "sarcelle"
               0.1176    0.5647    1.0000; ... % 19 - qc F for D profiles "dodgerblue"
               0.2000    0.2000    0.2000; ... % 20 - qc X for D profiles 
               255/255   69/255     0/255; ... % 21 - qc F (perprofyear)                "Orange Red"
               174/255   12/255     0/255];    % 22 - qc F (perprofyear) for D profiles "Mordant Red 19"


diary off
%% 
diary on
tStart = tic;
for i=1:n_param
    close all
 
    i_param=list_of_parameters_to_treat(i);
    
    if ismember("PSAL",list_of_parameters_to_treat)
        if (i_param == "TEMP" || i_param == "PRES") 
            % no need to output plot as mode is the same for TEMP, PRES,
            % PSAL
            continue
        else
            if i_param == "PSAL"
                i_param_str="CTD";
            else
                i_param_str=i_param;
            end
        end
    else
        i_param_str=i_param;
    end
    
    
    fprintf('Making Plots for %s \n',i_param_str)
    
    close all
           
    %%%%%%%%%%%%%%  Nb of float w.r.t DMQC status %%%%%%%%%%%%%% 
    disp('Nb of float w.r.t DMQC status')
    i_fig=1;
    figure(i_fig)
    % bigger figure
    set(gcf, 'Position', [200, 200, 1000, 600])
    % figure name
    set(gcf,'Name','Float DMQC status by country')

    
    stackData = cat(3,[nb_floats_DMQCed_per_country.(i_param) nb_floats_xage_DMQCed_per_country.(i_param)], ...
                      [(nb_floats_per_country.(i_param)-nb_floats_DMQCed_per_country.(i_param)) ...
                       (nb_floats_xage_per_country.(i_param)-nb_floats_xage_DMQCed_per_country.(i_param))]);
    hndl = plotBarStackGroups(stackData, nb_per_country_x);


    % FIGURE FORMAT
    hold on
    hp = plot(1,0,'w'); % for comment in legend
    % bar colors 
    set(hndl(1,1),'facecolor',bars_colors(3,:))
    set(hndl(1,2),'facecolor',bars_colors(1,:))
    set(hndl(2,1),'facecolor',bars_colors(4,:))
    set(hndl(2,2),'facecolor',bars_colors(2,:))
    % xlabels
    set(gca,'xtick',1:n_countries,'xticklabel', nb_per_country_x)
    if n_countries > 15
        set(gca,'XTickLabelRotation',90)
    end
    % title with update date
    title(['Float DMQC status for ' char(i_param_str) ' by country (updated ',update_date_str,')'], 'Interpreter', 'none')
    ylabel('Number of floats')
    % legend with total number
    if strcmp(project_name,'European_Fleet')
        legend([hndl(1,2), hndl(1,1), hndl(2,2), hndl(2,1), hp], ...
            {['Number of floats (Total: ',                       num2str(sum(nb_floats_per_country.(i_param),'omitnan')),')'],...
             ['Number of floats with D-profiles (Total: ',       num2str(sum(nb_floats_DMQCed_per_country.(i_param),'omitnan')),')'], ...
             ['Number of floats with profiles ' profile_old_b ' (Total: ',  num2str(sum(nb_floats_xage_per_country.(i_param),'omitnan')),')'],...
             ['Number of floats with D-profiles ' profile_old_b ' (Total: ',num2str(sum(nb_floats_xage_DMQCed_per_country.(i_param),'omitnan')),')'], ...
            prof_included})
    else
        legend([hndl(1,2), hndl(1,1), hndl(2,2), hndl(2,1), hp], ...
            {['Nb of floats to be managed at Euro-Argo (Total: ',                       num2str(sum(nb_floats_per_country.(i_param),'omitnan')),')'],...
             ['Nb of floats already addressed once in DMQC (Total: ',       num2str(sum(nb_floats_DMQCed_per_country.(i_param),'omitnan')),')'], ...
             ['Nb of floats for which DMQC can now be performed (' profile_old_b ') (Total: ',  num2str(sum(nb_floats_xage_per_country.(i_param),'omitnan')),')'],...
             ['Nb of floats already addressed once in DMQC (' profile_old_b ') (Total: ',num2str(sum(nb_floats_xage_DMQCed_per_country.(i_param),'omitnan')),')'], ...
            })
    end
    
        
        
    % background color
    set(gcf,'color','w');
    % grid in y axis
    ax = gca;
    ax.YGrid = 'on';
    
    ymax=max(nb_floats_per_country.(i_param));
    set(gca,'YLim',[0 ymax+ymax/4]);
    if n_countries <=15
        eps=0.02;
        for i_country =1:n_countries
            text(i_country-eps , nb_floats_per_country.(i_param)(i_country) + ymax/30, ...
                num2str(nb_floats_per_country.(i_param)(i_country)) , ...
                'HorizontalAlignment', 'right');
            text(i_country+eps, nb_floats_xage_per_country.(i_param)(i_country) + ymax/30, ...
                num2str(nb_floats_xage_per_country.(i_param)(i_country)), ...
                'HorizontalAlignment', 'left');
        end
    end

    % annotation: percentage of floats with R-profiles > 1 year
    floats_tobedone = (sum(nb_floats_xage_per_country.(i_param),'omitnan') - ... 
                       sum(nb_floats_xage_DMQCed_per_country.(i_param),'omitnan')) / ...
                           sum(nb_floats_xage_per_country.(i_param),'omitnan')*100;
    floats_tobedone_str = num2str(round(floats_tobedone));
    H = figure(1);
    set(H,'units','pix')
    annotation('textbox', [0.8, 0.01, .1, .1], 'string', ['Floats ' profile_old_b ' to be DMQCed: ',floats_tobedone_str,'%'],...
        'FitBoxToText','on','verticalalignment', 'bottom','HorizontalAlignment', 'right','FontWeight','bold')

    % save figure
    out_name = [output_plots_dir '/' sprintf('%02d',i_fig) '_' project_name '_' char(i_param_str) '_DMQC_status_nb_floats_by_country_' working_date];
    disp(['saving ' out_name])

%     export_fig([out_name '.png'])
    print('-dpng ', '-r100',[out_name '.png'])
    if print_svg == 1
        saveas(gcf,[out_name '.svg'])
    end

% end
% %%
% %for i=1:n_param
% for i=1:1
%     i_param=list_of_parameters_to_treat(i);
%     fprintf('Making Plots for %s \n',i_param)
%     
%     close all

    %%%%%%%%%%%%%% Profile DMQC status by country %%%%%%%%%%%%%%
    disp('Profile DMQC status by country')
    i_fig=2;
    figure(i_fig)
    % bigger figure
    set(gcf, 'Position', [200, 200, 1000, 600])
    % figure name
    set(gcf,'Name','Profile DMQC status by country')

    stackData = cat(3,[nb_profiles_DMQCed_per_country.(i_param) nb_profiles_xage_DMQCed_per_country.(i_param)], ...
                      [(nb_profiles_per_country.(i_param)-nb_profiles_DMQCed_per_country.(i_param)) ...
                       (nb_profiles_xage_per_country.(i_param)-nb_profiles_xage_DMQCed_per_country.(i_param))]);
    hndl = plotBarStackGroups(stackData, nb_per_country_x);


    % FIGURE FORMAT
    hold on
    hp = plot(1,0,'w'); % for comment in legend
    % bar colors 
    set(hndl(1,1),'facecolor',bars_colors(3,:))
    set(hndl(1,2),'facecolor',bars_colors(1,:))
    set(hndl(2,1),'facecolor',bars_colors(4,:))
    set(hndl(2,2),'facecolor',bars_colors(2,:))
    % xlabels
    set(gca,'xtick',1:n_countries,'xticklabel', nb_per_country_x)
    if n_countries > 15
        set(gca,'XTickLabelRotation',90)
    end
    % title with update date
    title(['Profile DMQC status for ' char(i_param_str) ' by country (updated ',update_date_str,')'], 'Interpreter', 'none')
    ylabel('Number of profiles')
    % legend with total number
    legend([hndl(1,2), hndl(1,1), hndl(2,2), hndl(2,1), hp], ...
        {['Number of profiles (Total: ',num2str(sum(nb_profiles_per_country.(i_param),'omitnan')),')'],...
         ['Number of D-profiles (Total: ',num2str(sum(nb_profiles_DMQCed_per_country.(i_param),'omitnan')),')'], ...
         ['Number of profiles ' profile_old_b ' (Total: ',num2str(sum(nb_profiles_xage_per_country.(i_param),'omitnan')),')'],...
         ['Number of D-profiles ' profile_old_b 'r (Total: ',num2str(sum(nb_profiles_xage_DMQCed_per_country.(i_param),'omitnan')),')'], ...
        prof_included})
    % background color
    set(gcf,'color','w');
    % grid in y axis
    ax = gca;
    ax.YGrid = 'on';
    
    ymax=max(nb_profiles_per_country.(i_param));
    set(gca,'YLim',[0 ymax+ymax/3]);

    % annotation: percentage of observations > 1 year with no DMQC done
    obs_tobedone = sum(nb_profiles_xage_per_country.(i_param)-nb_profiles_xage_DMQCed_per_country.(i_param),'omitnan')/sum(nb_profiles_xage_per_country.(i_param),'omitnan')*100;
    H = figure(i_fig);
    set(H,'units','pix')
    annotation('textbox', [0.8, 0.01, .1, .1], 'string', ['Profiles ' profile_old_b ' to be DMQCed: ',num2str(round(obs_tobedone)),'%'],...
        'FitBoxToText','on','verticalalignment', 'bottom','HorizontalAlignment', 'right','FontWeight','bold')

    % save figure
    out_name = [output_plots_dir '/' sprintf('%02d',i_fig) '_' project_name '_' char(i_param_str) '_DMQC_status_nb_profiles_by_country_' working_date];
    disp(['saving ' out_name])
%     export_fig([out_name '.png'])
    print('-dpng ', '-r100',[out_name '.png'])
    if print_svg == 1
        saveas(gcf,[out_name '.svg'])
    end
    
end
diary off
%% 
diary on
for i=1:n_param
    close all
 
    i_param=list_of_parameters_to_treat(i);
    
    if i_param == "PRES"
        % profile QC for PRES no yet available in index file
        continue
    end
    
    fprintf('Making Plots for %s \n',i_param)
    
    close all
     %%%%%%%%%%%%%% Data Quality %%%%%%%%%%%%%%
    disp('Data Quality')
    i_fig=3;
    figure(i_fig)
    % bigger figure
    set(gcf, 'Position', [200, 200, 1000, 600])
    % figure name
    set(gcf,'Name','Data Quality')

    colormap(parula)
    hold on

    for ix=1:7
        if i_group_AB_profQC == 1 && ix==1
            nb_prof_AB=nb_profiles_per_qc.(i_param)(1)+nb_profiles_per_qc.(i_param)(2);
            nb_prof_DMQCed_AB=nb_profiles_DMQCed_per_qc.(i_param)(1)+nb_profiles_DMQCed_per_qc.(i_param)(2);
            
            bT=bar(2,nb_prof_AB,0.3,'FaceColor',bars_colors(4+ix,:));
            bD=bar(2+0.3,nb_prof_DMQCed_AB,0.3,'FaceColor',bars_colors(13+ix,:));  
        else
            if i_group_AB_profQC == 0 || ix > 2
                bT=bar(ix,nb_profiles_per_qc.(i_param)(ix),0.3,'FaceColor',bars_colors(4+ix,:));
                bD=bar(ix+0.3,nb_profiles_DMQCed_per_qc.(i_param)(ix),0.3,'FaceColor',bars_colors(13+ix,:));
            end
        end
            
    end
    

    % FIGURE FORMAT
    hold on
    plot(1,0,'w'); % for comment in legend
    % xlabels
    if i_group_AB_profQC == 1
        set(gca,'xtick',1:7,'xticklabel', [''; 'A+B'; nb_per_qc_x(3:end-1); 'No QC'])
    else
        set(gca,'xtick',1:7,'xticklabel', [nb_per_qc_x(1:end-1); 'No QC'])
    end
    
    % title with update date
    title([ char(i_param) ' profile QC (updated ',update_date_str,')'], 'Interpreter', 'none')
    ylabel('Number of profiles')

    % figure labels
    if i_group_AB_profQC == 1
        ymax=max(max(nb_profiles_per_qc.(i_param)(1)+nb_profiles_per_qc.(i_param)(2),nb_profiles_per_qc.(i_param)(3:end)));
    else
        ymax=max(nb_profiles_per_qc.(i_param));
    end
    set(gca,'YLim',[0 ymax+ymax/4]);
    for ix =1:7
        
        if i_group_AB_profQC == 1 && ix==1
            text(2 +0.17 , nb_prof_AB + ymax/20, ...
                [num2str(nb_prof_AB) newline ...
                 num2str(round(100*nb_prof_AB/ ...
                           sum(nb_profiles_per_qc.(i_param)),1)) '%'], ...
                'HorizontalAlignment', 'right');

            text(2 +0.17, nb_prof_DMQCed_AB + ymax/20, ...
                ['D-' num2str(nb_prof_DMQCed_AB) newline ...
                      num2str(round(100*nb_prof_DMQCed_AB/ ...
                                sum(nb_profiles_DMQCed_per_qc.(i_param)),1)) '%'], ...
                'HorizontalAlignment', 'left');
        else
            if i_group_AB_profQC == 0 || ix > 2
                text(ix+0.17 , nb_profiles_per_qc.(i_param)(ix) + ymax/20, ...
                    [num2str(nb_profiles_per_qc.(i_param)(ix)) newline ...
                     num2str(round(100*nb_profiles_per_qc.(i_param)(ix)/ ...
                               sum(nb_profiles_per_qc.(i_param)),1)) '%'], ...
                    'HorizontalAlignment', 'right');

                text(ix+0.17, nb_profiles_DMQCed_per_qc.(i_param)(ix) + ymax/20, ...
                    ['D-' num2str(nb_profiles_DMQCed_per_qc.(i_param)(ix)) newline ...
                          num2str(round(100*nb_profiles_DMQCed_per_qc.(i_param)(ix)/ ...
                                    sum(nb_profiles_DMQCed_per_qc.(i_param)),1)) '%'], ...
                    'HorizontalAlignment', 'left');
            end
        end
        
        
    end
    % background color
    set(gcf,'color','w');
    % grid in y axis
    ax = gca;
    ax.YGrid = 'on';
    box on
    
    % annotation
    H = figure(i_fig);
    set(H,'units','pix')
    annotation('textbox', [0.28, 0.83, .1, .1], 'string', ['lighter color bars / xxxx : nb of profiles (Raw, Adjusted or Delayed mode) per profile QC' newline ...
        'darker color bars / D-xxxx : nb of D-profiles (Delayed mode only) per profile QC'],...
        'FitBoxToText','on','verticalalignment', 'bottom','HorizontalAlignment', 'left')

    % save figure
    out_name = [output_plots_dir '/' sprintf('%02d',i_fig) '_' project_name '_' char(i_param) '_profile_QC_' working_date];
    disp(['saving ' out_name])
%     export_fig([out_name '.png'])
    print('-dpng ', '-r100',[out_name '.png'])
    if print_svg == 1
        saveas(gcf,[out_name '.svg'])
    end
% 
% end
% 
% 
% %%
% for i=1:n_param
%     i_param=list_of_parameters_to_treat(i);
%     fprintf('Making Plots for %s \n',i_param)
%     
%     close all
     %%%%%%%%%%%%%% Data Quality %%%%%%%%%%%%%%
    disp('Data Quality Evolution per Launch date')
    i_fig=4;
    figure(i_fig)
    
    % bigger figure
    set(gcf, 'Position', [200, 200, 1000, 600])
    % figure name
    set(gcf,'Name','Data Quality evolution')
    % background color
    set(gcf,'color','w');
    
    % title with update date
    title([ char(i_param) ' profile QC evolution (updated ',update_date_str,')'], 'Interpreter', 'none')
    hold on
    
    %tiledlayout(2,1)
    % Top plot
    ax1 = subplot(211);
    title(ax1,[ char(i_param) ' profile QC evolution (updated ',update_date_str,')'], 'Interpreter', 'none')
    hold on
    if i_group_AB_profQC == 1
        plot(ax1,1:n_launch_years,percent_prof_QC_A_per_launch_year.(i_param)+percent_prof_QC_B_per_launch_year.(i_param),'color',bars_colors(5,:),'LineWidth',2)
    else
        plot(ax1,1:n_launch_years,percent_prof_QC_A_per_launch_year.(i_param),'color',bars_colors(5,:),'LineWidth',2)
        plot(ax1,1:n_launch_years,percent_prof_QC_B_per_launch_year.(i_param),'color',bars_colors(6,:),'LineWidth',2)
    end
    plot(ax1,1:n_launch_years,percent_prof_QC_C_per_launch_year.(i_param),'color',bars_colors(7,:),'LineWidth',2)
    plot(ax1,1:n_launch_years,percent_prof_QC_D_per_launch_year.(i_param),'color',bars_colors(8,:),'LineWidth',2)
    plot(ax1,1:n_launch_years,percent_prof_QC_E_per_launch_year.(i_param),'color',bars_colors(9,:),'LineWidth',2)
    plot(ax1,1:n_launch_years,percent_prof_QC_F_per_launch_year.(i_param),'color',bars_colors(10,:),'LineWidth',2)
    
    if i_group_AB_profQC == 1
        lgnd=legend(ax1,'QC A+B','QC C' ,'QC D' ,'QC E' ,'QC F' ,...
                    'location','east');
    else
        lgnd=legend(ax1,'QC A','QC B','QC C' ,'QC D' ,'QC E' ,'QC F' ,...
                    'location','east');
    end
    set(lgnd,'color','none');
    
    ylabel(ax1,'Percent of profiles [%]')
    %xlabel(ax1,'Launch year')
    set(ax1,'xtick',1:n_launch_years,'xticklabel', nb_per_launch_year_x)
    ax1.YGrid = 'on';
    ax1.XTickLabelRotation = 45;
    set(ax1, 'xlim',[1,n_launch_years+3])
    box on
    
    % Bottom plot
    ax2 = subplot(212);
    hold on
    plot(ax2,1:n_launch_years,nb_prof_per_launch_year.(i_param),'color','black','LineWidth',2)
    ylabel(ax2,'Number of profiles')
    xlabel(ax2,'Launch year')
    set(ax2,'xtick',1:n_launch_years,'xticklabel', nb_per_launch_year_x)
    ax2.YGrid = 'on';
    ax2.XTickLabelRotation = 45;
    set(ax2, 'xlim',[1,n_launch_years+3])
    box on
    
    % save figure
    out_name = [output_plots_dir '/' sprintf('%02d',i_fig) '_' project_name '_' char(i_param) '_profile_QC_evolution_' working_date];
    disp(['saving ' out_name])
%     export_fig([out_name '.png'])

    print('-dpng ', '-r100',[out_name '.png'])
    if print_svg == 1
        saveas(gcf,[out_name '.svg'])
    end
    
    
% end
% 
% %%
% for i=1:n_param
%     i_param=list_of_parameters_to_treat(i);
%     fprintf('Making Plots for %s \n',i_param)
%     
%     close all
     %%%%%%%%%%%%%% Data Quality %%%%%%%%%%%%%%
    disp('Data Quality Evolution for D-profiles per Launch date')
    i_fig=5;
    figure(i_fig)
    
    % bigger figure
    set(gcf, 'Position', [200, 200, 1000, 600])
    % figure name
    set(gcf,'Name','Data Quality evolution')
    % background color
    set(gcf,'color','w');
    
    % title with update date
    title([ char(i_param) ' D-profile QC evolution (updated ',update_date_str,')'], 'Interpreter', 'none')
    hold on
    
    %tiledlayout(2,1)
    % Top plot
    ax1 = subplot(211);
    title(ax1,[ char(i_param) ' D-profile QC evolution (updated ',update_date_str,')'], 'Interpreter', 'none')
    hold on
    if i_group_AB_profQC == 1
        plot(ax1,1:n_launch_years,percent_prof_DMQCed_QC_A_per_launch_year.(i_param)+percent_prof_DMQCed_QC_B_per_launch_year.(i_param),'color',bars_colors(5,:),'LineWidth',2)
    else
        plot(ax1,1:n_launch_years,percent_prof_DMQCed_QC_A_per_launch_year.(i_param),'color',bars_colors(5,:),'LineWidth',2)
        plot(ax1,1:n_launch_years,percent_prof_DMQCed_QC_B_per_launch_year.(i_param),'color',bars_colors(6,:),'LineWidth',2)
    end
    plot(ax1,1:n_launch_years,percent_prof_DMQCed_QC_C_per_launch_year.(i_param),'color',bars_colors(7,:),'LineWidth',2)
    plot(ax1,1:n_launch_years,percent_prof_DMQCed_QC_D_per_launch_year.(i_param),'color',bars_colors(8,:),'LineWidth',2)
    plot(ax1,1:n_launch_years,percent_prof_DMQCed_QC_E_per_launch_year.(i_param),'color',bars_colors(9,:),'LineWidth',2)
    plot(ax1,1:n_launch_years,percent_prof_DMQCed_QC_F_per_launch_year.(i_param),'color',bars_colors(10,:),'LineWidth',2)
    
    if i_group_AB_profQC == 1
        lgnd=legend(ax1,'QC A+B','QC C' ,'QC D' ,'QC E' ,'QC F' ,...
                    'location','east');
    else
        lgnd=legend(ax1,'QC A','QC B','QC C' ,'QC D' ,'QC E' ,'QC F' ,...
                    'location','east');
    end  
    set(lgnd,'color','none');
    
    ylabel(ax1,'Percent of profiles [%]')
    %xlabel(ax1,'Launch year')
    set(ax1,'xtick',1:n_launch_years,'xticklabel', nb_per_launch_year_x)
    ax1.XTickLabelRotation = 45;
    set(ax1, 'xlim',[1,n_launch_years+3])
    ax1.YGrid = 'on';
    box on
    
    % Bottom plot
    ax2 = subplot(212);
    hold on
    plot(ax2,1:n_launch_years,nb_prof_DMQCed_per_launch_year.(i_param),'color','black','LineWidth',2)
    ylabel(ax2,'Number of profiles')
    xlabel(ax2,'Launch year')
    set(ax2,'xtick',1:n_launch_years,'xticklabel', nb_per_launch_year_x)
    ax2.XTickLabelRotation = 45;
    set(ax2, 'xlim',[1,n_launch_years+3])
    ax2.YGrid = 'on';
    box on
    
    % save figure
    out_name = [output_plots_dir '/' sprintf('%02d',i_fig) '_' project_name '_' char(i_param) '_Dprofile_QC_evolution_' working_date];
    disp(['saving ' out_name])
%     export_fig([out_name '.png'])
    print('-dpng ', '-r100',[out_name '.png'])
    if print_svg == 1
        saveas(gcf,[out_name '.svg'])
    end
    
    
end
diary off
%% 
diary on
for i=1:n_param
    close all
    
    i_param=list_of_parameters_to_treat(i);
    fprintf('Making Plots for %s \n',i_param)
    
    close all

     %%%%%%%%%%%%%% exclusion list  %%%%%%%%%%%%%%
    disp('exclusion list status')
    i_fig=6;
    figure(i_fig)
    % bigger figure
    set(gcf, 'Position', [200, 200, 1000, 600])
    % figure name
    set(gcf,'Name','exclusion list')



    bar(1, nb_float_tot.(i_param),              'FaceColor', bars_colors(1,:))
    hold on
    bar(2, nb_float_xage_tot.(i_param),         'FaceColor', bars_colors(2,:))
    bar(3, nb_float_DMQCed_tot.(i_param),       'FaceColor', bars_colors(3,:))
    bar(4, nb_floats_xage_DMQCed_tot.(i_param), 'FaceColor', bars_colors(4,:))
    bar(5, nb_wmo_ope_in_exclusionlist_with_qc_3_or_4.(i_param),       'FaceColor', bars_colors(11,:))
    bar(6, nb_wmo_inactiveR_in_exclusionlist_with_qc_3_or_4.(i_param), 'FaceColor', bars_colors(12,:))
    bar(7, nb_wmo_inactiveD_in_exclusionlist_with_qc_3_or_4.(i_param),  'FaceColor', bars_colors(13,:))

    % FIGURE FORMAT
    % xlabels
    set(gca,'XTick',[])
    % title with update date
    title(['Floats status and exclusion list for ' char(i_param) ' (updated ',update_date_str,')'], 'Interpreter', 'none')
    ylabel('Number of floats')
    % legend with total number
    legend(['All floats ('          num2str(nb_float_tot.(i_param)) ')'],...
           ['Floats ' profile_old_b ' ('       num2str(nb_float_xage_tot.(i_param)) ')'],...
           ['Floats DMQC at least once ('              num2str(nb_float_DMQCed_tot.(i_param)) ')'],...
           ['Floats ' profile_old_b ' and DMQC at least once ('    num2str(nb_floats_xage_DMQCed_tot.(i_param)) ')'],...
           ['Active Floats in exclusion list ('     num2str(nb_wmo_ope_in_exclusionlist_with_qc_3_or_4.(i_param)) ...
              '/' num2str(nb_operational_floats.(i_param)) ' ->' ...
              num2str(round(100*nb_wmo_ope_in_exclusionlist_with_qc_3_or_4.(i_param)/nb_operational_floats.(i_param),1)) '% of active floats)'],...
           ['Inactive Floats with R or A profiles in exclusion list (' num2str(nb_wmo_inactiveR_in_exclusionlist_with_qc_3_or_4.(i_param)) ')'], ...
           ['Inactive Floats only with D profiles in exclusion list (' num2str(nb_wmo_inactiveD_in_exclusionlist_with_qc_3_or_4.(i_param)) ')'], ...
        'Location','southoutside')

    % background color
    set(gcf,'color','w');
    % grid in y axis
    ax = gca;
    ax.YGrid = 'on';
    
    % figure labels
    ymax=max(nb_float_tot.(i_param));
    set(gca,'YLim',[0 ymax+ymax/5]);
    text(1, nb_float_tot.(i_param) + ymax/20, num2str(nb_float_tot.(i_param)), ...
        'HorizontalAlignment', 'center');
    text(2, nb_float_xage_tot.(i_param) + ymax/20, num2str(nb_float_xage_tot.(i_param)), ...
        'HorizontalAlignment', 'center');
    text(3, nb_float_DMQCed_tot.(i_param) + ymax/20, num2str(nb_float_DMQCed_tot.(i_param)), ...
        'HorizontalAlignment', 'center');
    text(4, nb_floats_xage_DMQCed_tot.(i_param) + ymax/20, num2str(nb_floats_xage_DMQCed_tot.(i_param)), ...
        'HorizontalAlignment', 'center');
    text(5, nb_wmo_ope_in_exclusionlist_with_qc_3_or_4.(i_param) + ymax/20, num2str(nb_wmo_ope_in_exclusionlist_with_qc_3_or_4.(i_param)), ...
        'HorizontalAlignment', 'center');
    text(6, nb_wmo_inactiveR_in_exclusionlist_with_qc_3_or_4.(i_param) + ymax/20, num2str(nb_wmo_inactiveR_in_exclusionlist_with_qc_3_or_4.(i_param)), ...
        'HorizontalAlignment', 'center');
    text(7, nb_wmo_inactiveD_in_exclusionlist_with_qc_3_or_4.(i_param) + ymax/20, num2str(nb_wmo_inactiveD_in_exclusionlist_with_qc_3_or_4.(i_param)), ...
        'HorizontalAlignment', 'center');


    % save figure
    out_name = [output_plots_dir '/' sprintf('%02d',i_fig) '_' project_name '_' char(i_param) '_DMQC_status_and_exclusion_list_' working_date];
    disp(['saving ' out_name])
%     export_fig([out_name '.png'])
    print('-dpng ', '-r100',[out_name '.png'])
    if print_svg == 1
        saveas(gcf,[out_name '.svg'])
    end

end
diary off
%% 
diary on
for i=1:n_param
    close all
    
    i_param=list_of_parameters_to_treat(i);
    
    if ismember("PSAL",list_of_parameters_to_treat)
        if (i_param == "TEMP" || i_param == "PRES") 
            % no need to output plot as mode is the same for TEMP, PRES,
            % PSAL
            continue
        else
            if i_param == "PSAL"
                i_param_str="CTD";
            else
                i_param_str=i_param;
            end
        end
    else
        i_param_str=i_param;
    end
    
    fprintf('Making Plots for %s \n',i_param_str)
    
    close all
    %%%%%%%%%%%%%% DMQC CTD status per profile year %%%%%%%%%%%%%%
    disp('DMQC status per profile year')
    i_fig=7;
    figure(i_fig)

    % bigger figure
    set(gcf, 'Position', [200, 200, 1000, 600])
    % figure name
    set(gcf,'Name','DMQC status per profile year')

    if (nb_per_year_x(1) == string(1980))
        nb_per_year_x_cr=nb_per_year_x(2:end);
        nb_profiles_DMQCed_per_year_cr=nb_profiles_DMQCed_per_year.(i_param)(2:end);
        nb_profiles_per_year_cr=nb_profiles_per_year.(i_param)(2:end);
    else
        nb_per_year_x_cr=nb_per_year_x;
        nb_profiles_DMQCed_per_year_cr=nb_profiles_DMQCed_per_year.(i_param);
        nb_profiles_per_year_cr=nb_profiles_per_year.(i_param);
    end
    
    hndl = bar(1:length(nb_per_year_x_cr), ...
           [nb_profiles_per_year_cr nb_profiles_DMQCed_per_year_cr], 1);

    % FIGURE FORMAT
    hold on
    %plot(1,0,'w'); % for comment in legend
    % bars colors
    set(hndl(1),'facecolor',bars_colors(1,:))
    set(hndl(2),'facecolor',bars_colors(3,:))
    % xlabels
    %set(gca,'xtick',str2double(nb_per_year_x))
    set(gca,'xtick',1:length(nb_per_year_x_cr),'xticklabel', nb_per_year_x_cr)
    xtickangle(90)
    % title with update date
    title(['DMQC status for ' char(i_param_str) ' per profile year (updated ',update_date_str,')'], 'Interpreter', 'none')
    ylabel('Number of profiles')
    % legend with total number
    legend(['Number of profiles (Total: ',num2str(sum(nb_profiles_per_year_cr,'omitnan')),')'],...
           ['Number of D-profiles (Total: ',num2str(sum(nb_profiles_DMQCed_per_year_cr,'omitnan')),')'], ...
        prof_included,'Location','southoutside')
    % background color
    set(gcf,'color','w');
    % grid in y axis
    ax = gca;
    ax.YGrid = 'on';

    % figure labels
    ymax=max(nb_profiles_per_year_cr);
    set(gca,'YLim',[0 ymax+ymax/10]);
    for iyear = 1: length(nb_per_year_x_cr)
        if nb_profiles_per_year_cr(iyear) > 0
            percent_str = [num2str(round(nb_profiles_DMQCed_per_year_cr(iyear)/nb_profiles_per_year_cr(iyear)*100,1)) ' %'];
            h=text(iyear + 0.15, nb_profiles_DMQCed_per_year_cr(iyear) + ymax/30 ,percent_str);
            set(h,'Rotation',90);
        end
    end

    % save figure
    out_name = [output_plots_dir '/' sprintf('%02d',i_fig) '_' project_name '_' char(i_param_str) '_prof_DMQCstatus_byyear_' working_date];
    disp(['saving ' out_name])
%     export_fig([out_name '.png'])
    print('-dpng ', '-r100',[out_name '.png'])
    if print_svg == 1
        saveas(gcf,[out_name '.svg'])
    end

    
    
end
diary off
%% 
diary on
%for i=1:3
for i=1:n_param
    close all
    
    i_param=list_of_parameters_to_treat(i);
    
    if ismember("PSAL",list_of_parameters_to_treat)
        if (i_param == "PRES") 
            % no need to output plot as mode fro PRES is the same as mode 
            % for PSAL eand TEMP and profile_qc unavailable in index.
            continue
        else
            i_param_str=i_param;
        end
    else
        i_param_str=i_param;
    end
    
    fprintf('Making Plots for %s \n',i_param_str)
    
    close all
    %%%%%%%%%%%%%% DMQC status per profile year %%%%%%%%%%%%%%
    disp('DMQC and QC-F status per profile year')
    i_fig=8;
    figure(i_fig)

    % bigger figure
    set(gcf, 'Position', [200, 200, 1000, 600])
    % figure name
    set(gcf,'Name','DMQC and F status per profile year')

    if (nb_per_year_x(1) == string(1980))
        nb_per_year_x_cr=nb_per_year_x(2:end);
        nb_profiles_DMQCed_per_year_cr=nb_profiles_DMQCed_per_year.(i_param)(2:end);
        nb_profiles_per_year_cr=nb_profiles_per_year.(i_param)(2:end);
        nb_profiles_QCF_per_year_cr=nb_profiles_QCF_per_year.(i_param)(2:end);
        nb_profiles_DMQCed_QCF_per_year_cr=nb_profiles_DMQCed_QCF_per_year.(i_param)(2:end);
    else
        nb_per_year_x_cr=nb_per_year_x;
        nb_profiles_DMQCed_per_year_cr=nb_profiles_DMQCed_per_year.(i_param);
        nb_profiles_per_year_cr=nb_profiles_per_year.(i_param);
        nb_profiles_QCF_per_year_cr=nb_profiles_QCF_per_year.(i_param);
        nb_profiles_DMQCed_QCF_per_year_cr=nb_profiles_DMQCed_QCF_per_year.(i_param);
    end
    
    stackData = cat(3,[(nb_profiles_per_year_cr-nb_profiles_QCF_per_year_cr) ...
                       (nb_profiles_DMQCed_per_year_cr-nb_profiles_DMQCed_QCF_per_year_cr)], ...
                      [ nb_profiles_QCF_per_year_cr ...
                        nb_profiles_DMQCed_QCF_per_year_cr]);
                   
    hndl = plotBarStackGroups(stackData, 1:length(nb_per_year_x_cr));
    

    % FIGURE FORMAT
    hold on
    %plot(1,0,'w'); % for comment in legend
    % bars colors
    set(hndl(1,1),'facecolor',bars_colors(4,:))
    set(hndl(1,2),'facecolor',bars_colors(21,:))
    set(hndl(2,1),'facecolor',bars_colors(3,:))
    set(hndl(2,2),'facecolor',bars_colors(22,:))
    % xlabels
    set(gca,'xtick',1:length(nb_per_year_x_cr),'xticklabel', nb_per_year_x_cr)
    xtickangle(90)
    % title with update date
    title(['DMQC and profile-F status for ' char(i_param_str) ' per profile year (updated ',update_date_str,')'], 'Interpreter', 'none')
    ylabel('Number of profiles')
    % legend with total number
    legend(['Number of profiles with prof_QC<>F(Total: ',num2str(sum(nb_profiles_per_year_cr,'omitnan')),')'],...
           ['Number of profiles with prof_QC=F (Total: ',num2str(sum(nb_profiles_QCF_per_year_cr,'omitnan')),')'],...
           ['Number of D-profiles with prof_QC<>F (Total: ',num2str(sum(nb_profiles_DMQCed_per_year_cr,'omitnan')),')'], ...
           ['Number of D-profiles with prof_QC=F (Total: ',num2str(sum(nb_profiles_DMQCed_QCF_per_year_cr,'omitnan')),')'], ...
        prof_included,'Location','southoutside', 'Interpreter', 'none')
    % background color
    set(gcf,'color','w');
    % grid in y axis
    ax = gca;
    ax.YGrid = 'on';

    % figure labels
    ymax=max(nb_profiles_per_year_cr);
    set(gca,'YLim',[0 ymax+ymax/10]);
%     for iyear = 1: length(nb_per_year_x_cr)
%         if nb_profiles_per_year_cr(iyear) > 0
%             percent_str = [num2str(round(nb_profiles_DMQCed_per_year_cr(iyear)/nb_profiles_per_year_cr(iyear)*100,1)) ' %'];
%             h=text(iyear + 0.15, nb_profiles_DMQCed_per_year_cr(iyear) + ymax/30 ,percent_str);
%             set(h,'Rotation',90);
%         end
%     end

    % save figure
    out_name = [output_plots_dir '/' sprintf('%02d',i_fig) '_' project_name '_' char(i_param_str) '_prof_DMQC-and-F_status_byyear_' working_date];
    disp(['saving ' out_name])
%     export_fig([out_name '.png'])
    print('-dpng ', '-r100',[out_name '.png'])
    if print_svg == 1
        saveas(gcf,[out_name '.svg'])
    end

    
    
end
diary off
%% 
diary on
for i=1:n_param
    i_param=list_of_parameters_to_treat(i);
    fprintf('Making Plots for %s \n',i_param)
    
    close all

     %%%%%%%%%%%%%% DMQC status by profile age histogram %%%%%%%%%%%%%%
    disp('DMQC status by profile age histogram')
    i_fig=9;
    figure(i_fig)

    % bigger figure
    set(gcf, 'Position', [200, 200, 1000, 600])
    % figure name
    set(gcf,'Name','R/A-Profiles age')

    hndl = histogram(IndexData.profile_age(i_R_or_A_profile.(i_param))/365);

    % FIGURE FORMAT
    hold on
    plot(ones(1,hndl.NumBins),hndl.Values,'--k','LineWidth',2.5); % for comment in legend
    % xlabels
    xlabel('Age of profiles [years]')
    % title with update date
    title(['Age of ' char(i_param_str) ' profiles with no DMQC (updated ',update_date_str,')'], 'Interpreter', 'none')
    ylabel('Number of profiles')
    % background color
    set(gcf,'color','w');
    % grid in y axis
    ax = gca;
    ax.YGrid = 'on';

    % save figure
    out_name = [output_plots_dir '/' sprintf('%02d',i_fig) '_' project_name '_' char(i_param_str) '_prof_DMQCstatus_agehist_' working_date];
    disp(['saving ' out_name])
    %     export_fig([out_name '.png'])
    print('-dpng ', '-r100',[out_name '.png'])
    if print_svg == 1
        saveas(gcf,[out_name '.svg'])
    end
    
end
tEnd = toc(tStart);

diary off
%% 
diary on
% Additionnal graphs from previous get_DMQC_adjustment.m routine.
% Thanks to the update of the index file with psal adjustment information
% We can switch the making of the graphs to this routine. 
% There are 2 new kinds of graphs:
%  - R/A/D profiles per parameter
%  - a series of graphs by wmo (R/A/D, QC, PSAL_adj) (and by parameter).

disp('Number of R/A/D profiles by parameter')
i_fig=10;
close all
figure(i_fig)

% bigger figure
set(gcf, 'Position', [200, 200, 1000, 600])
% figure name
set(gcf,'Name','R/A/D profiles nb by parameter')

set(gca, 'units', 'normalized'); %Just making sure it's normalized
Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
                                 %[Left Bottom Right Top] spacing
NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3) 1-Tight(2)-Tight(4)-0.1]; %New plot position [X Y W H]
set(gca, 'Position', NewPos);

% To avoid the 3 same information for CTD:
if list_of_parameters_to_treat(1)=="TEMP" && ...
   list_of_parameters_to_treat(2)=="PRES" && ...
   list_of_parameters_to_treat(3)=="PSAL"

    if i_bgc == 1
        loc_list_of_parameters_to_treat=list_of_parameters_to_treat(4:end);
        loc_n_param=size(loc_list_of_parameters_to_treat,1);
        ctd_first=0;
    else
        loc_list_of_parameters_to_treat=list_of_parameters_to_treat(3:end);
        loc_n_param=size(loc_list_of_parameters_to_treat,1);
        ctd_first=1;
    end
else
    loc_list_of_parameters_to_treat=list_of_parameters_to_treat;
    loc_n_param=n_param;
    ctd_first=0;
end

stackData = NaN(loc_n_param,3);

for i=1:loc_n_param
    i_param=loc_list_of_parameters_to_treat(i);

    stackData(i,:) = cat(3,[sum(string(IndexData.param.(i_param).mode) == "R")'...
                            sum(string(IndexData.param.(i_param).mode) == "A")'...
                            sum(string(IndexData.param.(i_param).mode) == "D")']);
end

hndl = plotBarStackGroups(stackData, loc_list_of_parameters_to_treat');

% FIGURE FORMAT
hold on
hp = plot(1,0,'w'); % for comment in legend
% bar colors 
set(hndl(1),'facecolor',bars_colors(15,:))
set(hndl(2),'facecolor',bars_colors(2,:))
set(hndl(3),'facecolor',bars_colors(3,:))
% xlabels
if ctd_first == 1
    loc_list_of_parameters_to_treat(1)="CTD";
end
loc_list_of_parameters_to_treat_wo__=strrep(loc_list_of_parameters_to_treat,"_"," ");
set(gca,'xtick',1:loc_n_param,'xticklabel', loc_list_of_parameters_to_treat_wo__','XTickLabelRotation',45)


% title with update date
title(['Number of R/A/D - profiles by variable (updated ',update_date_str,')' newline project_name ], 'Interpreter', 'none')
ylabel('Number of profiles')

% legend with total number
legend([hndl(1), hndl(2), hndl(3), hp], ...
    {['R-profiles (Total: ',num2str(sum(stackData(:,1))),')'],...
     ['A-profiles (Total: ',num2str(sum(stackData(:,2))),')'],...
     ['D-profiles (Total: ',num2str(sum(stackData(:,3))),')'],...
    prof_included},'Location','southoutside')

% background color
set(gcf,'color','w');
% grid in y axis
ax = gca;
ax.YGrid = 'on';

 % save figure
out_name = [output_plots_dir '/' sprintf('%02d',i_fig) '_' project_name '_prof_RAD_mode_per_param_' working_date];
disp(['saving ' out_name])
%     export_fig([out_name '.png'])
print('-dpng ', '-r100',[out_name '.png'])
if print_svg == 1
    saveas(gcf,[out_name '.svg'])
end

diary off
%% 
diary on
if output_graphs_per_float ==1
    
    disp('R/A/D status by float and by cycle number')
    i_fig=11;
    for i=1:n_param
        close all
        i_param=list_of_parameters_to_treat(i);
        
        
        if ismember("PSAL",list_of_parameters_to_treat)
            if (i_param == "TEMP" || i_param == "PRES") 
                % no need to output plot as mode is the same for TEMP, PRES,
                % PSAL
                continue
            else
                if i_param == "PSAL"
                    i_param_str="CTD";
                else
                    i_param_str=i_param;
                end
            end
        else
            i_param_str=i_param;
        end
        
        fprintf('Making Plots for %s \n',i_param_str)
                
        subset_prof_WMO=IndexData.profile_WMO(IndexData.param.(i_param).presence==1);
        subset_prof_cycle=IndexData.cycle(IndexData.param.(i_param).presence==1);
        subset_prof_mode=IndexData.param.(i_param).mode(IndexData.param.(i_param).presence==1);
        
        list_wmo_for_graph=unique(subset_prof_WMO);
        n_wmo_param=size(list_wmo_for_graph,1);
        n_graph=ceil(n_wmo_param/n_max_float_per_graph);
        
        
        for i_graph=1:n_graph
            
            i_wmo_min=1+(i_graph-1)*n_max_float_per_graph;
            i_wmo_max=min(n_wmo_param,i_graph*n_max_float_per_graph);
            list_wmo_for_i_graph=list_wmo_for_graph(i_wmo_min:i_wmo_max);
            n_wmo_for_i_graph=size(list_wmo_for_i_graph,1);
            
            prof_subset=ismember(subset_prof_WMO,list_wmo_for_i_graph);
            
            xx=double(subset_prof_cycle((prof_subset==1)));
            yy=subset_prof_WMO((prof_subset==1));
            zz=string(subset_prof_mode((prof_subset==1)));
            
            i_R=(zz=="R");
            i_A=(zz=="A");
            i_D=(zz=="D");
            
            close all
            figure(i_fig)
            

            % bigger figure
            set(gcf, 'Position', [200, 200, 1200, 600])
            % figure name
            set(gcf,'Name','R/A/D status by float and by cycle number')
            hold on
            
            % title with update date
            title(['R/A/D status by float and by cycle for ' char(i_param_str) ' (updated ',update_date_str,')' ...
                  newline project_name ' - batch ' num2str(i_graph) '/' num2str(n_graph)], 'Interpreter', 'none')
            xlabel('Cycle Number')
            
            xx_iR=xx(i_R);
            [~,yy_iR]=ismember(yy(i_R),list_wmo_for_i_graph);
            
            xx_iA=xx(i_A);
            [~,yy_iA]=ismember(yy(i_A),list_wmo_for_i_graph);
            
            xx_iD=xx(i_D);
            [~,yy_iD]=ismember(yy(i_D),list_wmo_for_i_graph);
            
            scatter(xx_iR,yy_iR,10,'o','filled','MarkerFaceColor',bars_colors(15,:))
            scatter(xx_iA,yy_iA,10,'o','filled','MarkerFaceColor',bars_colors(2,:))
            scatter(xx_iD,yy_iD,10,'o','filled','MarkerFaceColor',bars_colors(3,:))
            
            
            set(gca,'ytick',1:n_wmo_for_i_graph,'yticklabel', list_wmo_for_i_graph')
            % legend
            lh = legend(['R-profiles (',num2str(size(yy_iR,1)),' profiles)'],...
                        ['A-profiles (',num2str(size(yy_iA,1)),' profiles)'],...
                        ['D-profiles (',num2str(size(yy_iD,1)),' profiles)'],...
                      'Location','northeastoutside');
            set(lh,'FontSize',10);
            
%             xlimits=xlim();
%             xlim([xlimits(1) xlimits(2)*1.3]);

            % background color
            set(gcf,'color','w');
            % grid in y axis
            ax = gca;
            ax.YGrid = 'on';

             % save figure
            out_name = [output_plots_dir '/' sprintf('%02d',i_fig) '_' project_name '_' char(i_param_str) '_RAD_mode_per_wmo_per_cycle_' sprintf('%03d',i_graph) '_' working_date];
            disp(['saving ' out_name])
            %     export_fig([out_name '.png'])
            print('-dpng ', '-r100',[out_name '.png'])
            if print_svg == 1
                saveas(gcf,[out_name '.svg'])
            end
            
        end        

    end  
    
end
diary off
%% 
diary on
if output_graphs_per_float ==1
    
    disp('Profile QC status by float and by cycle number')
    list_QCs=["A";"B";"C";"D";"E";"F";"X"];
    i_fig=12;

    for i=1:n_param
        close all
        i_param=list_of_parameters_to_treat(i);
        
        if i_param == "PRES" 
            % PRES profie QC information is not yet available in the index
            % file
            continue
        end

        
        
        fprintf('Making Plots for %s \n',i_param)
                
        subset_prof_WMO=IndexData.profile_WMO(IndexData.param.(i_param).presence==1);
        subset_prof_cycle=IndexData.cycle(IndexData.param.(i_param).presence==1);
        subset_prof_qc=IndexData.param.(i_param).qc(IndexData.param.(i_param).presence==1);
        
        list_wmo_for_graph=unique(subset_prof_WMO);
        n_wmo_param=size(list_wmo_for_graph,1);
        n_graph=ceil(n_wmo_param/n_max_float_per_graph);
        
        
        for i_graph=1:n_graph
            
            i_wmo_min=1+(i_graph-1)*n_max_float_per_graph;
            i_wmo_max=min(n_wmo_param,i_graph*n_max_float_per_graph);
            list_wmo_for_i_graph=list_wmo_for_graph(i_wmo_min:i_wmo_max);
            n_wmo_for_i_graph=size(list_wmo_for_i_graph,1);
            
            prof_subset=ismember(subset_prof_WMO,list_wmo_for_i_graph);
            
            xx=double(subset_prof_cycle((prof_subset==1)));
            yy=subset_prof_WMO((prof_subset==1));
            zz=string(subset_prof_qc((prof_subset==1)));
            
            close all
            figure(i_fig)
            

            % bigger figure
            set(gcf, 'Position', [200, 200, 1200, 600])
            % figure name
            set(gcf,'Name','Profile QC status by float and by cycle number')
            hold on
            
            % title with update date
            title(['Profile QC status by float and by cycle for ' char(i_param) ' (updated ',update_date_str,')' ...
                  newline project_name ' - batch ' num2str(i_graph) '/' num2str(n_graph)], 'Interpreter', 'none')
            xlabel('Cycle Number')
            
            for iQC=1:size(list_QCs,1)
                
                i_QC=list_QCs(iQC);
                % find profiles indices where QC = i_QC
                ii_QC=(zz==i_QC);
                
                % retrieve corresponding wmo and cycle
                xx_iQC=xx(ii_QC);
                [~,yy_iQC]=ismember(yy(ii_QC),list_wmo_for_i_graph);
                
                % plot using scatter plot
                scatter(xx_iQC,yy_iQC,10,'o','filled','MarkerFaceColor',bars_colors(iQC+4,:))
                
                % prepare the legend
                if i_QC ~= "X"
                    lgd.(i_QC)=['profile QC ' char(i_QC) ' (' sprintf('%d',size(yy_iQC,1)) ' profiles)'];
                else
                    lgd.(i_QC)=['no profile QC (' sprintf('%d',size(yy_iQC,1)) ' profiles)'];
                end
                
            end
            
            set(gca,'ytick',1:n_wmo_for_i_graph,'yticklabel', list_wmo_for_i_graph')
            % legend
            lh = legend(lgd.A,lgd.B,lgd.C,lgd.D,lgd.E,lgd.F,lgd.X, ...
                      'Location','northeastoutside');
            set(lh,'FontSize',10);

%             xlimits=xlim();
%             xlim([xlimits(1) xlimits(2)*1.3]);
            
            % background color
            set(gcf,'color','w');
            % grid in y axis
            ax = gca;
            ax.YGrid = 'on';

             % save figure
            out_name = [output_plots_dir '/' sprintf('%02d',i_fig) '_' project_name '_' char(i_param) '_profile_QC_per_wmo_per_cycle_' sprintf('%03d',i_graph) '_' working_date];
            disp(['saving ' out_name])
            %     export_fig([out_name '.png'])
            print('-dpng ', '-r100',[out_name '.png'])
            if print_svg == 1
                saveas(gcf,[out_name '.svg'])
            end
            
        end        

    end 
end


diary off
%% 
diary on
if output_graphs_per_float ==1
    
    disp('PSAL adjustment by float and by cycle number')
    
    i_fig=13;
    
    i_param = "PSAL";            
    subset_prof_WMO=IndexData.profile_WMO(IndexData.param.(i_param).presence==1);
    subset_prof_cycle=IndexData.cycle(IndexData.param.(i_param).presence==1);
    subset_prof_adj=IndexData.param.(i_param).adj(IndexData.param.(i_param).presence==1);
    subset_prof_qc=IndexData.param.(i_param).qc(IndexData.param.(i_param).presence==1);
    subset_prof_mode=IndexData.param.(i_param).mode(IndexData.param.(i_param).presence==1);

    list_wmo_for_graph=unique(subset_prof_WMO);
    n_wmo_param=size(list_wmo_for_graph,1);
    n_graph=ceil(n_wmo_param/n_max_float_per_graph);


    for i_graph=1:n_graph

        i_wmo_min=1+(i_graph-1)*n_max_float_per_graph;
        i_wmo_max=min(n_wmo_param,i_graph*n_max_float_per_graph);
        list_wmo_for_i_graph=list_wmo_for_graph(i_wmo_min:i_wmo_max);
        n_wmo_for_i_graph=size(list_wmo_for_i_graph,1);

        prof_subset=ismember(subset_prof_WMO,list_wmo_for_i_graph);

        xx=double(subset_prof_cycle((prof_subset==1)));
        yy=subset_prof_WMO((prof_subset==1));
        zz=double(subset_prof_adj((prof_subset==1)));
        zz_qc=string(subset_prof_qc((prof_subset==1)));
        zz_mode=string(subset_prof_mode((prof_subset==1)));
        close all
        figure(i_fig)


        % bigger figure
        set(gcf, 'Position', [200, 200, 1000, 600])
        % figure name
        set(gcf,'Name','Profile QC status by float and by cycle number')
        hold on

        % title with update date
        title(['PSAL adjustment by float and by cycle for ' char(i_param) ' (updated ',update_date_str,')' ...
              newline project_name ' - batch ' num2str(i_graph) '/' num2str(n_graph)], 'Interpreter', 'none')
        xlabel('Cycle Number')
        
        
        % Plot R and A-profile in light grey
        % find profiles indices where QC = i_QC
        ii_RA=( ((zz_mode=="R") | (zz_mode=="A")) & (zz_qc~="F") );
        % retrieve corresponding wmo and cycle
        xx_iRA=xx(ii_RA);
        [~,yy_iRA]=ismember(yy(ii_RA),list_wmo_for_i_graph);
        % plot using scatter plot
        scatter(xx_iRA,yy_iRA,10,'o','filled','MarkerFaceColor',[0.78 0.78 0.78])
        % output the legend
        lgd1='R and A profiles (not QC F)';

        % Plot R and A-profile QC F in medium grey
        % find profiles indices where QC = i_QC
        ii_RAF=( ((zz_mode=="R") | (zz_mode=="A")) & (zz_qc=="F"));
        % retrieve corresponding wmo and cycle
        xx_iRAF=xx(ii_RAF);
        [~,yy_iRAF]=ismember(yy(ii_RAF),list_wmo_for_i_graph);
        % plot using scatter plot
        scatter(xx_iRAF,yy_iRAF,10,'o','filled','MarkerFaceColor',[0.4 0.4 0.4])
        % output the legend
        lgd2='R and A profiles with QC F';
  
        
        % Plot D-profile QC F in black
        % find profiles indices where QC = i_QC
        ii_DF=((zz_mode=="D") & (zz_qc=="F"));
        % retrieve corresponding wmo and cycle
        xx_iDF=xx(ii_DF);
        [~,yy_iDF]=ismember(yy(ii_DF),list_wmo_for_i_graph);
        % plot using scatter plot
        scatter(xx_iDF,yy_iDF,10,'o','filled','MarkerFaceColor','black')
        % output the legend
        lgd3='D profiles with QC F';
        [lh,~] = legend(lgd1,lgd2,lgd3,'Location','northeast');
        set(lh,'FontSize',10);
        set(lh,'Position',[0.6 0.85 0.23 0.09]);
        

        [~,yy_plot]=ismember(yy,list_wmo_for_i_graph);
        % plot using scatter plot
        scatter(xx,yy_plot,10,zz,'o','filled')
        colormap(gca,flipud(jet))
        colorbar()
        caxis([-0.07 0.07])
        
        % As information can be overlaid (both adj and R, both QC F and R):
        % rearrangement of the order of apparition
        scatter(xx_iRA,yy_iRA,10,'o','filled','MarkerFaceColor',[0.78 0.78 0.78])
        scatter(xx_iRAF,yy_iRAF,10,'o','filled','MarkerFaceColor',[0.4 0.4 0.4])
        scatter(xx_iDF,yy_iDF,10,'o','filled','MarkerFaceColor','black')
        

        set(gca,'ytick',1:n_wmo_for_i_graph,'yticklabel', list_wmo_for_i_graph')
        
        ylimits=ylim();
        ylim([ylimits(1) ylimits(2)*1.15]);
        
        % background color
        set(gcf,'color','w');
        % grid in y axis
        ax = gca;
        ax.YGrid = 'on';
        ax.XGrid = 'on';
        
        % Add minor ticks/grid along the x axis.
        xt=ax.XTick;                   % and the current tick values
        nT=length(xt);                 % how many ticks are there???
        nMinorT=5;                     % set how many minor tick divisions wanted
        ax.XAxis.MinorTickValues=linspace(xt(1),xt(end),(nT-1)*nMinorT+1); % set those values

        ax.XMinorTick = 'on';
        ax.XMinorGrid = 'on';
        

         % save figure
        out_name = [output_plots_dir '/' sprintf('%02d',i_fig) '_' project_name '_' char(i_param) '_PSAL_adj_per_wmo_per_cycle_' sprintf('%03d',i_graph) '_' working_date];
        disp(['saving ' out_name])
        %     export_fig([out_name '.png'])
        print('-dpng ', '-r100',[out_name '.png'])

        if print_svg == 1
            saveas(gcf,[out_name '.svg'])
        end
            
        

    end 
end

close all

diary off
%% 
diary on

% Create output directory for text analyses:
output_syntheses_dir = [output_dir '/Syntheses'];

if ~exist(output_syntheses_dir, 'dir')
    disp('creating the output_directory')
    mkdir(output_syntheses_dir)
end
disp(['  Synthetic analyses will be saved in '  output_syntheses_dir])

diary off
%% 
diary on
% Warnings
% 
for i=1:n_param
    i_param=list_of_parameters_to_treat(i);
    fprintf('\nOutputing warnings for %s \n',i_param)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % File 4: warnings by wmo (numbered 4 for legacy reasons)
    fprintf('\nRecording warnings by WMOs for %s \n',i_param)

    outfile = [output_syntheses_dir '/' 'DMQC_status_warnings_per_wmo_for_' char(i_param) '_' IndexData.index_update(1:8) '.txt'];
    fprintf('\nSaving results in %s ...\n',outfile)

    fid=fopen(outfile,'w');
    
    fprintf(fid, ['# Warnings for ' char(i_param) '\n']);
    fprintf(fid, '# Project : %s\n', project_name);
    fprintf(fid, '# Update date : %s\n', IndexData.index_update);    
        

    % floats with F cycles
    wmos_with_F_profiles.(i_param)=unique(IndexData.profile_WMO(IndexData.param.(i_param).presence==1 & ...
                                          IndexData.param.(i_param).qc=="F"));

    % floats with "X" cycles
    wmos_with_X_profiles.(i_param)=unique(IndexData.profile_WMO(IndexData.param.(i_param).presence==1 & ...
                                          IndexData.param.(i_param).qc=="X"));

    % WMO not found in index (reasons may be : 
    %                 - no ascending profile 
    %                 - no profile at all
    %                 - core profile but missing bgc profile in bgc case) 
    % => the list corresponds to the ones for which no profile are found
    % (see table 1 of the text output.

    % old profiles never processed in delayed mode

    wmos_with_D_profile.(i_param)=unique(IndexData.profile_WMO(IndexData.param.(i_param).mode == 'D'));
    wmos_older_than_sage.(i_param)=unique(IndexData.profile_WMO(i_old_profiles & IndexData.param.(i_param).presence==1));
    i_list=~ismember(wmos_older_than_sage.(i_param),wmos_with_D_profile.(i_param));
    tmp=wmos_older_than_sage.(i_param);
    wmos_older_than_sage_never_DMQCed.(i_param)=tmp(i_list);


    % Active floats in exclusion list
    wmos_with_RA_profiles_in_exclusion_list.(i_param) = string(exclusion_list.PLATFORMCODE(i_wmo_ope_in_exclusionlist_with_qc_3_or_4.(i_param) | ...
                                                          i_wmo_inactiveR_in_exclusionlist_with_qc_3_or_4.(i_param),:));
    
    % list all WMOs in warnings
    wmos_in_warning.(i_param) = unique([wmos_with_F_profiles.(i_param); ...
                                        wmos_with_X_profiles.(i_param); ...
                                        wmos_older_than_sage_never_DMQCed.(i_param); ...
                                        wmos_with_RA_profiles_in_exclusion_list.(i_param)]);
    
    n_wmos_in_warning = length(wmos_in_warning.(i_param));
    
    if n_wmos_in_warning > 0
    
        Param_name=repmat(i_param,n_wmos_in_warning,1);
        [~,loc]=ismember(wmos_in_warning.(i_param),string(Floats_list.WMO));
        warning_country = string(Floats_list.COUNTRYCODE(loc,:));
        warning_Fprof   = strings(n_wmos_in_warning,1);
        warning_Xprof   = strings(n_wmos_in_warning,1);
        warning_notDM   = strings(n_wmos_in_warning,1);
        warning_excll   = strings(n_wmos_in_warning,1);
        warning_Fprof(contains(wmos_in_warning.(i_param), wmos_with_F_profiles.(i_param))) = 'X';
        warning_Xprof(contains(wmos_in_warning.(i_param), wmos_with_X_profiles.(i_param))) = 'X';
        warning_notDM(contains(wmos_in_warning.(i_param), wmos_older_than_sage_never_DMQCed.(i_param))) = 'X';
        warning_excll(contains(wmos_in_warning.(i_param), wmos_with_RA_profiles_in_exclusion_list.(i_param))) = 'X';
        
        header4=['PARAM,','WMO,','Country,','prof_QC=F_once,',...
                'prof_QC_not_filled_once,','No_DMQC_',profile_old_c,',',...
                'exclusion_listed'];
        
        fprintf(fid, '%s \n', header4);
            
        %table4 = [Param_name wmos_in_warning.(i_param)  warning_country  warning_Fprof  warning_Xprof  warning_notDM  warning_excll]';
        table4 = [cellstr(Param_name), ...
                  cellstr(wmos_in_warning.(i_param)), ...
                  cellstr(warning_country), ...
                  cellstr(warning_Fprof),...
                  cellstr(warning_Xprof),...
                  cellstr(warning_notDM),...
                  cellstr(warning_excll)]';


        fprintf(fid,'%s,%s,%s,%s,%s,%s,%s\n',table4{:});
        
    end
end


diary off
%% 
diary on

for i=1:n_param


    i_param=list_of_parameters_to_treat(i);
    
    
    % File 1: Statistics per float
    fprintf('\nComputing and recording by float stats for %s \n',i_param)

    outfile = [output_syntheses_dir '/' 'DMQC_status_per_wmo_for_' char(i_param) '_' IndexData.index_update(1:8) '.txt'];
    fprintf('\nSaving results in %s ...\n',outfile)

    fid=fopen(outfile,'w');

    % File header
    fprintf(fid, '# DMQC status and Data quality control statistics\n');
    fprintf(fid, '# Project : %s\n', project_name);
    fprintf(fid, '# Update date : %s\n', IndexData.index_update);
    fprintf(fid, ['# Statistics per float for parameter ' char(i_param) '\n']);

    Floats_list.(i_param).WMO=unique(IndexData.profile_WMO(IndexData.param.(i_param).presence==1));
    n_floats_par=size(Floats_list.(i_param).WMO,1);
    
    fprintf('For %s there are %d floats\n',i_param,n_floats_par)
    
    [~,loc]=ismember(string(Floats_list.(i_param).WMO),Floats_list.WMO);
    Floats_list.(i_param).COUNTRYCODE = Floats_list.COUNTRYCODE(loc,:);
    Floats_list.(i_param).PROGRAM = Floats_list.PROGRAM(loc,:);
    Floats_list.(i_param).LAUNCHDATE = Floats_list.LAUNCHDATE(loc,:);
    
    [~,loc]=ismember(string(Floats_list.(i_param).WMO),IndexData.profile_WMO);
    Floats_list.(i_param).DAC = IndexData.profile_dac(loc,:);
    
    
    Param_name=repmat(i_param,n_floats_par,1);
    
    % Initialisations:
    Floats_list.(i_param).more_1_yr = zeros(n_floats_par,1);
    Floats_list.(i_param).wmo_dm_done = zeros(n_floats_par,1);
    Floats_list.(i_param).exclusion_list = zeros(n_floats_par,1);
    Floats_list.(i_param).prof_number = zeros(n_floats_par,1);
    Floats_list.(i_param).prof_older_than_sage_number = zeros(n_floats_par,1);
    Floats_list.(i_param).prof_last_date = strings(n_floats_par,1);
    Floats_list.(i_param).prof_DMQCed_number = zeros(n_floats_par,1);
    Floats_list.(i_param).prof_DMQCed_last_update_date = strings(n_floats_par,1);
    Floats_list.(i_param).prof_noDM_older_than_sage_number = zeros(n_floats_par,1);
    Floats_list.(i_param).prof_last_DM_date = strings(n_floats_par,1);
    Floats_list.(i_param).prof_DMQCed_percent = zeros(n_floats_par,1);
    Floats_list.(i_param).prof_QC_A = zeros(n_floats_par,1);
    Floats_list.(i_param).prof_DMQCed_QC_A = zeros(n_floats_par,1);
    Floats_list.(i_param).prof_QC_B = zeros(n_floats_par,1);
    Floats_list.(i_param).prof_DMQCed_QC_B = zeros(n_floats_par,1);
    Floats_list.(i_param).prof_QC_C = zeros(n_floats_par,1);
    Floats_list.(i_param).prof_DMQCed_QC_C = zeros(n_floats_par,1);
    Floats_list.(i_param).prof_QC_D = zeros(n_floats_par,1);
    Floats_list.(i_param).prof_DMQCed_QC_D = zeros(n_floats_par,1);
    Floats_list.(i_param).prof_QC_E = zeros(n_floats_par,1);
    Floats_list.(i_param).prof_DMQCed_QC_E = zeros(n_floats_par,1);
    Floats_list.(i_param).prof_QC_F = zeros(n_floats_par,1);
    Floats_list.(i_param).prof_DMQCed_QC_F = zeros(n_floats_par,1);
    Floats_list.(i_param).prof_QC_X = zeros(n_floats_par,1);
    Floats_list.(i_param).prof_DMQCed_QC_X = zeros(n_floats_par,1);
    
    %'float_age,'
    Floats_list.(i_param).wmo_age = datenum(IndexData.index_update,'yyyymmddHHMMSS') - datenum(Floats_list.(i_param).LAUNCHDATE(:,1:10),'yyyy/mm/dd');

    %'float_more1year,'
    [~,loc]=ismember(wmos_older_than_sage.(i_param),string(Floats_list.(i_param).WMO));
    Floats_list.(i_param).more_1_yr(loc) = 1; 

    %'DM_done,'
    [~,loc]=ismember(wmos_with_D_profile.(i_param),string(Floats_list.(i_param).WMO));
    Floats_list.(i_param).wmo_dm_done(loc) = 1;

    %'exclusionlist,'
    [~,loc]=ismember(wmos_with_RA_profiles_in_exclusion_list.(i_param),string(Floats_list.(i_param).WMO));
    Floats_list.(i_param).exclusion_list(loc) = 1;

    %'profile_number,'
    [xx,~,ic]=unique(IndexData.profile_WMO(IndexData.param.(i_param).presence == 1));
    val=accumarray(ic,1);
    [~,loc]=ismember(xx,string(Floats_list.(i_param).WMO));
    Floats_list.(i_param).prof_number(loc)=val;

    %'prof_older_than_sage_number,'
    [xx,~,ic]=unique(IndexData.profile_WMO(i_old_profiles & IndexData.param.(i_param).presence == 1));
    val=accumarray(ic,1);
    [~,loc]=ismember(xx,string(Floats_list.(i_param).WMO));
    Floats_list.(i_param).prof_older_than_sage_number(loc)=val;
    
    %'date_last_prof,'
    loc_wmos  = IndexData.profile_WMO(IndexData.param.(i_param).presence == 1);
    loc_date = datenum(IndexData.date(IndexData.param.(i_param).presence == 1),'yyyymmddHHMMSS');
    % Index the positions of each wmo (here 'stable' is very important for
    % the accumarray to work as expected)
    [xx,~,ic] = unique(loc_wmos,'stable');
    % compute the max for each wmo.
    val=datestr(accumarray(ic,loc_date,[],@max),'yyyy/mm/dd HH:MM:SS'); 
    [~,loc]=ismember(xx,string(Floats_list.(i_param).WMO));
    Floats_list.(i_param).prof_last_date(loc)=string(val);
    
    

    %'number_DMprof,'
    [xx,~,ic]=unique(IndexData.profile_WMO(IndexData.param.(i_param).mode == 'D'));
    val=accumarray(ic,1);
    [~,loc]=ismember(xx,string(Floats_list.(i_param).WMO));
    Floats_list.(i_param).prof_DMQCed_number(loc)=val;

    %'obs_more1year_noDM,'
    [xx,~,ic]=unique(IndexData.profile_WMO(IndexData.param.(i_param).mode ~= 'D' & ...
                                           IndexData.param.(i_param).presence == 1 & ...
                                           i_old_profiles));
    val=accumarray(ic,1);
    [~,loc]=ismember(xx,string(Floats_list.(i_param).WMO));
    Floats_list.(i_param).prof_noDM_older_than_sage_number(loc)=val;

    %'date_last_DMprof,'
    loc_wmos  = IndexData.profile_WMO(IndexData.param.(i_param).mode == 'D');
    loc_date = datenum(IndexData.date(IndexData.param.(i_param).mode == 'D'),'yyyymmddHHMMSS');
    % Index the positions of each wmo (here 'stable' is very important for
    % the accumarray to work as expected)
    [xx,~,ic] = unique(loc_wmos,'stable');
    % compute the max for each wmo.
    val=datestr(accumarray(ic,loc_date,[],@max),'yyyy/mm/dd HH:MM:SS'); 
    [~,loc]=ismember(xx,string(Floats_list.(i_param).WMO));
    Floats_list.(i_param).prof_last_DM_date(loc)=string(val);
    
    % update_date_last_DMprof
    i_index_loc= IndexData.param.(i_param).presence == 1 & ...
                 IndexData.param.(i_param).mode == 'D';
    if (sum(i_index_loc) > 0) 
        loc_wmos  = IndexData.profile_WMO(i_index_loc);
        loc_date = datenum(IndexData.prof_update_date(i_index_loc),'yyyymmddHHMMSS');
        % Index the positions of each wmo (here 'stable' is very important for
        % the accumarray to work as expected)
        [xx,~,ic] = unique(loc_wmos,'stable');
        % compute the max for each wmo.
        val=datestr(accumarray(ic,loc_date,[],@max),'yyyy/mm/dd'); 
        [~,loc]=ismember(xx,string(Floats_list.(i_param).WMO));
        Floats_list.(i_param).prof_DMQCed_last_update_date(loc)=string(val);
    end
    
    %'percentage_DMobs,'
    i_ok=(Floats_list.(i_param).prof_number ~= 0);
    Floats_list.(i_param).prof_DMQCed_percent(i_ok) = 100 * Floats_list.(i_param).prof_DMQCed_number(i_ok) ./ Floats_list.(i_param).prof_number(i_ok);
    
    %'profil_QC_A,'
    [xx,~,ic]=unique(IndexData.profile_WMO(IndexData.param.(i_param).qc == "A"));
    val=accumarray(ic,1);
    [~,loc]=ismember(xx,string(Floats_list.(i_param).WMO));
    Floats_list.(i_param).prof_QC_A(loc)=val;
    
    [xx,~,ic]=unique(IndexData.profile_WMO(IndexData.param.(i_param).qc == "A" & ...
                                           IndexData.param.(i_param).mode == 'D')) ;
    val=accumarray(ic,1);
    [~,loc]=ismember(xx,string(Floats_list.(i_param).WMO));
    Floats_list.(i_param).prof_DMQCed_QC_A(loc)=val;
    
    %'profil_QC_B,'
    [xx,~,ic]=unique(IndexData.profile_WMO(IndexData.param.(i_param).qc == "B"));
    val=accumarray(ic,1);
    [~,loc]=ismember(xx,string(Floats_list.(i_param).WMO));
    Floats_list.(i_param).prof_QC_B(loc)=val;
    
    [xx,~,ic]=unique(IndexData.profile_WMO(IndexData.param.(i_param).qc == "B" & ...
                                           IndexData.param.(i_param).mode == 'D')) ;
    val=accumarray(ic,1);
    [~,loc]=ismember(xx,string(Floats_list.(i_param).WMO));
    Floats_list.(i_param).prof_DMQCed_QC_B(loc)=val;
    
    %'profil_QC_C,'
    [xx,~,ic]=unique(IndexData.profile_WMO(IndexData.param.(i_param).qc == "C"));
    val=accumarray(ic,1);
    [~,loc]=ismember(xx,string(Floats_list.(i_param).WMO));
    Floats_list.(i_param).prof_QC_C(loc)=val;
    
    [xx,~,ic]=unique(IndexData.profile_WMO(IndexData.param.(i_param).qc == "C" & ...
                                           IndexData.param.(i_param).mode == 'D')) ;
    val=accumarray(ic,1);
    [~,loc]=ismember(xx,string(Floats_list.(i_param).WMO));
    Floats_list.(i_param).prof_DMQCed_QC_C(loc)=val;
    
    
     %'profil_QC_D,'
    [xx,~,ic]=unique(IndexData.profile_WMO(IndexData.param.(i_param).qc == "D"));
    val=accumarray(ic,1);
    [~,loc]=ismember(xx,string(Floats_list.(i_param).WMO));
    Floats_list.(i_param).prof_QC_D(loc)=val;
    
    [xx,~,ic]=unique(IndexData.profile_WMO(IndexData.param.(i_param).qc == "D" & ...
                                           IndexData.param.(i_param).mode == 'D')) ;
    val=accumarray(ic,1);
    [~,loc]=ismember(xx,string(Floats_list.(i_param).WMO));
    Floats_list.(i_param).prof_DMQCed_QC_D(loc)=val;

    %'profil_QC_E,'
    [xx,~,ic]=unique(IndexData.profile_WMO(IndexData.param.(i_param).qc == "E"));
    val=accumarray(ic,1);
    [~,loc]=ismember(xx,string(Floats_list.(i_param).WMO));
    Floats_list.(i_param).prof_QC_E(loc)=val;
    
    [xx,~,ic]=unique(IndexData.profile_WMO(IndexData.param.(i_param).qc == "E" & ...
                                           IndexData.param.(i_param).mode == 'D')) ;
    val=accumarray(ic,1);
    [~,loc]=ismember(xx,string(Floats_list.(i_param).WMO));
    Floats_list.(i_param).prof_DMQCed_QC_E(loc)=val;
    
    %'profil_QC_F,'
    [xx,~,ic]=unique(IndexData.profile_WMO(IndexData.param.(i_param).qc == "F"));
    val=accumarray(ic,1);
    [~,loc]=ismember(xx,string(Floats_list.(i_param).WMO));
    Floats_list.(i_param).prof_QC_F(loc)=val;
    
    [xx,~,ic]=unique(IndexData.profile_WMO(IndexData.param.(i_param).qc == "F" & ...
                                           IndexData.param.(i_param).mode == 'D')) ;
    val=accumarray(ic,1);
    [~,loc]=ismember(xx,string(Floats_list.(i_param).WMO));
    Floats_list.(i_param).prof_DMQCed_QC_F(loc)=val;
    
    %'profil_QC_X,'
    [xx,~,ic]=unique(IndexData.profile_WMO(IndexData.param.(i_param).qc == "X" & ...
                                           IndexData.param.(i_param).presence ==1)); % maybe review the replace "" by "X" ...
    val=accumarray(ic,1);
    [~,loc]=ismember(xx,string(Floats_list.(i_param).WMO));
    Floats_list.(i_param).prof_QC_X(loc)=val;
    
    [xx,~,ic]=unique(IndexData.profile_WMO(IndexData.param.(i_param).qc == "X" & ...
                                           IndexData.param.(i_param).presence ==1 & ...
                                           IndexData.param.(i_param).mode == 'D')) ;
    val=accumarray(ic,1);
    [~,loc]=ismember(xx,string(Floats_list.(i_param).WMO));
    Floats_list.(i_param).prof_DMQCed_QC_X(loc)=val;


    % header
    header1 = ['PARAM,' 'WMO,' 'DAC,' 'COUNTRY,' 'PROGRAM,' 'launch_date,' ...
               'prof_last_date,' 'DM_done,'  ...
               'float_age,' 'float_' profile_old_c ',' 'exclusionlist,' 'nb_prof,' ...
               'nb_prof_' profile_old_c ',' 'nb_prof_DM,' 'nb_prof_noDM_' profile_old_c ',' ...
               'prof_last_DM_date,' 'prof_last_DM_update_date,' 'percentage_DM_prof,'...
               'nb_prof_QC_A,' 'nb_prof_DM_QC_A,' 'nb_prof_QC_B,' 'nb_prof_DM_QC_B,' ...
               'nb_prof_QC_C,' 'nb_prof_DM_QC_C,' 'nb_prof_QC_D,' 'nb_prof_DM_QC_D,' ...
               'nb_prof_QC_E,' 'nb_prof_DM_QC_E,' 'nb_prof_QC_F,' 'nb_prof_DM_QC_F,' ...
               'nb_prof_QC_X,' 'nb_prof_DM_QC_X'];
    fprintf(fid, '%s \n', header1);
    % data
    table1 = [cellstr(Param_name), ...
              cellstr(Floats_list.(i_param).WMO), ...
              cellstr(Floats_list.(i_param).DAC), ...
              cellstr(Floats_list.(i_param).COUNTRYCODE),...
              cellstr(Floats_list.(i_param).PROGRAM),...
              cellstr(Floats_list.(i_param).LAUNCHDATE),...
              cellstr(Floats_list.(i_param).prof_last_date),...
              num2cell(Floats_list.(i_param).wmo_dm_done),...
              num2cell(Floats_list.(i_param).wmo_age),...
              num2cell(Floats_list.(i_param).more_1_yr),...
              num2cell(Floats_list.(i_param).exclusion_list),...
              num2cell(Floats_list.(i_param).prof_number),...
              num2cell(Floats_list.(i_param).prof_older_than_sage_number),...
              num2cell(Floats_list.(i_param).prof_DMQCed_number),...
              num2cell(Floats_list.(i_param).prof_noDM_older_than_sage_number),...
              cellstr(Floats_list.(i_param).prof_last_DM_date),...
              cellstr(Floats_list.(i_param).prof_DMQCed_last_update_date),...
              num2cell(Floats_list.(i_param).prof_DMQCed_percent),...
              num2cell(Floats_list.(i_param).prof_QC_A),...
              num2cell(Floats_list.(i_param).prof_DMQCed_QC_A),...
              num2cell(Floats_list.(i_param).prof_QC_B),...
              num2cell(Floats_list.(i_param).prof_DMQCed_QC_B),...
              num2cell(Floats_list.(i_param).prof_QC_C),...
              num2cell(Floats_list.(i_param).prof_DMQCed_QC_C),...
              num2cell(Floats_list.(i_param).prof_QC_D),...
              num2cell(Floats_list.(i_param).prof_DMQCed_QC_D),...
              num2cell(Floats_list.(i_param).prof_QC_E),...
              num2cell(Floats_list.(i_param).prof_DMQCed_QC_E),...
              num2cell(Floats_list.(i_param).prof_QC_F),...
              num2cell(Floats_list.(i_param).prof_DMQCed_QC_F),...
              num2cell(Floats_list.(i_param).prof_QC_X),...
              num2cell(Floats_list.(i_param).prof_DMQCed_QC_X)]';
              

    fprintf(fid, '%s,%s,%s,%s,"%s",%s,%s,%d,%.0f,%d,%d,%d,%d,%d,%d,%s,%s,%.2f,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n',table1{:});
    
    fclose(fid);
 
end

diary off
%% 
diary on
for i=1:n_param
    i_param=list_of_parameters_to_treat(i);
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % File 2 (and 3_total): Statistics per country (data from figures 1 and 2)
    fprintf('\nComputing and recording by country stats for %s \n',i_param)

    outfile = [output_syntheses_dir '/' 'DMQC_status_per_country_for_' char(i_param) '_' IndexData.index_update(1:8) '.txt'];
    fprintf('\nSaving results in %s ...\n',outfile)

    fid=fopen(outfile,'w');

    % File header
    fprintf(fid, '# DMQC status and Data quality control statistics\n');
    fprintf(fid, '# Project : %s\n', project_name);
    fprintf(fid, '# Update date : %s\n', IndexData.index_update);
    fprintf(fid, ['# Statistics per country for parameter ' char(i_param) '\n']);
    % 
    
    Param_name=repmat(i_param,n_countries,1);

    % header
    header2 = ['Param_name,' 'Country,' ...
                'nb_floats,' 'nb_floats_' profile_old_c ',' ...
                'Nb_floats_DM_done,' 'nb_floats_' profile_old_c '_DM_done,' ...
                'nb_profiles,' 'nb_profiles_' profile_old_c ',' ...
                'nb_profiles_DM_done,' 'nb_profiles_' profile_old_c '_DM_done,'];
    fprintf(fid, '%s \n', header2);
    % data
    table2 = [cellstr(Param_name),...
              cellstr(nb_per_country_x), ...
              num2cell(nb_floats_per_country.(i_param)), ...
              num2cell(nb_floats_xage_per_country.(i_param)),...
              num2cell(nb_floats_DMQCed_per_country.(i_param)),...
              num2cell(nb_floats_xage_DMQCed_per_country.(i_param)),...
              num2cell(nb_profiles_per_country.(i_param)), ...
              num2cell(nb_profiles_xage_per_country.(i_param)),...
              num2cell(nb_profiles_DMQCed_per_country.(i_param)),...
              num2cell(nb_profiles_xage_DMQCed_per_country.(i_param))]';
    
    fprintf(fid, '%s,%s,%d,%d,%d,%d,%d,%d,%d,%d\n',table2{:});
 
    
    table3 = [cellstr(i_param),...
              cellstr('All countries'), ...
              num2cell(sum(nb_floats_per_country.(i_param))), ...
              num2cell(sum(nb_floats_xage_per_country.(i_param))),...
              num2cell(sum(nb_floats_DMQCed_per_country.(i_param))),...
              num2cell(sum(nb_floats_xage_DMQCed_per_country.(i_param))),...
              num2cell(sum(nb_profiles_per_country.(i_param))), ...
              num2cell(sum(nb_profiles_xage_per_country.(i_param))),...
              num2cell(sum(nb_profiles_DMQCed_per_country.(i_param))),...
              num2cell(sum(nb_profiles_xage_DMQCed_per_country.(i_param)))]';
          
    fprintf(fid, '%s,%s,%d,%d,%d,%d,%d,%d,%d,%d\n',table3{:});   
          
    fclose(fid);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

diary off
%% 
diary on
for i=1:n_param
    i_param=list_of_parameters_to_treat(i);
    
%     if (nb_per_year_x(1) == string(1980))
%         nb_per_year_x_cr=nb_per_year_x(2:end);
%         nb_profiles_DMQCed_per_year_cr=nb_profiles_DMQCed_per_year.(i_param)(2:end);
%         nb_profiles_per_year_cr=nb_profiles_per_year.(i_param)(2:end);
%         nb_profiles_QCF_per_year_cr=nb_profiles_QCF_per_year.(i_param)(2:end);
%         nb_profiles_DMQCed_QCF_per_year_cr=nb_profiles_DMQCed_QCF_per_year.(i_param)(2:end);
%     else
        nb_per_year_x_cr=nb_per_year_x;
        nb_profiles_DMQCed_per_year_cr=nb_profiles_DMQCed_per_year.(i_param);
        nb_profiles_per_year_cr=nb_profiles_per_year.(i_param);
        nb_profiles_QCF_per_year_cr=nb_profiles_QCF_per_year.(i_param);
        nb_profiles_DMQCed_QCF_per_year_cr=nb_profiles_DMQCed_QCF_per_year.(i_param);
%     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Statistics per profile year (data for figures 7 and 8)
    fprintf('\nComputing and recording by profile year stats for %s \n',i_param)

    outfile = [output_syntheses_dir '/' 'DMQC_status_per_profile_year_for_' char(i_param) '_' IndexData.index_update(1:8) '.txt'];
    fprintf('\nSaving results in %s ...\n',outfile)

    fid=fopen(outfile,'w');

    % File header
    fprintf(fid, '# DMQC status and Data quality control statistics\n');
    fprintf(fid, '# Project : %s\n', project_name);
    fprintf(fid, '# Update date : %s\n', IndexData.index_update);
    fprintf(fid, ['# Statistics per profile year for parameter ' char(i_param) '\n']);
    % 
    
    Param_name=repmat(i_param,n_years,1);

    % header
    header2 = ['Param_name,' 'Year,' ...
                'nb_profiles,' 'nb_profiles_DMQCed,' ...
                'nb_profiles_QCF,' 'nb_profiles_DMQCed_QCF'];
    fprintf(fid, '%s \n', header2);
    % data
    table = [cellstr(Param_name),...
              cellstr(nb_per_year_x_cr), ...
              num2cell(nb_profiles_per_year_cr), ...
              num2cell(nb_profiles_DMQCed_per_year_cr),...
              num2cell(nb_profiles_QCF_per_year_cr),...
              num2cell(nb_profiles_DMQCed_QCF_per_year_cr)]';
    
    fprintf(fid, '%s,%s,%d,%d,%d,%d\n',table{:});

          
    fclose(fid);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

diary off
%% 
diary on

fprintf('\nAdditional information\n')

outfile = [output_syntheses_dir '/' 'Additional_Info_' IndexData.index_update(1:8) '.txt'];
fprintf('\nSaving results in %s ...\n',outfile)

fid=fopen(outfile,'w');
fprintf(fid, '\n');

if i_bgc == 1
    if size(IndexData.ParamNotFound_In_BGC_Index_Requested,1) - 1 == 0
        fprintf(fid, '# All requested BGC parameters were found in the synthetic index file\n');
        fprintf(fid, '\n');
    else
        fprintf(fid, '# List of requested BGC parameters that were not found in the synthetic index file\n');
        table5 = cellstr(IndexData.ParamNotFound_In_BGC_Index_Requested)';
        fprintf(fid, '%s\n',table5{:});   
        fprintf(fid, '\n');   
    end

    if size(IndexData.ParamFound_In_BGC_Index_notRequested,1) - 1 == 0
        fprintf(fid, '# All parameters found in the synthetic index file for the wmo list were requested\n');
        fprintf(fid, '\n');
    else
        fprintf(fid, '# List of parameters found in the synthetic index file (and associated to the float list) that were not requested\n');
        table6 = cellstr(IndexData.ParamFound_In_BGC_Index_notRequested)';
        fprintf(fid, '%s\n',table6{:});   
        fprintf(fid, '\n'); 
    end
end


if i_bgc ==1
    unfound_CTD_WMOs=IndexData.WMOs_unfound_in_IndexFile.CTD;
else
    unfound_CTD_WMOs=IndexData.WMOs_unfound_in_IndexFile;
end

if size(unfound_CTD_WMOs,1)==0
    fprintf(fid, '# All requested WMOs were found in the index file (CTD)\n');
    fprintf(fid, '\n');
else
    fprintf(fid, '# List of WMOs not found in the index file (CTD) \n');
    table7 = cellstr(unfound_CTD_WMOs)';
    fprintf(fid, '%s\n',table7{:});   
    fprintf(fid, '\n'); 
end

if i_bgc == 1
    if size(IndexData.WMOs_unfound_in_IndexFile.BGC,1)==0
        fprintf(fid, '# All requested WMOs were found in the synthetic index file (BGC)\n');
        fprintf(fid, '\n');
    else
        fprintf(fid, '# List of WMOs not found in the synthetic index file (BGC) \n');
        table7 = cellstr(IndexData.WMOs_unfound_in_IndexFile.BGC)';
        fprintf(fid, '%s\n',table7{:});   
        fprintf(fid, '\n'); 
    end
end

fclose(fid);  
    
diary off
%% 
diary on
% Zipping the local index file:
disp('- Zipping the local copy of the index file ...')

if exist(local_index_file,'file')
    gzip(local_index_file)
    delete(local_index_file)
end

try
    if exist(local_index_synthetic_file,'file') 
        gzip(local_index_synthetic_file)
        delete(local_index_synthetic_file)
    end
catch
    % the local_index_synthetic_file variable does not exist
end

disp('- Zipping ended ...')
diary off
