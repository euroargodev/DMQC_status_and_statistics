function [IndexData] = get_data_from_index(index_file, nb_header_lines,float_ref, i_descending_profile,i_bgc,list_of_parameters_to_treat)
% get_data_from_index gets relevant variables from index file for later
% processing
%
% Inputs
% - index_file           : index file path
% - float_ref            : cell array with float WMO (exemple: {'3901875','3901872'})
% - i_descending_profile : (optional) if 1 descent profiles are consider, if 0 descent profiles
%                          are not consider. Default value: 1.
% - i_bgc                : (optional) if the synthetic index is to be read (1) or not (0)
% - list_of_parameters_to_treat: (optional) ["TEMP"; "PRES"; "PSAL"; "DOXY];
%
% Outputs
% - IndexData: struct with all variables in file as fields.
%
% WARNING : Profile_QC for PRES information is not yet available. IndexData
% is filled with qc="X" for the moment, and plots related to pres profile
% qc are skipped.
%
% NOTE
% the archistructure has been reviewed in release 2.0: the
% information is now recorded by profile/parameter. When needed, the 
% reconstruction by wmo is done in get_DMQC_stats.m through "sql-like" 
% Matlab functions.
%
% Author: Euro-Argo ERIC (contact@euro-argo.eu)
%
% Version: 3.1 (2023/07/13)
%
%
% Historic:
% V1.0  (2022): 
%        - This script originally created by Andrea Garcia Juan and Romain Cancouët.
% v2.0  (2023/06/19):
%        - The script architecture was reviewed on 2023/06/19 by Delphine Dobler 
%          to include the processing of BGC floats and to enhance performances.
% v2.01 (2023/07/07): 
%        - correct typo in list_of_parameters corrected in if nargin ==2 loop.
%        - extract the list of parameter to treat by subsetting the list_of_parameters_to_treat
%           with the parameters effectively found in the index for float_ref list
%           (and output information notice)
% v2.1  (2023/07/07): 
%        - add update_date information
%        - add DAC fetching from the index instead of requiring an input.
%        - record unfound wmos in the output for file recording.
%        - initiate list_of_parameters output in the core index case
% v3.1  (2023/07/13):
%        - add ad_psal_adjustement_mean to the function output (i_bgc=0
%          case)
%        - fill prof_pres_qc with X instead of prof_temp_qc values, not to
%        mislead the graph reader (prof_pres_qc not available in index
%        yet).
%        - add number of header lines in the function argument (cases when 
%        index file optimization is not possible).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% no descending profiles defaut value
if nargin == 3
    i_descending_profile = 0;
    i_bgc=0;
    list_of_parameters_to_treat=["TEMP"; "PRES"; "PSAL"];
end

%% extracting the number of floats from the list of wmos provided in the input file.

float_ref_char = char(float_ref);

% Lets try a different approach, maybe more efficient way of reading a structured
% file :). importdata is interesting however the separation into text and
% double type data is not customizable which is annoying. readtable seems
% to be more what is wanted.
% there seems to be a size limit. 2,8m of lines seems too many for
% importdata, the same for 606160 (European subset)
% the readtable routine seems more robust, and was able to tackle the
% 606160 lines of the European subset in 1,7 s.


disp('---- Loading the index file in a Matlab table and convert to string format')
RT_OUT = readtable(index_file,'Headerlines',nb_header_lines);
RT_OUT_column_name= RT_OUT.Properties.VariableNames;
RT_OUT = string(table2cell(RT_OUT));


disp('---- Let find the column index for relevant piece of information')
ic_file=(string(RT_OUT_column_name)=='file');
ic_profile_date=(string(RT_OUT_column_name)=='date');
ic_longitude=(string(RT_OUT_column_name)=='longitude');
ic_latitude=(string(RT_OUT_column_name)=='latitude');
ic_prof_update_date=(string(RT_OUT_column_name)=='date_update');
if i_bgc==0 
    %ic_pres_profqc=(string(index_column_name)=='profile_pres_qc');% not
    %present in index file ... yet: record temp qc in the meanwhile
    ic_pres_profqc=(string(RT_OUT_column_name)=='profile_temp_qc');
    ic_temp_profqc=(string(RT_OUT_column_name)=='profile_temp_qc');
    ic_psal_profqc=(string(RT_OUT_column_name)=='profile_psal_qc');
    ic_psal_adj=(string(RT_OUT_column_name)=='ad_psal_adjustment_mean');
else
    ic_param_names=(string(RT_OUT_column_name)=='parameters');
    ic_param_modes=(string(RT_OUT_column_name)=='parameter_data_mode');
    ic_param_profqc=(string(RT_OUT_column_name)=='parameter_quality');
end


disp('---- Extract index_wmos, index_dac and index_cycles:')

% index_wmos
RT_OUT_profile_files_fullpath = RT_OUT(:,ic_file);
RT_OUT_profile_path=split(RT_OUT_profile_files_fullpath,"/");
RT_OUT_wmos=RT_OUT_profile_path(:,2);

% index_dac
RT_OUT_dacs=RT_OUT_profile_path(:,1);

% index_cycles
tmp=split(RT_OUT_profile_path(:,4),"_");
tmp=tmp(:,2);
tmp=split(tmp,".");
RT_OUT_cycles=tmp(:,1);
clearvars tmp

disp('---- Finding ascending and descending profiles and indexing relevant lines from the input file')
% categorize ascending/desceding cycles
% i_descending_cycles=contains(index_cycles,"D");
i_ascending_cycles=~contains(RT_OUT_cycles,"D");

% subsetting the index information for lines concerning wmos provided in input 
% and accounting for ascending or descending profiles option.
i_index_wmos_in_wmo_list=ismember(RT_OUT_wmos,float_ref_char);
if i_descending_profile
    i_relevant_index_lines=i_index_wmos_in_wmo_list;
else
    i_relevant_index_lines=i_ascending_cycles & i_index_wmos_in_wmo_list;
end

disp('---- Recording wmo, cycle number, date, lat, lon of relevant profiles')
IndexData.profile_WMO = RT_OUT_wmos(i_relevant_index_lines);
IndexData.profile_dac = RT_OUT_dacs(i_relevant_index_lines);
IndexData.cycle       = RT_OUT_cycles(i_relevant_index_lines);
IndexData.date        = RT_OUT(i_relevant_index_lines,ic_profile_date);
IndexData.latitude    = RT_OUT(i_relevant_index_lines,ic_latitude);
IndexData.longitude   = RT_OUT(i_relevant_index_lines,ic_longitude);
IndexData.prof_update_date   = RT_OUT(i_relevant_index_lines,ic_prof_update_date);
if i_descending_profile == 0
    IndexData.comment='only ascending profiles';
else
    IndexData.comment='including ascending and descending profiles';
end


disp('---- finding and recording mode and QC depending on index file source')


RT_OUT_profile_files_char=char(RT_OUT_profile_path(:,4));

if i_bgc==0 
   
    %even if the mode is the same, in analogy with other BGC parameters, it
    %is recorded. Maybe changed later on.
    IndexData.param.PSAL.mode = RT_OUT_profile_files_char(i_relevant_index_lines,1);
    IndexData.param.TEMP.mode = RT_OUT_profile_files_char(i_relevant_index_lines,1);
    IndexData.param.PRES.mode = RT_OUT_profile_files_char(i_relevant_index_lines,1);
    IndexData.param.PSAL.qc   = RT_OUT(i_relevant_index_lines,ic_psal_profqc);
    IndexData.param.TEMP.qc   = RT_OUT(i_relevant_index_lines,ic_temp_profqc);
    IndexData.param.PRES.qc   = RT_OUT(i_relevant_index_lines,ic_pres_profqc);
    % PRES_profile_qc not in the index yet. (To remove when index is
    % updated).
    IndexData.param.PRES.qc(:)   = "X";
    IndexData.param.PSAL.adj   = RT_OUT(i_relevant_index_lines,ic_psal_adj);
    IndexData.param.PSAL.nb_profile = sum(i_relevant_index_lines);
    IndexData.param.TEMP.nb_profile = sum(i_relevant_index_lines);
    IndexData.param.PRES.nb_profile = sum(i_relevant_index_lines);
    IndexData.param.PSAL.presence = ones(sum(i_relevant_index_lines),1);
    IndexData.param.TEMP.presence = ones(sum(i_relevant_index_lines),1);
    IndexData.param.PRES.presence = ones(sum(i_relevant_index_lines),1);
    
    
    
    % Those 3 are not necessary for core case but need to be initiated for 
    % later common processing
    IndexData.ParamNotFoundInIndexList="";
    IndexData.ParamFoundInIndex_notRequested="";
    IndexData.ParamList=list_of_parameters_to_treat;
    
else
    % in BGC index file case, parameters are grouped in a same field
    
    RT_OUT_param=RT_OUT(i_relevant_index_lines,ic_param_names);
    RT_OUT_param_ws=strings(size(RT_OUT_param,1),1);
    
    RT_OUT_param_mode=char(RT_OUT(i_relevant_index_lines,ic_param_modes));
    RT_OUT_param_profqc=char(RT_OUT(i_relevant_index_lines,ic_param_profqc));
    
    %-------------------------------------------------------------------------------
    % Here is a painful task to treat in an interpreted manner without looping:
    % the split function will only work with the
    % same number of delimiters + the arrangements in variables varies
    % within the index file ...
    % first count number of whitespaces for each index entry:
    RT_OUT_param_nb_whitespaces=count(RT_OUT_param," ");
    uniq_nb_whitespaces=unique(RT_OUT_param_nb_whitespaces);
    max_nb_of_ws=max(uniq_nb_whitespaces);
    %then for each case of number of parameters: fill in with whitespaces to complete.
    for i = 1:size(uniq_nb_whitespaces)
        l_nbws=uniq_nb_whitespaces(i);
        % compute the number of whitespaces to add
        nb_ws_to_add=max_nb_of_ws-l_nbws;
        
        
        %find entries with such number of ws:
        i_lines_to_complete=(RT_OUT_param_nb_whitespaces==l_nbws);
        nb_lines_to_complete=sum(i_lines_to_complete);
        
        % if it is already the max: only save in output matrix and go to
        % next
        if nb_ws_to_add==0
            RT_OUT_param_ws(i_lines_to_complete)=RT_OUT_param(i_lines_to_complete);
            continue
        end
        % else create a string array with the missing ws:
        if nb_ws_to_add >1
            ws_str_array=join(repmat(" ",nb_lines_to_complete,nb_ws_to_add),"");
        else
            ws_str_array=repmat(" ",nb_lines_to_complete,nb_ws_to_add);
        end
        
        % and add whitespaces at the end:
        RT_OUT_param_ws(i_lines_to_complete)=join(cat(2,RT_OUT_param(i_lines_to_complete),ws_str_array),"");
        
        %RT_OUT_param_nb_whitespaces_b=count(RT_OUT_param_ws(i_lines_to_complete)," ");
        
    end
    % now split using whitespace delimiter
    RT_OUT_PARAM_split=split(RT_OUT_param_ws," ");
    
    %replace DOXY2 by DOX2 and DOXY3 by DOX3, otherise contains method will
    %find them. Matlab is not an easy tool for whole-word search.
    RT_OUT_PARAM_split=replace(RT_OUT_PARAM_split,"DOXY2","DOX2");
    RT_OUT_PARAM_split=replace(RT_OUT_PARAM_split,"DOXY3","DOX3");
    %---------------------------------------------------------------------------------
    
    % now, let's find the list of parameters in the file:
    list_of_parameters=unique(reshape(RT_OUT_PARAM_split,size(RT_OUT_PARAM_split,1)*size(RT_OUT_PARAM_split,2),1));
    
    IndexData.ParamList=intersect(list_of_parameters,list_of_parameters_to_treat);
    IndexData.ParamNotFoundInIndexList=list_of_parameters_to_treat(~ismember(list_of_parameters_to_treat,list_of_parameters));
    IndexData.ParamFoundInIndex_notRequested=list_of_parameters(~ismember(list_of_parameters,list_of_parameters_to_treat));

    fprintf('%d requested parameters were not found in the index file:\n',size(IndexData.ParamNotFoundInIndexList,1))
    disp(IndexData.ParamNotFoundInIndexList)
    
    fprintf('%d parameters found in the index file (and associated to the float list) were not requested:\n',size(IndexData.ParamFoundInIndex_notRequested,1)-1)
    disp(IndexData.ParamFoundInIndex_notRequested)
    
    
    n_lines=size(RT_OUT_PARAM_split,1);
    for i=1:size(list_of_parameters_to_treat,1)
        
        i_param=list_of_parameters_to_treat(i);
        fprintf('Reading %s qc and mode information\n',i_param);
        
        %Initialise empty array
        IndexData.param.(i_param).qc=strings(n_lines,1);
        IndexData.param.(i_param).mode=strings(n_lines,1);
        
        %find the param position in the index file (both lines and position in
        %param field):
        i_pos_param=contains(RT_OUT_PARAM_split,i_param);
        
        %check if there is no anomaly:
        i_line_with_multiple_param=(sum(i_pos_param,2)>1);
        if sum(i_line_with_multiple_param)>0
            fprintf('Issue with param %s that has multiple name occurrences found in index lines \n',i_param)
            find(i_line_with_multiple_param)
        end
        
        
        % find the concerned lines:
        %/!\ the indexing function goes in the sense of Matlab array
        % reading sense (read by column then by line): The matrixes must be 
        % transposed and it only works with one parameter occurrence in the 
        % line (see above warning).
        i_pos_param_tr=i_pos_param';
        RT_OUT_param_profqc_tr=RT_OUT_param_profqc';
        RT_OUT_param_mode_tr=RT_OUT_param_mode';
        i_line_with_param=(sum(i_pos_param,2)==1);
        IndexData.param.(i_param).qc(i_line_with_param)    = RT_OUT_param_profqc_tr(i_pos_param_tr);
        IndexData.param.(i_param).mode(i_line_with_param)  = RT_OUT_param_mode_tr(i_pos_param_tr);
        IndexData.param.(i_param).nb_profile = sum(i_line_with_param);
        IndexData.param.(i_param).presence   = sum(i_pos_param,2);

    end

end

disp('---- List of unfound WMOs in the index file or without any ascending profile if conf i_descending_profile=0')
i_unfound_wmos=~ismember(float_ref_char,IndexData.profile_WMO);
disp(float_ref_char(i_unfound_wmos,:))
IndexData.WMOs_unfound_in_IndexFile=string(float_ref_char(i_unfound_wmos,:));

