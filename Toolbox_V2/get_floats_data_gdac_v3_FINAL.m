function [Floats, WMO2delete_index] = get_floats_data_gdac_v3_FINAL(Floats_list, variables_list, where_file, gdac_path, mc_code)
% EXAMPLE: [Floats] = get_floats_data_gdac_v3(Floats_list, variables_list, where_file, gdac_path, mc_code)
% gets data from tech, traj, aux or meta (config) file
% Difference with get_floats_data_gdac : no format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
% Floats_list: struc with at least two fields:
%              Floats_list.WMO: list of floats reference number
%              Floats_list.RT or Floats_list.DAC: name of dac folder in
%              gdac or dac name in database 
% variables_list: list of variables we want to get
% where_file: where to find each variable in variables_list. Only six 
%             options: traj, tech, aux, config, meta and index
% gdac_path: path to gdac folders
% mc_code: MEASUREMENT_CODE for traj variable. If no mc code, we get all
% mc_codes
%
% OUTPUTS
% Floats: struct where fields are variables_list, WMO, dac and subfields 
%         are data and cycle. Each field contains n_floats cells and each 
%         cell a different number of cycles. Not found floats are not 
%         included in the struct. 
%
% Auxiliary functions:
%    database2dac
%    get_traj_param
%    get_tech_param
%    get_configparam_meta
%    get_data_from_index
%
% NOTES
% (1) Data is formated not formated dont each variable could have different
% number of cycles for each float.
% (2) If JULD or JULD_ADJUSTED is in variables_list, REFERENCE_DATE_TIME is
% also in output struct
% (3) Trajectory variables are searched in delayed mode file ('_Dtraj.nc').
% If D file does not exit, real time file ('_Rtraj.nc') is used
% (5) Index variables are searched in ar_index_global_prof.txt file. This
% research can last some minutes
%
% AUTHOR: Andrea Garcia Juan, Euro-Argo ERIC
%         (andrea.garcia.juan@euro-argo.eu)
%
% Modified on 2020/03/13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% TODO: check if mc_code when traj files
if nargin < 5
    mc_code = [];
end


% variables to look for
tech_variables = variables_list(contains(where_file,'tech'));
aux_variables = variables_list(contains(where_file,'aux'));
traj_variables  = variables_list(strcmp(where_file,'traj'));
% check length mc_code
if ~isempty(traj_variables) && length(mc_code)~= length(traj_variables)
    error('mc_code should be same length as traj variables')
end
trajc_variables = variables_list(contains(where_file,'trajc'));
config_variables = variables_list(contains(where_file,'config'));
index_variables = variables_list(contains(where_file,'index'));
meta_variables = variables_list(contains(where_file,'meta'));
% not index variables in variables list (out of for loop)
variables_list(contains(where_file,'index')) = [];
where_file(contains(where_file,'index')) = [];
% paths
dac_path = [gdac_path 'dac/'];
aux_path = [gdac_path '/aux/'];


% initialisations
itotal = length(Floats_list.WMO);
Floats.WMO.data = cellstr(Floats_list.WMO);
Floats.dac.data = Floats_list.DAC;
ifloat = 1;
n_floats = length(Floats.WMO.data);



while ifloat <= n_floats % floats loop
    
    disp(' ')
    disp(Floats.WMO.data{ifloat})
    
    
    %% TRAJ file
    % check if we need variables from traj file
    if ~isempty(traj_variables)
        
        traj_file = [dac_path char(Floats.dac.data{ifloat}) '/' Floats.WMO.data{ifloat} '/' Floats.WMO.data{ifloat} '_Dtraj.nc'];
        
        % check if file exits
        if exist(traj_file, 'file') ~= 2  % file dont exist
            % try RT file
             traj_file = [dac_path char(Floats.dac.data{ifloat}) '/' Floats.WMO.data{ifloat} '/' Floats.WMO.data{ifloat} '_Rtraj.nc'];
            if exist(traj_file, 'file') ~= 2  % file dont exist
                disp(ifloat)
                WMO2delete_index(ifloat) = ifloat;
                %test3(ifloat) = Floats.WMO(ifloat);
                fprintf(2,'          traj file not found\n')
            else
                WMO2delete_index(ifloat) = nan;
            end
        end 
    
        % get traj data
        try
            [Data.traj] = get_traj_param(traj_file, traj_variables, mc_code); % 703: surface locations during surface drift
                 
            if contains('JULD',variables_list) || contains('JULD_ADJUSTED',variables_list)
                Floats.REFERENCE_DATE_TIME.data{ifloat} = ncread(traj_file,'REFERENCE_DATE_TIME');
            end
            
        catch e
            fprintf(2,'          %s\n', e.message)
            for ivar = 1:length(traj_variables)
                Data.traj.(traj_variables{ivar}).data = [];
                Data.traj.(traj_variables{ivar}).cycle = [];
                Data.traj.(traj_variables{ivar}).mc_code = [];
            end
        end
             
        
    end % end traj file
    
    
    %% TRAJ file (for n_cycle variables)
    % check if we need variables from traj file
    if ~isempty(trajc_variables)
        
        traj_file = [dac_path char(Floats.dac.data{ifloat}) '/' Floats.WMO.data{ifloat} '/' Floats.WMO.data{ifloat} '_Dtraj.nc'];
        
        % check if file exits
        if exist(traj_file, 'file') ~= 2  % file dont exist
            % try RT file
             traj_file = [dac_path char(Floats.dac.data{ifloat}) '/' Floats.WMO.data{ifloat} '/' Floats.WMO.data{ifloat} '_Rtraj.nc'];
            if exist(traj_file, 'file') ~= 2  % file dont exist
                
                fprintf(2,'          traj file not found\n')
            end
        end 
    
        % get trajc data
        for ivar = 1:length(trajc_variables)
            try
                [Data.trajc.(trajc_variables{ivar}).data] = ncread(traj_file,char(trajc_variables{ivar}));
                [Data.trajc.(trajc_variables{ivar}).cycle] = ncread(traj_file,'CYCLE_NUMBER_INDEX');
            catch e
                fprintf(2,'       %s\n', e.message)
                Data.trajc.(trajc_variables{ivar}).data = [];
                Data.trajc.(trajc_variables{ivar}).cycle = [];
            end
        end
                
    end % end traj file
    
    
    %% TECH file
    % check if we need variables from tech file
    if ~isempty(tech_variables)

        tech_file = [dac_path char(Floats.dac.data{ifloat}) '/' Floats.WMO.data{ifloat} '/' Floats.WMO.data{ifloat} '_tech.nc'];
        
        % check if file exist
        if exist(tech_file, 'file') ~= 2  % file dont exist
            fprintf(2,'          tech file not found\n')
        end
        
        % get data
        try
            [Data.tech] = get_tech_param(tech_file,tech_variables);
        catch e
            fprintf(2,'        %s\n', e.message)
            for ivar = 1:length(tech_variables)
                Data.tech.(tech_variables{ivar}).data = [];
                Data.tech.(tech_variables{ivar}).cycle = [];
            end
        end
        
 
    end % tech variables
    
    
    %% AUX file
     if ~isempty(aux_variables)
        
        aux_file = [aux_path char(Floats.dac.data{ifloat}) '/' Floats.WMO.data{ifloat} '/' Floats.WMO.data{ifloat} '_tech_aux.nc'];
        
        % check if file exist
        if exist(aux_file, 'file') ~= 2  % file dont exist
            fprintf(2,'          aux file not found\n')
        end
        
        % get data
        try
            [Data.aux] = get_tech_param(aux_file,aux_variables);
        catch e
            fprintf(2,'        %s\n', e.message)
            for ivar = 1:length(aux_variables)
                Data.aux.(aux_variables{ivar}).data = [];
                Data.aux.(aux_variables{ivar}).cycle = [];
            end
        end
        
        
    end % aux variables
    
    
    %% CONFIG variable
    
    if ~isempty(config_variables)
        
        config_file = [dac_path char(Floats.dac.data{ifloat}) '/' Floats.WMO.data{ifloat} '/' Floats.WMO.data{ifloat} '_meta.nc'];
        
        % check if file exist
        if exist(config_file, 'file') ~= 2  % file dont exist
            fprintf(2,'          meta file not found\n')
        end
        
        % get data
        try
            [Data.config] = get_configparam_meta(config_file, config_variables);          
        catch e
            % tech param do not exits
            fprintf(2,'        %s\n', e.message)
            for ivar = 1:length(config_variables)
                Data.config.(config_variables{ivar}).data = [];
                Data.config.(config_variables{ivar}).cycle = [];
                Data.config.(config_variables{ivar}).mission = [];
            end
        end
        
    end % config variables

    
    %% META variable
    
    if ~isempty(meta_variables)
        
        meta_file = [dac_path char(Floats.dac.data{ifloat}) '/' Floats.WMO.data{ifloat} '/' Floats.WMO.data{ifloat} '_meta.nc'];
        
        % check if file exist
        if exist(config_file, 'file') ~= 2  % file dont exist
            fprintf(2,'          meta file not found\n')
        end
        
        % get data
        for ivar = 1: length(meta_variables)
            try
                Data.meta.(meta_variables{ivar}).data = ncread(meta_file, meta_variables{ivar});
            catch e
                % meta param do not exits
                fprintf(2,'        %s\n', e.message)
                Data.meta.(meta_variables{ivar}).data = [];
            end

        end
        
      
    end % meta variables
   
    
    %% creating cells
   
    for ivar = 1:length(variables_list)
        variables_list{ivar} = erase(variables_list{ivar},'^');
        variables_list{ivar} = erase(variables_list{ivar},'/');
        
        Floats.(variables_list{ivar}).data{ifloat} = Data.(where_file{ivar}).(variables_list{ivar}).data;
        Floats.(variables_list{ivar}).cycle{ifloat} = Data.(where_file{ivar}).(variables_list{ivar}).cycle;
        
        if strcmp(where_file{ivar},'traj')
            Floats.(variables_list{ivar}).m_code{ifloat} = Data.(where_file{ivar}).(variables_list{ivar}).m_code;
        elseif strcmp(where_file{ivar},'config')
            Floats.(variables_list{ivar}).mission{ifloat} = Data.(where_file{ivar}).(variables_list{ivar}).mission;
        end
        
    end
    
        
    ifloat = ifloat +1;

    
end % floats loop

%% INDEX variables
    
if ~isempty(index_variables)
    
    disp(' ')
    index_file_dir = [ gdac_path '/ar_index_global_prof.txt']; 
    [Data.index] = get_data_from_index(index_file_dir, cellstr(Floats.WMO.data), index_variables, 0);
        
    % format
    for ifloat = 1:size(Data.index.WMO.data,1)
        
        % create strct fields
        for ivar = 1:length(index_variables)
            Floats.(index_variables{ivar}).data{ifloat} = Data.index.(index_variables{ivar}).data{ifloat};
        end
    
    end
        
end
   
        
    
