function [IndexData] = get_data_from_index(index_file_dir, float_ref,Dprof)
% EXEMPLE: [IndexData] = get_data_from_index(index_file_dir, float_ref)
% get_data_from_index gets all variables from index file for given floats 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% index_file_dir: index file path
% float_ref: cell array with float WMO (exemple: {'3901875','3901872'})
% Dprof: (optional) if 1 descent profiles are consider, if 0 descent profiles
% are not consider. Default value: 1.
%
% Output
% IndexData: struct with all variables in file as fields. Values are char
% (fill value : '---') or double (fill value : NaN)
%
% NOTES
% (1) Data is stored in data fields (ex: IndexData.latitude.data) and 
% dimensions in dim field (ex: IndexData.latitude.dim). For some variables 
% there is a comment field (ex: IndexData.n_obs.comment)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TODO get more variables from file path
% TODO: if float 1 is not found and if other than date are 14 string size
% TODO : problems with profile_type as num

% INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%float_ref = {'3901872','3901875'} % VERY IMPORTANT: cell array
%index_file_dir = 'argo_profile_detailled_index.txt'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% no descent profiles defaut value
if nargin == 2
    Dprof = 1;
end

% open file
%tic
disp('Reading index file ...')
% open file
fid = fopen(index_file_dir);
% float reference
IndexData.WMO.data = char(float_ref);
IndexData.WMO.dim = "n_float x string7";
[a,b] = size(float_ref);
n_floats = max([a,b]);

% read file line by line
tline = fgetl(fid);

% reading header
while tline(1) == '#'
    matches = contains(tline,'Date of update');
    if all(matches)~=0
        IndexData.index_update.data = tline(20:end);
    end
    tline = fgetl(fid);
end
IndexData.index_update.comment = 'Date of index file last update';

% get variables names
variables = strsplit(tline,',');
% count number of columns
n_var = length(variables);
% format string for textscan
string_format = repmat('%s ',1,n_var);

% get data
tline = fgetl(fid);
%initialisations
n_obs = zeros(1,n_floats);
i=0;
while ischar(tline) % lines loop

        % Find float ref in line
        matches = contains(tline,float_ref);
        
        if all(matches)~=0
            % read line
%             disp(tline)
            C = textscan(tline,string_format,'Delimiter',',');
            
            % if descent profile
            if (contains(C{1},'D.nc') && Dprof == 0)
                tline = fgetl(fid);
                continue
            end 
            
            
            % discriminate diferent floats
            p = textscan(char(C{1}),'%s %s %s %s','Delimiter','/');
            float_str = char(p{2});
            i_new = find(contains(float_ref,float_str)); % i columns (floats name)
            if i_new ~= i 
                %disp("diferent float")
                k=0; % initialise rows (cycle number)
                i = i_new;
            end
                
            k=k+1; % k rows (cycle number)
             
            % save in struc
            for ivar=1:n_var
                IndexData.(variables{ivar}).data(k,1:length(char(C{ivar})),i) = char(C{ivar});  
            end
            
            % get more variables from file path (eventually add more)
            file_name = strsplit(char(p{4}),'/');
            cycle = strsplit(char(file_name),'_');
            cycle = char(erase(cycle(2),'.nc'));
            IndexData.cycle_number.data(k,1:length(cycle),i) = cycle;
            
            % number of observations for each float
            n_obs(i) = k;
            
        end
     
    tline = fgetl(fid);
end

fclose(fid);

% Warning not found floats
inot = find(~n_obs);
not_found = IndexData.WMO.data(inot,:);
if ~isempty(not_found)
    fprintf(2,'Floats not found: %s \n',strjoin(cellstr(not_found),', '))
end



% FORMATING
% add new variables from file path 
variables{end+1}= 'cycle_number';
% recalculate number of variables
n_var = length(variables);


% check if variable is a number or not
% TODO: if float 1 is not found and if other than date are 14 string size
% TODO : problems with profile_type as num
variable_num = zeros(1,n_var);
for ivar = 1:n_var
    
    var_value = IndexData.(variables{ivar}).data(2,:,1);
    var_value(find(var_value <= ' ' | isspace(var_value), 1)) = [];
    % 0 if char values, 1 if numerical values (dates are considered char)
    variable_num(ivar) = ~(length(var_value) == 14 || isempty(str2num(var_value)));

end


% fill with '---' or NaN when different number of cycles and dimensions
max_cycle = max(n_obs);
for ivar=1:n_var % variables loop
    
    if variable_num(ivar) == 1 % varibles is a number
        
        IndexData.(variables{ivar}).data = str2double(string(IndexData.(variables{ivar}).data));
        IndexData.(variables{ivar}).dim = "n_cycle x n_float"; % struct dimensions
        
    else % varibles is a char array
        
        string_length = length(IndexData.(variables{ivar}).data(2,:,1));
        for i = 1:n_floats % floats loop
            
            in_range = IndexData.(variables{ivar}).data(1:n_obs(i),:,i);
            repmat('-',max_cycle-n_obs(i),string_length);
            IndexData.(variables{ivar}).data(:,:,i) = [in_range ; repmat('-',max_cycle-n_obs(i),string_length)];
            IndexData.(variables{ivar}).dim = "n_cycle x string_size x n_float"; % struct dimensions
            
        end
   
    end

end

% number of observations field
IndexData.n_obs.data = n_obs;
IndexData.n_obs.dim = "n_floats";
IndexData.n_obs.comment = 'number of observations per float';
if Dprof == 0
    IndexData.n_obs.comment=[IndexData.n_obs.comment,' including D profiles'];
end


disp('Index file closed')
%toc