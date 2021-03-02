function [FloatsPaths] = get_floats_filespath(gdac_dir, float_ref, output_file)
% EXEMPLE: [FloatsPaths] = get_floats_filespath(gdac_dir, float_ref, output_file)
% get_floats_filespath gets files paths from ar_index_global_meta.txt file 
% for given floats and optionally creates a .txt file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% gdac_dir: gdac path
% float_ref: cell array with float WMO (exemple: {'3901875','3901872'})
% output_file: (optional) output file name
%
% Output
% FloatsPaths: struct containing WMO, DAC, float path and index update date
% txt file (optional): three columns: WMO, RT and files_path
%
% NOTES
% (1) Results in FloatsPaths struct are cell values
%
% AUTHOR: Andrea Garcia Juan, Euro-Argo ERIC
%         (andrea.garcia.juan@euro-argo.eu)
%
% Modified on 2020-02-25 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% float_ref = {'3901872','3901875'} % VERY IMPORTANT: cell array
% gdac_dir = '/home/ref-argo/gdac'
% output_file = 'testfile.txt'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% open file
%tic
disp('Looking for floats files paths in index file ...')
% open file
index_file_dir = [gdac_dir '/ar_index_global_meta.txt'];
fid = fopen(index_file_dir);

% read file line by line
tline = fgetl(fid);

% reading header
while tline(1) == '#'
    matches = contains(tline,'Date of update');
    if all(matches)~=0
        FloatsPaths.index_update = tline(20:end);
    end
    tline = fgetl(fid);
end


% get variables names
variables = strsplit(tline,',');
% count number of columns
n_var = length(variables);
% format string for textscan
string_format = repmat('%s ',1,n_var);

% get data
tline = fgetl(fid);
%initialisations
ifloat = 0;
while ischar(tline) % lines loop

        % Find float ref in line
        C = textscan(tline,string_format,'Delimiter',',');
        matches = contains(C{1},float_ref);
        
        if all(matches)~=0
            % counter
            ifloat= ifloat+1;
            
            % read line
            % disp(tline)
            % C = textscan(tline,string_format,'Delimiter',',');
                        
            % save in struc
            % file path
            file_path = char(C{1});
            file_path = file_path(1:end-15);
            FloatsPaths.files_path{ifloat} = file_path;
            % found float reference
            p = textscan(char(C{1}),'%s %s %s %s','Delimiter','/');
            FloatsPaths.WMO{ifloat} = char(p{2});
            % DAC
            FloatsPaths.DAC{ifloat} = char(p{1});
            
        end
     
    tline = fgetl(fid);
end

fclose(fid);

% Warning not found floats
not_found = float_ref(~ismember(float_ref,FloatsPaths.WMO));
if ~isempty(not_found)
    fprintf(2,'Floats not found: %s \n',strjoin(cellstr(not_found),', '))
end


disp('Index file closed')

%% output txt file

if nargin == 3 % creates txt file
    
    fprintf('\nSaving results in %s ...\n',output_file)
    fid=fopen(output_file,'w');

    % header
    header = 'WMO,DAC,files_path';
    fprintf(fid, '%s \n', header);
    % data
    table = [FloatsPaths.WMO; FloatsPaths.DAC; FloatsPaths.files_path];
    fprintf(fid, '%s,%s,%s\n',table{:});
    
    fclose(fid);
    fprintf('Results saved\n')
    
end
%toc