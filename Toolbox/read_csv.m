function [Metadata] = read_csv(file_dir,delimiter)
% EXEMPLE: [Metadata] = read_csv(file_dir,delimiter)
% read_csv read a csv file and generates an struct with file variables
% (values are in char)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% file_dir: csv file path
% delimiter : delimiter use in file (string). Exemples: ',',';'
%
% Output
% Metadata: struct with char values
%
% NOTES
% (1) csv file format: first line with variables, next lines with values
% (2) If the file contains special characters (é,à, etc) you should delate  
% them (Matlab don't understant special characters)
% (3) Be careful with csv from sql: some times for one fields there are more
% than one value separated by comma, for exemple: for PR_EXPERIMENT_ID field
% one value is MOCCA,MOCCA-NETHERLANDS. In such cases it is better to export
% the table in txt format with delimiter ' '.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TODO find delimiter automatically

% INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%file_dir = '/home1/datahome/co_dev/agarciaj/MOCCA_list.csv'
%delimiter = '",' % if " is use it should be in delimiter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% open file
fprintf('\nReading %s ...\n',file_dir)
fid = fopen(file_dir);

% read first line
tline = fgetl(fid);
tline = strrep(tline, ';', delimiter); % sometimes delimiters are not coherent in all file
% split line
variables = strsplit(tline,delimiter);
variables = erase(variables,'"'); % in the case there are "
%count number of columns
n_var = length(variables);

% delate blanks in variables names
for i = 1:n_var % erase blanks
    char_var = char(variables{i});
    char_var =  char_var(~isspace(char_var));
    char_var =  char_var(isletter(char_var));
    variables{i} = char_var;
end

tline = fgetl(fid);
tline = strrep(tline, ';', delimiter);
cnt=0;
while ischar(tline) % lines loop
    
    tline = strrep(tline, ';', delimiter);
    % counter
    cnt = cnt+1;
    C = regexp(tline,delimiter,'split');
    C = erase(C,'"'); % in the case there are "
    for i=1:n_var % variables loop
        Metadata.(variables{i})(cnt,1:length(C{i})) = C{i};
    end
    
    % next line
    tline = fgetl(fid);
    
end
 
fclose(fid);

disp('file closed')
