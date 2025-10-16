
function [ pt_input_parameters ] = load_configuration( ps_filename )

%
% loadConfiguration
%
% loadConfiguration loads a set of properties from a file
% and returns the results in a structure variable. This routine 
% copied from OWC tool.
%
% The file is formatted as value pairs (one line per value):
%
%    key_name=key_value
%
% Also, any line starting with a '%' will be ignored.
%
% Note: all values of key_name must be valid matlab
%       variable names.
%
%       In the case of duplicate entries, the last key
%       to appear in the file will override previous
%       definitions.
%
%       Invalid lines are disregarded.
%
%  Jason Fabritz, 2001
%

lh_file = fopen( ps_filename ) ;

pt_input_parameters = struct('conf', ps_filename ) ;

while not(feof(lh_file)) 
    ls_line = fgetl(lh_file);
    ls_line(strfind(ls_line, ' ')) = [];
    if ~isempty(ls_line)
        if ls_line(1:1) ~= '%'
        ln_equals_index = strfind(ls_line, '=');
            if ~isempty(ln_equals_index)
                if ln_equals_index(1) > 1
                    sKey   = ls_line(1:ln_equals_index(1) - 1);
                    sValue = ls_line(ln_equals_index(1) + 1:length(ls_line));
                    pt_input_parameters.(sKey)=sValue;
                end
            end
        end
    end   
end

fclose( lh_file ) ;

return

