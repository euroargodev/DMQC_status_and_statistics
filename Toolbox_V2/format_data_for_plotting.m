function [Floats, notused_floats] = format_data_for_plotting(Floats, remove_cycles)
% format data for plotting: remove empty floats, all NaN floats and
% chosen cycles
% EXAMPLE: [Floats, notused_floats] = format_data_for_plotting(Floats, notused_floats, remove_cycles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
% Floats: struc with floats data
% notused_floats: struc from get_floats_data_gdac with floats which are not
%     going to be used 
% remove_cycles: cycles to be remove
%
% OUTPUT
% Floats: struc with formated data
% notused_floats: struc from with deleted floats
%
% AUTHOR: Andrea Garcia Juan, Euro-Argo ERIC
%         (andrea.garcia.juan@euro-argo.eu)
%
% Modified on 2020/03/12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


variables_list = fieldnames(Floats);
variables_list(contains(variables_list,'dac')) = [];
variables_list(contains(variables_list,'WMO')) = [];

notused_floats.WMO = cell(0,0);
notused_floats.why = cell(0,0);
notused_floats.dac = cell(0,0);
    
% remove empty floats
for ivar = 1:length(variables_list)
    index_empty = cellfun(@(x) isempty(x),Floats.(variables_list{ivar}).data);
    for ivar2 = 1:length(variables_list)
        Floats.(variables_list{ivar2}).data(index_empty) = [];
        Floats.(variables_list{ivar2}).cycle(index_empty) = [];
    end
    
    notused_floats.WMO(end+1:end + sum(index_empty)) = cellstr(Floats.WMO.data(index_empty,:))';
    notused_floats.why(end+1:end + sum(index_empty)) = cellstr(repmat('empty variable', sum(index_empty), 1));
    notused_floats.dac(end+1:end + sum(index_empty)) = Floats.dac.data(index_empty); 
    Floats.dac.data(index_empty) = [];
    Floats.WMO.data(index_empty,:) = [];
end

% remove all NaN floats
for ivar = 1:length(variables_list)
    index_nan = cellfun(@(x) all(isnan(x)),Floats.(variables_list{ivar}).data);
    for ivar2 = 1:length(variables_list)
        Floats.(variables_list{ivar2}).data(index_nan) = [];
        Floats.(variables_list{ivar2}).cycle(index_nan) = [];
    end
    
    notused_floats.WMO(end+1:end + sum(index_nan)) = cellstr(Floats.WMO.data(index_nan,:))';
    notused_floats.why(end+1:end + sum(index_nan)) = cellstr(repmat('all NaN', sum(index_nan), 1));
    notused_floats.dac(end+1:end + sum(index_nan)) = Floats.dac.data(index_nan); 
    Floats.dac.data(index_nan) = [];
    Floats.WMO.data(index_nan,:) = [];
end

n_floats = size(Floats.WMO.data,1);


% remove cycles
for ifloat = 1:n_floats
    
    % remove cycles
    for ivar = 1:length(variables_list)
       for irem = 1: length(remove_cycles)
           Floats.(variables_list{ivar}).data{ifloat}(Floats.(variables_list{ivar}).cycle{ifloat} == remove_cycles(irem)) = [];
           Floats.(variables_list{ivar}).cycle{ifloat}(Floats.(variables_list{ivar}).cycle{ifloat} == remove_cycles(irem)) = []; 
       end
   end
           
end


% remove empty floats
for ivar = 1:length(variables_list)
    index_empty = cellfun(@(x) isempty(x),Floats.(variables_list{ivar}).data);
    for ivar2 = 1:length(variables_list)
        Floats.(variables_list{ivar2}).data(index_empty) = [];
        Floats.(variables_list{ivar2}).cycle(index_empty) = [];
    end
    
    notused_floats.WMO(end+1:end + sum(index_empty)) = cellstr(Floats.WMO.data(index_empty,:))';
    notused_floats.why(end+1:end + sum(index_empty)) = cellstr(repmat('empty variable', sum(index_empty), 1));
    notused_floats.dac(end+1:end + sum(index_empty)) = Floats.dac.data(index_empty); 
    Floats.dac.data(index_empty) = [];
    Floats.WMO.data(index_empty,:) = [];
end

% remove all NaN floats
for ivar = 1:length(variables_list)
    index_nan = cellfun(@(x) all(isnan(x)),Floats.(variables_list{ivar}).data);
    for ivar2 = 1:length(variables_list)
        Floats.(variables_list{ivar2}).data(index_nan) = [];
        Floats.(variables_list{ivar2}).cycle(index_nan) = [];
    end
    
    notused_floats.WMO(end+1:end + sum(index_nan)) = cellstr(Floats.WMO.data(index_nan,:))';
    notused_floats.why(end+1:end + sum(index_nan)) = cellstr(repmat('all NaN', sum(index_nan), 1));
    notused_floats.dac(end+1:end + sum(index_nan)) = Floats.dac.data(index_nan); 
    Floats.dac.data(index_nan) = [];
    Floats.WMO.data(index_nan,:) = [];
end