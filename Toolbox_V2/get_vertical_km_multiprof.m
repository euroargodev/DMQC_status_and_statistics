function [vertical_km, vertical_km_mean] = get_vertical_km_multiprof(Floats, dac_dir) 
% EXAMPLE: [vertical_km, vertical_km_mean] = get_vertical_km_multiprof(Floats, dac_dir)
% calculates vertical km using multi profile files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
% Floats: struct with al least two fields: WMO.data and dac.data 
% dac_dir: path to gdac 
%
% OUTPUT
% vertical_km: sum of vertical km (up and down) travelled by one float 
%      during all its life time
% vertical_km_mean: mean vertical km (up and down) travelled by a float 
%      in one cycle
%
% NOTES:
% (1) Vertical km are calculated using PRES variable in multiprofile file
%
% AUTHOR: Andrea Garcia Juan, Euro-Argo ERIC
%         (andrea.garcia.juan@euro-argo.eu)
%
% Modified on 2020/03/13 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp('Calculating vertical km ...')

n_floats = length(Floats.WMO.data);

vertical_km = NaN(n_floats,1);
vertical_km_mean = NaN(n_floats,1);
n_cycle = NaN(n_floats,1);

%% floats loop
for ifloat = 1: n_floats
    
    disp(' ')
    fprintf('%s\n', Floats.WMO.data{ifloat})
    % floats directory path string
    dac = Floats.dac.data{ifloat};
    % profiles file dir
    prof_file = [dac_dir '/dac/' char(dac) '/' Floats.WMO.data{ifloat} '/' Floats.WMO.data{ifloat} '_prof.nc'];
    
    % get pres matrix from prof file
    try
        pres = ncread(prof_file, 'PRES');
    catch e
        fprintf(2,'        %s\n', e.message)
        prof_file = NaN;
    end
   
    vertical_m_vector = 2* max(pres);
    
    % sum all cycles
    vertical_km(ifloat) = sum(vertical_m_vector,'omitnan')/1000; % km
    vertical_km_mean(ifloat) = sum(vertical_m_vector,'omitnan')/length(vertical_m_vector); % km
    
end
