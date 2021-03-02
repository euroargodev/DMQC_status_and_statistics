function [vertical_km, vertical_km_mean, floats_age, last_cycle] = get_verticalkm_age_multiprof(Floats, dac_dir) 
% EXAMPLE: [vertical_km, vertical_km_mean, floats_age, last_cycle] = get_verticalkm_age_multiprof(Floats, dac_dir)
% calculates vertical km, number of cycles and float age using multi profile files
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
% floats_age: float age in years from cycle 0
% last_cycle: number of cycles performed by float
%
% NOTES:
% (1) Vertical km are calculated using PRES variable in multiprofile file
%
% AUTHOR: Andrea Garcia Juan, Euro-Argo ERIC
%         (andrea.garcia.juan@euro-argo.eu)
%
% Modified on 2020/03/20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ')
disp('Calculating vertical km and floats age...')

n_floats = length(Floats.WMO);

vertical_km = NaN(n_floats,1);
vertical_km_mean = NaN(n_floats,1);
floats_age = NaN(n_floats,1);
last_cycle = NaN(n_floats,1);

%% floats loop
for ifloat = 1: n_floats
    
    disp(' ')
    fprintf('%s\n', Floats.WMO{ifloat})
    % floats directory path string
    dac = Floats.DAC{ifloat};
    % profiles file dir
    prof_file = [dac_dir '/dac/' char(dac) '/' Floats.WMO{ifloat} '/' Floats.WMO{ifloat} '_prof.nc'];
    
    % get pres matrix from prof file
    try
        pres = ncread(prof_file, 'PRES');
        juld = ncread(prof_file, 'JULD');
        cycle_number = ncread(prof_file, 'CYCLE_NUMBER');
    catch e
        fprintf(2,'        %s\n', e.message)
        pres = NaN;
        juld = NaN;
        cycle_number = NaN;
    end
   
    vertical_m_vector = 2* max(pres);
    
    % sum all cycles
    vertical_km(ifloat) = sum(vertical_m_vector,'omitnan')/1000; % km
    vertical_km_mean(ifloat) = sum(vertical_m_vector,'omitnan')/length(vertical_m_vector); % km
    floats_age(ifloat) = (max(juld) - min(juld))/365;
    last_cycle(ifloat) = cycle_number(end);
    
end
