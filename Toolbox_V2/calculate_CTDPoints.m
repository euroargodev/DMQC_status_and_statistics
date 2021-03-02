Floats_lifunction [Floats] = calculate_CTDPoints(Floats, config_param, drift_points)
% calculates CTD points
% EXAMPLE: [Floats] = calculate_CTDPoints(Floats, config_param, drift_points)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
% Floats: struc with floats data
% config_param: configuration parameters to be used for calculating CTD
%     points
% drift_points: 0 (CTD points mesured during drift are not taken in to
%     account) or 1 (CTD points mesured during drift are taken in to account)
%
% OUTPUT
% Floats: struc with recalculated cycle period in hours or days
%
% AUTHOR: Andrea Garcia Juan, Euro-Argo ERIC Office
%         (andrea.garcia.juan@euro-argo.eu)
%
% Modified on 2020-02-14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

n_floats = size(Floats.WMO,1);
n_fields = length(config_param);

for ifloat = 1:n_floats
        
    % format
    for ifield = 1:n_fields
        if isempty(Floats.(config_param{ifield}).data{ifloat})
            Floats.(config_param{ifield}).data{ifloat} = NaN(length(Floats.CONFIG_ProfilePressure_dbar.data{ifloat}),1);
        end
    end
        
    % calculation
    if drift_points == 1
    % with points during drift
        Floats.CONFIG_CTDPoints_NUMBER.data{ifloat} = Floats.CONFIG_PressureThresholdDataReductionShallowToIntermediate_dbar.data{ifloat}./Floats.CONFIG_ProfileSurfaceSlicesThickness_dbar.data{ifloat} + ...
          (Floats.CONFIG_PressureThresholdDataReductionIntermediateToDeep_dbar.data{ifloat} - Floats.CONFIG_PressureThresholdDataReductionShallowToIntermediate_dbar.data{ifloat})./Floats.CONFIG_ProfileIntermediateSlicesThickness_dbar.data{ifloat} + ...
          (Floats.CONFIG_ProfilePressure_dbar.data{ifloat} - Floats.CONFIG_PressureThresholdDataReductionIntermediateToDeep_dbar.data{ifloat})./Floats.CONFIG_ProfileBottomSlicesThickness_dbar.data{ifloat} + ...
          (Floats.CONFIG_CycleTime_hours.data{ifloat}./24).*(24./Floats.CONFIG_ParkSamplingPeriod_hours.data{ifloat});
        Floats.CONFIG_CTDPoints_NUMBER.cycle{ifloat} = Floats.CONFIG_CycleTime_hours.cycle{ifloat};
    elseif drift_points == 0
    % without points during drift
        Floats.CONFIG_CTDPoints_NUMBER.data{ifloat} = Floats.CONFIG_PressureThresholdDataReductionShallowToIntermediate_dbar.data{ifloat}./Floats.CONFIG_ProfileSurfaceSlicesThickness_dbar.data{ifloat} + ...
          (Floats.CONFIG_PressureThresholdDataReductionIntermediateToDeep_dbar.data{ifloat} - Floats.CONFIG_PressureThresholdDataReductionShallowToIntermediate_dbar.data{ifloat})./Floats.CONFIG_ProfileIntermediateSlicesThickness_dbar.data{ifloat} + ...
          (Floats.CONFIG_ProfilePressure_dbar.data{ifloat} - Floats.CONFIG_PressureThresholdDataReductionIntermediateToDeep_dbar.data{ifloat})./Floats.CONFIG_ProfileBottomSlicesThickness_dbar.data{ifloat};
        Floats.CONFIG_CTDPoints_NUMBER.cycle{ifloat} = Floats.CONFIG_CycleTime_hours.cycle{ifloat};
    else
            disp('drift_points should be 0 or 1')
    end
        
    
end