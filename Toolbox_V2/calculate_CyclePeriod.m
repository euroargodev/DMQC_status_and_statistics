function [Floats] = calculate_CyclePeriod(Floats,option)
% Recalculate cycle period in hours or days
% EXAMPLE: [Floats] = calculate_CyclePeriod(Floats,option)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
% Floats: struc with floats data
% option: 'days' or 'hours'
%
% OUTPUT
% Floats: struc with recalculated cycle period in hours or days
%
% AUTHOR: Andrea Garcia Juan, Euro-Argo ERIC Office
%         (andrea.garcia.juan@euro-argo.eu)
%
% Modified on 2020/02/14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_floats = size(Floats.WMO.data,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONVERSION TO DAYS %%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(option,'days')
% change to days
    for ifloat = 1:n_floats
    
        if isempty(Floats.CONFIG_CycleTime_days.data{ifloat})
            WMO = Floats.WMO.data{ifloat};
            fprintf("WMO : %s  => Cycle Time period not in days ! ", WMO);
            
            if isempty(Floats.CONFIG_CycleTime_hours.data{ifloat})
                disp('Conversion from minutes')
                Floats.CONFIG_CycleTime_days.data{ifloat} = Floats.CONFIG_CycleTime_minutes.data{ifloat}/1440;
                Floats.CONFIG_CycleTime_days.cycle{ifloat} = Floats.CONFIG_CycleTime_minutes.cycle{ifloat};
            else
                disp('Conversion from hours')
                Floats.CONFIG_CycleTime_days.data{ifloat} = Floats.CONFIG_CycleTime_hours.data{ifloat}/24;
                Floats.CONFIG_CycleTime_days.cycle{ifloat} = Floats.CONFIG_CycleTime_hours.cycle{ifloat};
            end
        end
    
    end
    % remove field hours
    Floats = rmfield(Floats,'CONFIG_CycleTime_hours');
    Floats = rmfield(Floats,'CONFIG_CycleTime_minutes');
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONVERSION TO HOURS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

elseif strcmp(option,'hours')
    % change to hours
    for ifloat = 1:n_floats
    
        if isempty(Floats.CONFIG_CycleTime_hours.data{ifloat})
            WMO = Floats.WMO.data{ifloat};
            fprintf("WMO : %s  => Cycle Time period not in hours ! ", WMO);
         
            if isempty(Floats.CONFIG_CycleTime_days.data{ifloat})
                disp('Conversion from minutes')
                Floats.CONFIG_CycleTime_hours.data{ifloat} = Floats.CONFIG_CycleTime_minutes.data{ifloat}/60;
                Floats.CONFIG_CycleTime_hours.cycle{ifloat} = Floats.CONFIG_CycleTime_minutes.cycle{ifloat};
            else
                disp('Conversion from days')
                Floats.CONFIG_CycleTime_hours.data{ifloat} = Floats.CONFIG_CycleTime_days.data{ifloat}*24;
                Floats.CONFIG_CycleTime_hours.cycle{ifloat} = Floats.CONFIG_CycleTime_days.cycle{ifloat};
            end
        end
    
    end
    % remove field days
    Floats = rmfield(Floats,'CONFIG_CycleTime_minutes');
    Floats = rmfield(Floats, 'CONFIG_CycleTime_days');

    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONVERSION TO MINUTES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 elseif strcmp(option,'minutes')
    % change to minutes
     for ifloat = 1:n_floats
         
         if isempty(Floats.CONFIG_CycleTime_minutes.data{ifloat})
             WMO = Floats.WMO.data{ifloat};
             fprintf("WMO : %s  => Cycle Time period not in minutes ! ", WMO);
             
             if isempty(Floats.CycleTime_days.data{ifloat})
                 disp('Conversion from hours')
                 Floats.CONFIG_CycleTime_minutes.data{ifloat} = Floats.CONFIG_CycleTime_hours.data{ifloat}*60;
                 Floats.CONFIG_CycleTime_minutes.cycle{ifloat} = Floats.CONFIG_CycleTime_hours_cycle{ifloat};
             else
                 disp('Conversion from days')
                 Floats.CONFIG_CycleTime_minutes.data{ifloat} = Floats.CONFIG_CycleTime_days.data{ifloat}*1440;
                 Floats.CONFIG_CycleTime_minutes.cycle{ifloat} = Floats.CONFIG_CycleTime_days_cycle{ifloat};
             end
         end
     end
     % Remove field hours
     Floats = rmfield(Floats, 'CONFIG_Cycle_Time_hours');
     Floats = rmfield(Floats, 'CONFIG_Cycle_Time_days');
    
    
else
    disp('option input should be days, hours or minutes')
end

   