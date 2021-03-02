function [bar_matrix, analysis] = get_matrix_barplot(x_data, meta_data, units, option1, merge_little_values)
% generates a matrix for using in plot. config parameter is clasified with
% meta data
% option1: change or values or cycles


if strcmp(option1, 'change')
    
    CONFIG_values = cellfun(@(x) length(unique(x(~isnan(x)))), x_data);
    all_diff_values = unique(CONFIG_values);
    analysis.floats = zeros(1,length(all_diff_values));
    analysis.labels = cell(1,length(all_diff_values));
    for ivalue = 1: length(all_diff_values)
        analysis.floats(ivalue) = sum(CONFIG_values == all_diff_values(ivalue));
        analysis.meta(ivalue) = {meta_data(CONFIG_values == all_diff_values(ivalue),:)};
        analysis.labels(ivalue) = {[num2str(all_diff_values(ivalue)-1)]};
        %' ' units
    end
    
    % bar matrix
    all_meta = unique(cellstr(meta_data));
    for ivalue = 1: length(all_diff_values)
        for imeta = 1: length(all_meta)
            bar_matrix(imeta, ivalue) = sum(strcmp(cellstr(analysis.meta{ivalue}), all_meta(imeta)));
        end
    end

    
elseif strcmp(option1, 'values') 
    
    CONFIG_values = cellfun(@(x) x(1), x_data); % first value
    all_diff_values = sort(unique(cat(1,x_data{:})));
    all_diff_values(isnan(all_diff_values)) = [];

    for ivalue = 1: length(all_diff_values)
        analysis.floats(ivalue) = sum(CONFIG_values == all_diff_values(ivalue));
        analysis.meta(ivalue) = {meta_data(CONFIG_values == all_diff_values(ivalue),:)};
        analysis.labels(ivalue) = {[num2str(round(all_diff_values(ivalue),2)) ' ' units]};
    end
    
    % bar matrix
    all_meta = unique(cellstr(meta_data));
    for ivalue = 1: length(all_diff_values)
        for imeta = 1: length(all_meta)
            bar_matrix(imeta, ivalue) = sum(strcmp(cellstr(analysis.meta{ivalue}), all_meta(imeta)));
        end
    end
    
    
elseif strcmp(option1, 'cycles')
    
    CONFIG_values = cat(1,x_data{:}); % all values
    all_diff_values = sort(unique(CONFIG_values));
    all_diff_values(isnan(all_diff_values)) = [];
    
    analysis.cycles = NaN(length(x_data), length(all_diff_values));
    analysis.labels = cell(1,length(all_diff_values));
    for ivalue = 1: length(all_diff_values)
        %all_diff_values(ivalue)
        analysis.cycles(:,ivalue) = cellfun(@(x) sum(x == all_diff_values(ivalue)), x_data);
        %analysis.cycles(:,ivalue)
        analysis.meta(ivalue) = {meta_data(analysis.cycles(:,ivalue) ~= 0, :)};
        %meta_data(analysis.cycles(:,ivalue) ~= 0, :)
        analysis.labels(ivalue) = {[num2str(round(all_diff_values(ivalue),2)) ' ' units]};
    end
    
    
    % bar matrix
    all_meta = unique(cellstr(meta_data));
    for ivalue = 1: length(all_diff_values)
        %disp(' ')
        %all_diff_values(ivalue)
        for imeta = 1: length(all_meta)
            %all_meta(imeta)
            %strcmp(cellstr(meta_data), all_meta(imeta))'
            %analysis.cycles(strcmp(cellstr(meta_data), all_meta(imeta)),ivalue)'
            
            bar_matrix(imeta, ivalue) = sum(analysis.cycles(strcmp(cellstr(meta_data), all_meta(imeta)),ivalue));
        end
    end
    
    analysis.cycles = sum(analysis.cycles);
    
else 
    
    fprintf(1, 'Option1 should be "change", "values" or "cycles"')

end


%% merge little values

%little values in one bar
if merge_little_values == 1 && (strcmp(option1, 'change') || strcmp(option1, 'values'))
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    nfloats_limit = 2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    little_values_index = analysis.floats <= nfloats_limit;
    start_pos = strfind([0,little_values_index == 1],[0 1]);
    end_pos = strfind([little_values_index==1,0],[1 0]);
    n_groups = length(start_pos);
    for igroup = 1:n_groups
       sum_value(igroup) = sum(analysis.floats(start_pos(igroup):end_pos(igroup)));
       %Copy sum in last position
       analysis.floats(end_pos(igroup)) = sum_value(igroup);
       if strcmp(option1, 'change')
           analysis.labels(end_pos(igroup)) = {[num2str(round(all_diff_values(start_pos(igroup))-1,2)) '-' num2str(round(all_diff_values(end_pos(igroup))-1,2))]};
           % ' ' units
       else
           analysis.labels(end_pos(igroup)) = {[num2str(round(all_diff_values(start_pos(igroup)),2)) '-' num2str(round(all_diff_values(end_pos(igroup)),2))]};
            %' ' units
       end
       bar_matrix(:,end_pos(igroup)) = sum(bar_matrix(:,start_pos(igroup):end_pos(igroup)),2);
    end
    % erase all except last one
    little_values_index(end_pos) = 0;
    analysis.floats(little_values_index) = [];
    all_diff_values(little_values_index) = [];
    analysis.labels(little_values_index) = [];
    bar_matrix(:,little_values_index) = []; 
    
end

if merge_little_values == 1 && strcmp(option1, 'cycles')
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    ncycles_limit = 400;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    little_values_index = analysis.cycles <= ncycles_limit;
    start_pos = strfind([0,little_values_index == 1],[0 1]);
    end_pos = strfind([little_values_index==1,0],[1 0]);
    n_groups = length(start_pos);
    for igroup = 1:n_groups
       sum_value(igroup) = sum(analysis.cycles(start_pos(igroup):end_pos(igroup)));
       %Copy sum in last position 
       analysis.cycles(end_pos(igroup)) = sum_value(igroup);
       analysis.labels(end_pos(igroup)) = {[num2str(round(all_diff_values(start_pos(igroup)),2)) '-' num2str(round(all_diff_values(end_pos(igroup)),2)) ' ' units]};
       bar_matrix(:,end_pos(igroup)) = sum(bar_matrix(:,start_pos(igroup):end_pos(igroup)),2);
    end
    % erase all except last one
    little_values_index(end_pos) = 0;
    analysis.cycles(little_values_index) = [];
    all_diff_values(little_values_index) = [];
    analysis.labels(little_values_index) = [];
    bar_matrix(:,little_values_index) = []; 
    
end
% output

analysis.all_diff_values = all_diff_values;
analysis.all_meta = all_meta;
