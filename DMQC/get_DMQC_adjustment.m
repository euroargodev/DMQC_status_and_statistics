% get_DMQC_adjustment
% shows var adjustement for a list of floats
%
% Input (below)
% - Profiles files from dac (same file structure as dac)
% - floats list: csv file, with in header WMO, RT (dac in charge of real time
% processing) and DM (institution in charge of delayed mode processing), 
% separated by ';'. Exemple: floats_rt_dm.csv
% Paths to these files should be modified below
%
% Output
% -Figures: saved in (today date) folder (ex: 20181210). Look at format 
% options below.
% -DATA struc with results (it can be saved)
%
% Auxiliary functions needed:
%    read_csv
%    split_figure
%    export_fig package
%    legendflex
% 
% NOTES :
% (1) Getting data process may take some minutes (depending on number of floats)
% because we are checking all profile files
% (2) Descent profiles are not included
% (3) Format options (number of floats per figure and yaxis ticks size)
% can be modified below
%
% Modified on 20191126

close all                  
clear variables

% add paths (packages and auxiliary functions)
% aux_functions_path = [pwd '/aux_functions'];
% % aux_functions_path = '/home1/datahome/co_arg/agarciaj/DMQC_status/aux_functions';
% addpath(genpath(aux_functions_path))

% addpath /home1/datahome/co_arg/rcancoue/decodeur_matlab/work_Romain
addpath /home1/datahome/co_arg/larduini/Scripts/Toolbox
addpath /home1/datahome/co_arg/larduini/Scripts/Toolbox/flexLegend/legendflex
addpath /home1/datahome/co_arg/larduini/Scripts/Toolbox/flexLegend/setgetpos_V1.2
addpath /home1/datahome/co_arg/larduini/Scripts/Toolbox/export_fig-master % export a matlab figure



% INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
project_name = 'DMQC_adjustments_FSD_CTD'; % for png files names

floats_file_dir = '/home1/datahome/co_arg/larduini/Lists/European_CTD_FSD_SN_for_DMQC_status.csv'
dac_dir = '/home/ref-argo/gdac/dac'

variable='PSAL' % Variable to retrieve in the netcdf file. Examples: "PSAL" ; "TEMP"; "DOXY"
units = 'PSU' %other examples: °C or PSAL. Only used for figures titles and names

% -b file pattern
var_bio = 0;

update_date = datestr(now(),'yyyy-mm-dd'); % dac update date (if it is not 
% today because we are working with a snapshot, change it)
export_dir = '/home1/datahome/co_arg/larduini/Exports/DMQC/DMQC_adjustements'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FORMAT OPTIONS
% if so many floats, figures are divided. Number of floats in each figure 
% and yticks size can be change here:
floats_per_fig = 40; % number of floats per figure (80 is the very limit)
fontsize_wmo = 10; % yaxis ticks font size (8 is the limit)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TODO: VAR IMPUT: PRES var AND temp
% TODO : dimensions struc

% read floats floats list
[floats_list] = read_csv(floats_file_dir,';');
%floats_list.WMO = char(sort(cellstr(floats_list.WMO)));
DATA.WMO = floats_list.WMO;
n_floats = size(floats_list.WMO,1);

disp(' ')
disp('Getting data from profiles files ...')
% initialisations
n_cycle=NaN(n_floats,1);
for ifloat = 1:n_floats % floats loop
    
    fprintf('        %s\n', floats_list.WMO(ifloat,:))
    % floats directory path string
    dac = lower(cellstr(floats_list.RT(ifloat,:)));
    float_dir = [dac_dir,'/',char(dac),'/',floats_list.WMO(ifloat,:),'/profiles/'];
    
    % list files
    list = dir(float_dir);
    prof_list = {list.name}; % files names
    prof_list2{ifloat} = {list.name};
    
    % float not found
    if isempty(prof_list)
        disp(prof_list)
        fprintf(2,'Float not found\n') 
        cycle_number{ifloat}{:} = NaN;
        DM_done{ifloat}{:} = NaN;
        var_correction{ifloat}{:} = NaN;
        profile_var_qc{ifloat}{:} = '-';
%         var_qc_mode{ifloat}{:} = NaN;
%         profile_var_qc{ifloat}{:} = NaN;
        var_adjusted_bad{ifloat}{:} = NaN;
        n_cycle(ifloat) = 0;
        continue
    end  
    
    % not including descent profiles
    prof_list(contains(prof_list,'D.nc')) = [];
    % not including deep files
    prof_list(contains(prof_list,'SR')) = [];
    prof_list(contains(prof_list,'SD')) = [];
    prof_list(contains(prof_list,'MR')) = [];

    
    
% -b file pattern
    pattern_bio = "BR" + floats_list.WMO(ifloat,:);
    pattern_sfile = "S" + floats_list.WMO(ifloat,:); 

    if var_bio == 1
        prof_list = prof_list(contains(prof_list,pattern_bio));
%         prof_list = prof_list(3:end);
    else
        
    % not include -b files
    prof_list(contains(prof_list,'BR')) = [];
    prof_list(contains(prof_list,'BD')) = [];
    
    end
    
    pattern = "R" + floats_list.WMO(ifloat,:);
    pattern2 = "D" + floats_list.WMO(ifloat,:);
%     pattern2(ifloat) = "R" + floats_list.WMO(ifloat,:);
    
    if sum(contains(prof_list,pattern)) == 0
        disp("no -R nor -D file found")
        cycle_number{ifloat}{:} = NaN;
        DM_done{ifloat}{:} = NaN;
        var_correction{ifloat}{:} = NaN;
        profile_var_qc{ifloat}{:} = '-';
%         var_qc_mode{ifloat}{:} = NaN;
%         profile_var_qc{ifloat}{:} = NaN;
        var_adjusted_bad{ifloat}{:} = NaN;
        n_cycle(ifloat) = 0;
        continue
    end
    
    if sum(contains(prof_list,pattern)) > 0
        prof_list_filled = prof_list(contains(prof_list,pattern));
        prof_list_r{ifloat} = prof_list(contains(prof_list,pattern));
    end
    
    if sum(contains(prof_list,pattern2)) > 0
        prof_list_filled = prof_list(contains(prof_list,pattern2));
        prof_list_d{ifloat} = prof_list(contains(prof_list,pattern2));
        
        prof_list_filled = [prof_list_d{ifloat} prof_list_r{ifloat}];
        
    end
    prof_list_final{ifloat} = prof_list_filled;
    
    
    %         prof_list = prof_list(3:end);
    
    
%     prof_list(contains(prof_list,'BR')) = []; % comented because DOXY studied here
%     prof_list = prof_list(3:end);
    
    n_cycle(ifloat) = length(prof_list_filled);
      
        
    for icycle = 1:n_cycle(ifloat) % files loop (cycle loop)
        
        % get cycle number from file name
        cy = strsplit(char(prof_list_filled(icycle)),'_');
        cycle_number{ifloat}{icycle} = str2double(erase(char(cy(2)),'.nc'));
        
        % DM done or not
        DM_done{ifloat}{icycle} = double(contains(prof_list_filled(icycle),'D'));
        
        % file path
        A{ifloat}{icycle}= [float_dir,char(prof_list_filled(icycle))];
        file_path = [float_dir,char(prof_list_filled(icycle))];
        % get varibles var and var_adjusted and QC
        DATA_MODE{ifloat}{icycle} = ncread(file_path, 'DATA_MODE');
        
        file_info = ncinfo(file_path);
        len_var = size(file_info.Variables,2);
        name = cell(1,len_var);
        for f=1:len_var
            name{f} = file_info.Variables(f).Name;
        end
        
        if sum(contains(name,variable)) == 0
            
%             disp([file_path ' == variable not found'])
            var = {nan};
            var_adjusted = {nan};
            profile_qc = {nan};
            var_qc = {nan};
            profile_var_qc{ifloat}{icycle} = {'-'};
            var_correction{ifloat}{icycle} = {nan};
            var_qc_mode{ifloat}{icycle} = {nan};
            
            continue
        else


%             disp([file_path ' == variable found'])
            var = ncread(file_path,variable);
            var_adjusted = ncread(file_path,[variable,'_ADJUSTED']);
            profile_qc = ncread(file_path,['PROFILE_',variable,'_QC']);
            var_qc = ncread(file_path,[variable,'_ADJUSTED_QC']);

            A1{ifloat}{icycle} = ncread(file_path,variable);
            A2{ifloat}{icycle} = ncread(file_path,[variable,'_ADJUSTED']);
            A3{ifloat}{icycle} = ncread(file_path,['PROFILE_',variable,'_QC']);
            A4{ifloat}{icycle} = ncread(file_path,[variable,'_ADJUSTED_QC']);

            % if var_ADJUSTED_QC is fill value, we use var_QC
            if isspace(var_qc)
%                 var_qc = ncread(file_path,'var_QC');
                var_qc = ncread(file_path,[variable,'_QC']);
            end


            % delate second column
            [~,col_qc] = size(var_qc);
            if col_qc > 1
                var_qc(:,2) = [];
            end
            var_qc = str2num(var_qc);
            
  
            % calculate diference
            if isempty(var)
                var = 1;
            elseif isempty(var_adjusted)
                var_adjusted = 1;

            else
                var_corr = mean(var_adjusted - var,'omitnan');
                var_corr = mean(var_corr,'omitnan');
            
            % struc
                var_correction{ifloat}{icycle} = var_corr;
                profile_var_qc{ifloat}{icycle} = profile_qc(1);
            end

            % the most frequently occurring value in var_qc
            if isempty(var_qc)
                var_qc_mode{ifloat}{icycle} = NaN;
            else
                var_qc_mode{ifloat}{icycle} =  mode(var_qc);
            end
        end 
    end
    
    % var_correction(:,ifloat)'
    % var_adjusted_qc(:,ifloat)'
    
%     profile_var_qc = num2cell(profile_var_qc)
    
end


% Formatting results: convert cells to matrix struc
% fill with NaN or "-" when different number of cycles
max_cycle = max(n_cycle);
DATA.cycle_number = NaN(max_cycle,n_floats);
DATA.DM_done = NaN(max_cycle,n_floats);
DATA.var_correction = NaN(max_cycle,n_floats);
DATA.profile_var_qc = repmat('-',max_cycle,n_floats); % string array
% DATA.profile_var_qc = NaN(max_cycle,n_floats); % string array
DATA.var_qc_mode = NaN(max_cycle,n_floats); % percentage of 4 (bad data) in water column
%%
for ifloat = 1:n_floats
       disp(ifloat)
       DATA.cycle_number(1:n_cycle(ifloat),ifloat) = cat(1,cycle_number{ifloat}{:});
       DATA.DM_done(1:n_cycle(ifloat),ifloat) = cat(1,DM_done{ifloat}{:});
       
       % PROFILE_var_QC
       if any(cellfun(@isempty,profile_var_qc{ifloat}))
%        if any(cellfun(@isempty,profile_var_qc))    
           x = size(cat(1,profile_var_qc{ifloat}{:}),1);
           DATA.profile_var_qc(1:x,ifloat) = cat(1,profile_var_qc{ifloat}{:}); % string array
       else
           DATA.profile_var_qc(1:n_cycle(ifloat),ifloat) = cat(1,profile_var_qc{ifloat}{:}); % string array
       end
       
       % var_CORRECTION
       if any(cellfun(@isempty,var_correction{ifloat}))
           y = size(cat(1,var_correction{ifloat}{:}),1);
           DATA.var_correction(1:y,ifloat) = cat(1,var_correction{ifloat}{:});
       else
           DATA.var_correction(1:n_cycle(ifloat),ifloat) = cat(1,var_correction{ifloat}{:});
       end
       
       % var_QC_MODE 
       if isempty(var_qc_mode{ifloat})
           DATA.var_qc_mode(1:n_cycle(ifloat),ifloat) = cat(1,var_qc_mode{ifloat});
       elseif any(cellfun(@isempty,var_qc_mode{ifloat}))
           z = size(cat(1,var_qc_mode{ifloat}{:}),1);
           DATA.var_qc_mode(1:z,ifloat) = cat(1,var_qc_mode{ifloat}{:});
       else
           DATA.var_qc_mode(1:n_cycle(ifloat),ifloat) = cat(1,var_qc_mode{ifloat}{:});
       end    
end
% TODO : dimensions struc

DATA.profile_var_qc = num2cell(DATA.profile_var_qc);

%% Screen outputs
% get non zero adjustment floats reference
[~,c] = find(~isnan(DATA.var_correction) & DATA.var_correction ~= 0);
f_nonzero = DATA.WMO(c,:);
[~,irows] = unique(f_nonzero,'rows','stable');
f_nonzero = f_nonzero(irows,:);
disp(' ')
disp('Floats with non zero var correction:')
if isempty(f_nonzero)
    disp('(none)')
else
    disp(f_nonzero)
end

% floats with F cycles
% [f,c] = find(contains(cellstr(DATA.profile_var_qc), "F"));
[f,c] = find(contains(DATA.profile_var_qc, "F"));
f_Fcycles = DATA.WMO(c,:);
[~,irows] = unique(f_Fcycles,'rows');
f_Fcycles = f_Fcycles(irows,:);
disp(' ')
disp('Floats with bad cycles (F in profile_var_qc):')
if isempty(f_Fcycles)
    disp('(none)')
else
    disp(f_Fcycles)
end


%% Plot
close all

disp(' ')
disp('Making plots ...')

% check if output folder exits
working_date = erase(update_date,'-');
if ~exist(working_date, 'dir')
    mkdir(working_date)
end


% plots data

% Fig.1: var_correction diferent from 0 or not
% non zero correction
[~,ifloatN_1] = find(~isnan(DATA.var_correction) & DATA.var_correction ~= 0);
icycleN_1 = cat(1,DATA.cycle_number(~isnan(DATA.var_correction) & DATA.var_correction ~= 0));
sumN_1 = length(icycleN_1);
% zero correction
[~,ifloatZ_1] = find(DATA.var_correction == 0);
icycleZ_1 = cat(1,DATA.cycle_number(DATA.var_correction == 0));
sumZ_1 = length(icycleZ_1);

% Fig.3 : profile_var_qc
[~,ifloatA] = find(contains(DATA.profile_var_qc,'A')); % 100%
icycleA = cat(1,DATA.cycle_number(contains(DATA.profile_var_qc,'A')));
sumA = sum(sum(contains(DATA.profile_var_qc,'A')));
[~,ifloatB] = find(contains(DATA.profile_var_qc,'B')); % 100% - 75%
icycleB = cat(1,DATA.cycle_number(contains(DATA.profile_var_qc,'B')));
sumB = sum(sum(contains(DATA.profile_var_qc,'B')));
[~,ifloatC] = find(contains(DATA.profile_var_qc,'C')); % 75% - 50%
icycleC = cat(1,DATA.cycle_number(contains(DATA.profile_var_qc,'C')));
sumC = sum(sum(contains(DATA.profile_var_qc,'C')));
[~,ifloatD] = find(contains(DATA.profile_var_qc,'D')); % 50% - 25%
icycleD = cat(1,DATA.cycle_number(contains(DATA.profile_var_qc,'D')));
sumD = sum(sum(contains(DATA.profile_var_qc,'D')));
[~,ifloatE] = find(contains(DATA.profile_var_qc,'E')); % 25% - 0%
icycleE = cat(1,DATA.cycle_number(contains(DATA.profile_var_qc,'E')));
sumE = sum(sum(contains(DATA.profile_var_qc,'E')));
[~,ifloatF] = find(contains(DATA.profile_var_qc,'F')); % no good data
icycleF = cat(1,DATA.cycle_number(contains(DATA.profile_var_qc,'F')));
sumF = sum(sum(contains(DATA.profile_var_qc,'F')));

% Fig.4 DMQC done
[~,ifloatY_4] = find(DATA.DM_done == 1); % DM done
icycleY_4 = cat(1,DATA.cycle_number(DATA.DM_done == 1));
sumY_4 = length(icycleY_4);
[~,ifloatN_4] = find(DATA.DM_done == 0); % DM not done
icycleN_4 = cat(1,DATA.cycle_number(DATA.DM_done == 0));
sumN_4 = length(icycleN_4);

% Fig.5 : Most frequently occuring QC flac in var
[~,ifloat0] = find(DATA.var_qc_mode == 0); % 0: No QC was performed
icycle0 = cat(1,DATA.cycle_number(DATA.var_qc_mode == 0));
sum0 = sum(sum(DATA.var_qc_mode == 0));
[~,ifloat1] = find(DATA.var_qc_mode == 1); % 1: Good data
icycle1 = cat(1,DATA.cycle_number(DATA.var_qc_mode == 1));
sum1 = sum(sum(DATA.var_qc_mode == 1));
[~,ifloat2] = find(DATA.var_qc_mode == 2); % 2: Probably good data
icycle2 = cat(1,DATA.cycle_number(DATA.var_qc_mode == 2));
sum2 = sum(sum(DATA.var_qc_mode == 2));
[~,ifloat3] = find(DATA.var_qc_mode == 3); % 3: Bad data that are potentially correctable
icycle3 = cat(1,DATA.cycle_number(DATA.var_qc_mode == 3));
sum3 = sum(sum(DATA.var_qc_mode == 3));
[~,ifloat4] = find(DATA.var_qc_mode == 4); % 4: Bad data
icycle4 = cat(1,DATA.cycle_number(DATA.var_qc_mode == 4));
sum4 = sum(sum(DATA.var_qc_mode == 4));
[~,ifloat5] = find(DATA.var_qc_mode == 5); % 5: Value changed
icycle5 = cat(1,DATA.cycle_number(DATA.var_qc_mode == 5));
sum5 = sum(sum(DATA.var_qc_mode == 5));
[~,ifloat8] = find(DATA.var_qc_mode == 8); % 8: Estimated value
icycle8 = cat(1,DATA.cycle_number(DATA.var_qc_mode == 8));
sum8 = sum(sum(DATA.var_qc_mode == 8));
[~,ifloat9] = find(DATA.var_qc_mode == 9); % 9: Missing value
icycle9 = cat(1,DATA.cycle_number(DATA.var_qc_mode == 9));
sum9 = sum(sum(DATA.var_qc_mode == 9));



%% Fig.1 : var_correction diferent from 0 or not

if n_floats > floats_per_fig % if so many floats, figures are divided
    % split figure
    [hdl] = split_figure(floats_per_fig, DATA.WMO, {icycleY_4, ifloatY_4, icycleN_4, ifloatN_4, icycleZ_1, ifloatZ_1, icycleN_1, ifloatN_1});
else
    % make simple plot
    figure
    hdl{1,1} = scatter(icycleY_4,ifloatY_4,10,[0.8242 0.8242 0.8242],'filled');
    hold on
    hdl{1,3} = scatter(icycleN_4,ifloatN_4,10,[0.7945 1.0000 0.7945],'filled');
    hdl{1,5} = scatter(icycleZ_1,ifloatZ_1,10,'b','filled');
    hdl{1,7} = scatter(icycleN_1,ifloatN_1,10,'r','filled');
    % y axis format
    yticks(1:n_floats)
    yticklabels(cellstr(DATA.WMO(1:n_floats,:)))
    ylim([0,n_floats+0.5])
    
end

% figure format
n_fig = get(gcf,'Number');
for ifig = 1 : n_fig % figures loop
     
    figure(ifig)
    
    % full screen
    set(gcf, 'Position', get(0, 'Screensize'));
    % figure name
    set(gcf,'Name',['Non zero ' variable ' correction [Lot ' num2str(ifig) ']'])
    % axis labels size
    set(gca, 'FontSize', 12)
    %title
    if ifig == n_fig % last lot (diferent number of floats)
        title(['\fontsize{20}Non zero ' variable 'correction (updated ' update_date ')' newline ...
            '\rm\fontsize{18}[Lot ' num2str(ifig) ': ' num2str(n_floats - (ifig-1)*floats_per_fig) ' floats]'])
    else
        title(['\fontsize{20}Non zero ' variable 'correction (updated ' update_date ')' newline ...
            '\rm\fontsize{18}[Lot ' num2str(ifig) ': ' num2str(floats_per_fig) ' floats]'])
    end
     % colorbar
     lh = legend([hdl{ifig,5},hdl{ifig,7},hdl{ifig,3},hdl{ifig,1}],[variable ' correction = 0' newline '(',num2str(sumZ_1),' profiles)'],...
        [variable ' correction \neq 0' newline '(',num2str(sumN_1),' profiles)'],...
        'Real time cycle not corrected',...
        'DM cycle not corrected',...
        'Location','bestoutside');
        set(lh,'FontSize',10);
     % y ticks size
     yl = get(gca,'YLabel');
     ylFontSize = get(yl,'FontSize');
     yAY = get(gca,'YAxis');
     set(yAY,'FontSize', fontsize_wmo)
     set(yAY,'Color', 'k')
     ylabel('Float WMO','FontSize', 18)
     % x axis
     %xticks([1,10:10:max_cycle])
     xlabel('Cycle number','FontSize', 18)
     xlim([0,max_cycle+0.5])
     % background color
     set(gcf,'color','w');
     % box
     box on
     
     % save figure
     
     out_name = [pwd '/' working_date '/' project_name '_' variable 'CorrectionStatus_lot' num2str(ifig) '_' working_date '.png'];
     export_fig(out_name)
     
end  % figures loop



%% Fig.2 : var_correction for corrected floats

if ~isempty(f_nonzero) % if there is not nonzero correction cycles, this figure is not plotted
    %figure
    figure
    % full screen
    set(gcf, 'Position', get(0, 'Screensize'));
    % figure name
    set(gcf,'Name', ['Mean ' variable ' correction per cycle'])
    set(gcf,'Renderer', 'painters')
    hold on
        
    cnt = 0;
    max_cycles_fig = 0;
    for i=1:n_floats
    
        if contains(DATA.WMO(i,:),cellstr(f_nonzero)) % if correction is not 0
            
            cnt = cnt +1;
            max_cycles_fig = max([max_cycles_fig; DATA.cycle_number(:,i)]);
            %fprintf('        %s\n', DATA.WMO(i,:))
            % dm cycles
            index =ismember(ifloatY_4,i);
            h1 = scatter(icycleY_4(index),cnt*ones(1,sum(index)),8,'o','r');
            hold on
            
            % correction
            %DATA.var_correction(:,i)'
            h2 = scatter(DATA.cycle_number(:,i),cnt*ones(1,max_cycle), 50, DATA.var_correction(:,i),'filled');
            %pause()
            % real time cycles
            index =ismember(ifloatN_4,i); 
            h3 = scatter(icycleN_4(index),cnt*ones(1,sum(index)),8,'k','filled');
            
            % new line most frequent
            index =ismember(ifloat0,i);
            h4 = scatter(icycle0(index),(cnt-0.3)*ones(1,sum(index)),10,[0.8242 0.8242 0.8242],'filled'); % use grey
            index =ismember(ifloat1,i);
            h5 = scatter(icycle1(index),(cnt-0.3)*ones(1,sum(index)),10,'b','filled');
            index =ismember(ifloat2,i);
            h6 = scatter(icycle2(index),(cnt-0.3)*ones(1,sum(index)),10,'g','filled');
            index =ismember(ifloat3,i);
            h7 = scatter(icycle3(index),(cnt-0.3)*ones(1,sum(index)),10,[1.0000 0.5469 0],'filled');
            index =ismember(ifloat4,i);
            h8 = scatter(icycle4(index),(cnt-0.3)*ones(1,sum(index)),10,'r','filled');
            index =ismember(ifloat5,i);
            h9 = scatter(icycle5(index),(cnt-0.3)*ones(1,sum(index)),10,[0.5625 0.9297 0.5625],'filled');
            index =ismember(ifloat8,i);
            h10 = scatter(icycle8(index),(cnt-0.3)*ones(1,sum(index)),10,[0.4961 1.0000 0.8281],'filled');
            index =ismember(ifloat9,i);
            h11 = scatter(icycle9(index),(cnt-0.3)*ones(1,sum(index)),10,'k','filled');
            
        end
    
    end
    
   
    % figure format
    title(sprintf(['Mean ' variable ' correction per cycle (updated %s, only corrected floats)'],update_date),'FontSize', 20)
    % axis labels size
    set(gca, 'FontSize', 10)
    % colorbar
    h = colorbar;
    ylabel(h, [variable ' correction ' units],'FontSize', 12)
    % y axis
    yticks(1:length(f_nonzero))
    yticklabels(cellstr(f_nonzero))
    ylabel('Float WMO','FontSize', 12)
    ylim([0,cnt+1])
    % x axis
    %xticks([1,10:10:max_cycle])
    hlabel = xlabel('Cycle number','FontSize', 12);
    xlim([0,max_cycles_fig+0.5])
    % figure position (need space for legend)
    set(gca,'Units','pixels')
    set(gca,'Position',get(gca,'Position').*[1 2 1 0.8])
    % background color
    set(gcf,'color','w');
    % box
    box on
    
    % legend
    [hl(3).leg, hl(3).obj, hl(3).hout, hl(3).mout] = ...
    legendflex([h4,h5,h6,h7,h8,h9,h10,h11], {'0: No QC','1: Good data','2: Probably good data','3: Bad data potentially correctable',...
    '4: Bad data','5: Value changed','8: Estimated value','9: Missing value'},...
    'anchor', [5 3], ...
    'buffer', [0 -80], ...
    'fontsize', 10, ...
    'ncol', 2, ...%'box', 'off', ...
    'title', ['\bfMost frequent' variable ' QC flag \rm(2nd line)']);
    [hl(2).leg, hl(2).obj, hl(2).hout, hl(2).mout] = ...
    legendflex([h2,h3,h1], {[variable ' correction'],'Real time cycle','DM cycle not corrected'}, ...
    'ref', hl(3).leg, ...
    'anchor', [1 3], ...
    'buffer',  [-5 0], ...
    'fontsize',10, ...%'box', 'off', ...
    'title', '\bf Correction \rm(1st line)');

    % colormap in legend
    hold off
    % get the position of the legend, and calculate the place for the colormaps:
    pos = hl(2).leg.Position.*[1.02 1.55 0.18 0.21];
    % Create a 'picture' of what you want to appear in the legend:
    legax = axes('units','pixels','Position',pos); % place the new picture above the legend
    cmap = parula;
    imagesc(legax,1:max(size(parula))) % Create the picture
    colormap(cmap) % appy custom colormap
    axis off % remove all axes details


    % save figure
    out_name = [pwd '/' working_date '/' project_name '_' variable 'CorrectionValues_' working_date '.png'];
    export_fig(out_name,'-painters')

end



%% Fig.3 : profile_var_qc

start_fig = get(gcf,'Number') +1;

if n_floats > floats_per_fig % if so many floats, figures are divided
    % split figure
    split_figure(floats_per_fig, DATA.WMO, {icycleA,ifloatA,...
        icycleB,ifloatB,...
        icycleC,ifloatC,...
        icycleD,ifloatD,...
        icycleE,ifloatE,...
        icycleF,ifloatF});
else
    % make simple plot
    figure
    scatter(icycleA,ifloatA,10,'o','b','filled')
    hold on
    scatter(icycleB,ifloatB,10,'o','c','filled')
    scatter(icycleC,ifloatC,10,'o','g','filled')
    scatter(icycleD,ifloatD,10,'o','y','filled')
    scatter(icycleE,ifloatE,10,'o','m','filled')
    scatter(icycleF,ifloatF,10,'o','r','filled')
    % y axis format
    yticks(1:n_floats)
    yticklabels(cellstr(DATA.WMO(1:n_floats,:)))
    ylim([0,n_floats+0.5])
end

% figure format
n_fig = get(gcf,'Number');
for ifig = start_fig : n_fig % figures loop
     
    figure(ifig)
    lot = ifig - start_fig +1;
    % full screen
    set(gcf, 'Position', get(0, 'Screensize'));
    % figure name
    set(gcf,'Name',[variable ' quality control flag [Lot '  num2str(lot) ']'])
    % axis labels size
    set(gca, 'FontSize', 12)
    %title
    if ifig == n_fig % last lot (diferent number of floats)
        title(['\fontsize{20}' variable ' quality control flag (updated ' update_date ')' newline ...
            '\rm\fontsize{18}[Lot ' num2str(lot) ': ' num2str(n_floats - (lot-1)*floats_per_fig) ' floats]'])
    else
        title(['\fontsize{20}' variable ' quality control flag (updated ' update_date ')' newline ...
            '\rm\fontsize{18}[Lot ' num2str(lot) ': ' num2str(floats_per_fig) ' floats]'])
    end
     % colorbar
     lh = legend(['A: 100% good data' newline '(',num2str(sumA),' profiles)'],...
            ['B: 100% - 75% good data' newline '(',num2str(sumB),' profiles)'],...
            ['C: 75% - 50% good data' newline '(',num2str(sumC),' profiles)'],...
            ['D: 50% - 25% good data' newline '(',num2str(sumD),' profiles)'],...
            ['E: 25% - 0% good data' newline '(',num2str(sumE),' profiles)'],...
            ['F: no good data' newline '(',num2str(sumF),' profiles)'],...
            'Location','bestoutside');
     set(lh,'FontSize',10);
     % y ticks size
     yl = get(gca,'YLabel');
     ylFontSize = get(yl,'FontSize');
     yAY = get(gca,'YAxis');
     set(yAY,'FontSize', fontsize_wmo)
     set(yAY,'Color', 'k')
     ylabel('Float WMO','FontSize', 18)
     % x axis
     %xticks([1,10:10:max_cycle])
     xlabel('Cycle number','FontSize', 18)
     xlim([0,max_cycle+0.5])
     % background color
     set(gcf,'color','w');
     % box
     box on
     
     % save figure
     out_name = [pwd '/' working_date '/' project_name '_' variable 'QCscatter_lot' num2str(lot) '_' working_date '.png'];
     export_fig(out_name)
     
end  % figures loop



%% Fig.4 DMQC done

start_fig = get(gcf,'Number') +1;

if n_floats > floats_per_fig % if so many floats, figures are divided
    % split figure
    split_figure(floats_per_fig, DATA.WMO, {icycleN_4,ifloatN_4,icycleY_4,ifloatY_4});
else
    % make simple plot
    figure
    scatter(icycleN_4,ifloatN_4,10,'o','b','filled')
    hold on
    scatter(icycleY_4,ifloatY_4,10,'o','g','filled')
    % y axis format
    yticks(1:n_floats)
    yticklabels(cellstr(DATA.WMO(1:n_floats,:)))
    ylim([0,n_floats+0.5])
end

% figure format
n_fig = get(gcf,'Number');
for ifig = start_fig : n_fig % figures loop
     
    figure(ifig)
    lot = ifig - start_fig +1;
    % full screen
    set(gcf, 'Position', get(0, 'Screensize'));
    % figure name
    set(gcf,'Name',['DMQC status [Lot '  num2str(lot) ']'])
    % axis labels size
    set(gca, 'FontSize', 12)
    %title
    if ifig == n_fig % last lot (diferent number of floats)
        title(['\fontsize{20}DMQC status (updated ' update_date ')' newline ...
            '\rm\fontsize{18}[Lot ' num2str(lot) ': ' num2str(n_floats - (lot-1)*floats_per_fig) ' floats]'])
    else
        title(['\fontsize{20}DMQC status (updated ' update_date ')' newline ...
            '\rm\fontsize{18}[Lot ' num2str(lot) ': ' num2str(floats_per_fig) ' floats]'])
    end
     % colorbar
     lh = legend(['DMQC not done' newline '(',num2str(sumN_4),' profiles)'],...
            ['DMQC done' newline '(',num2str(sumY_4),' profiles)'],...
            'Location','bestoutside');
     set(lh,'FontSize',10);
     % y ticks size
     yl = get(gca,'YLabel');
     ylFontSize = get(yl,'FontSize');
     yAY = get(gca,'YAxis');
     set(yAY,'FontSize', fontsize_wmo)
     set(yAY,'Color', 'k')
     ylabel('Float WMO','FontSize', 18)
     % x axis
     %xticks([1,10:10:max_cycle])
     xlabel('Cycle number','FontSize', 18)
     xlim([0,max_cycle+0.5])
     % background color
     set(gcf,'color','w');
     % box
     box on
     
     % save figure
     out_name = [pwd '/' working_date '/' project_name '_DMQCstatusScatter_lot' num2str(lot) '_' working_date '.png'];
     export_fig(out_name)
     
end  % figures loop



%% Fig.5 Most frequently occuring QC flag in var

start_fig = get(gcf,'Number') +1;

if n_floats > floats_per_fig % if so many floats, figures are divided
    % split figure
    split_figure(floats_per_fig, DATA.WMO, {icycle0,ifloat0,...
        icycle1,ifloat1,...
        icycle2,ifloat2,...
        icycle3,ifloat3,...
        icycle4,ifloat4,...
        icycle5,ifloat5,...
        icycle8,ifloat8,...
        icycle9,ifloat9,...
        });
else
    % make simple plot
    figure
    scatter(icycle0,ifloat0,10,[0.8242 0.8242 0.8242],'filled')
    hold on
    scatter(icycle1,ifloat1,10,'b','filled')
    scatter(icycle2,ifloat2,10,'g','filled')
    scatter(icycle3,ifloat3,10,[1.0000 0.5469 0],'filled')
    scatter(icycle4,ifloat4,10,'r','filled')
    scatter(icycle5,ifloat5,10,[0.5625 0.9297 0.5625],'filled')
    scatter(icycle8,ifloat8,10,[0.4961 1.0000 0.8281],'filled')
    scatter(icycle9,ifloat9,10,'k','filled')
    % y axis format
    yticks(1:n_floats)
    yticklabels(cellstr(DATA.WMO(1:n_floats,:)))
    ylim([0,n_floats+0.5])
end

% figure format
n_fig = get(gcf,'Number');
for ifig = start_fig : n_fig % figures loop
     
    figure(ifig)
    lot = ifig - start_fig +1;
    % full screen
    set(gcf, 'Position', get(0, 'Screensize'));
    % figure name
    set(gcf,'Name',['Most frequent QC flag [Lot '  num2str(lot) ']'])
    % axis labels size
    set(gca, 'FontSize', 12)
    %title
    if ifig == n_fig % last lot (diferent number of floats)
        title(['\fontsize{20}Most frequently occurring QC flag in ' variable ' (updated ' update_date ')' newline ...
            '\rm\fontsize{18}[Lot ' num2str(lot) ': ' num2str(n_floats - (lot-1)*floats_per_fig) ' floats]'])
    else
        title(['\fontsize{20}Most frequently occurring QC flag in ' variable ' (updated ' update_date ')' newline ...
            '\rm\fontsize{18}[Lot ' num2str(lot) ': ' num2str(floats_per_fig) ' floats]'])
    end
     % colorbar
     lh = legend(['0: No QC' newline '(',num2str(sum0),' profiles)'],...
            ['1: Good data' newline '(',num2str(sum1),' profiles)'],...
            ['2: Probably good data' newline '(',num2str(sum2),' profiles)'],...
            ['3: Bad data potentially correctable' newline '(',num2str(sum3),' profiles)'],...
            ['4: Bad data' newline '(',num2str(sum4),' profiles)'],...
            ['5: Value changed' newline '(',num2str(sum5),' profiles)'],...
            ['8: Estimated value' newline '(',num2str(sum8),' profiles)'],...
            ['9: Missing value' newline '(',num2str(sum9),' profiles)'],...
            'Location','bestoutside');
     set(lh,'FontSize',10);
     % y ticks size
     yl = get(gca,'YLabel');
     ylFontSize = get(yl,'FontSize');
     yAY = get(gca,'YAxis');
     set(yAY,'FontSize', fontsize_wmo)
     set(yAY,'Color', 'k')
     ylabel('Float WMO','FontSize', 18)
     % x axis
     %xticks([1,10:10:max_cycle])
     xlabel('Cycle number','FontSize', 18)
     xlim([0,max_cycle+0.5])
     % background color
     set(gcf,'color','w');
     % box
     box on
     
     % save figure
     out_name = [pwd '/' working_date '/' project_name '_mostfreq' variable 'QC_lot' num2str(lot) '_' working_date '.png'];
     export_fig(out_name)
     
end  % figures loop

fprintf('Figures saved in %s folder\n',working_date)
