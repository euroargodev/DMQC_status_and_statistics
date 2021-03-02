% get_DMQC_adjustment
% shows psal adjustement for a list of floats
%
% Input (below)
% - Profiles files from dac (same file structure as dac)
% - floats list: csv file, with in header WMO, RT (dac in charge of real time
% processing) and DM (institution in charge of delayed mode processing), 
% separated by ';'. Exemple: MOCCA_rt_dm.csv
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

                  
clear variables

% add paths (packages and auxiliary functions)
aux_functions_path = [pwd '/aux_functions'];
% aux_functions_path = '/home1/datahome/co_arg/agarciaj/DMQC_status/aux_functions';
addpath(genpath(aux_functions_path))


% INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
project_name = 'MOCCA'; % for png files names
%MOCCA_file_dir = '/home1/datahome/co_arg/agarciaj/DMQC_status/MOCCA_rt_dm_trafic.csv'
MOCCA_file_dir = '/home1/datahome/co_arg/rcancoue/ROMAIN/ANDREA/DMQC_status/trueMOCCA_rt_dm.csv'
dac_dir = '/home/ref-argo/gdac/dac'
update_date = datestr(now(),'yyyy-mm-dd'); % dac update date (if it is not 
% today because we are working with a snapshot, change it)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FORMAT OPTIONS
% if so many floats, figures are divided. Number of floats in each figure 
% and yticks size can be change here:
floats_per_fig = 40; % number of floats per figure (80 is the very limit)
fontsize_wmo = 10; % yaxis ticks font size (8 is the limit)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TODO: VAR IMPUT: PRES PSAL AND temp
% TODO : dimensions struc

% read MOCCA floats list
[MOCCA_list] = read_csv(MOCCA_file_dir,';');
%MOCCA_list.WMO = char(sort(cellstr(MOCCA_list.WMO)));
DATA.WMO = MOCCA_list.WMO;
n_floats = size(MOCCA_list.WMO,1);

disp(' ')
disp('Getting data from profiles files ...')
% initialisations
n_cycle=NaN(n_floats,1);
for ifloat = 1:n_floats % floats loop
    
    fprintf('        %s\n', MOCCA_list.WMO(ifloat,:))
    % floats directory path string
    dac = lower(cellstr(MOCCA_list.RT(ifloat,:)));
    float_dir = [dac_dir,'/',char(dac),'/',MOCCA_list.WMO(ifloat,:),'/profiles/'];
    
    % list files
    list = dir(float_dir);
    prof_list = {list.name}; % files names
    % float not found
    if isempty(prof_list)
        fprintf(2,'Float not found\n') 
        cycle_number{ifloat}{:} = NaN;
        DM_done{ifloat}{:} = NaN;
        psal_correction{ifloat}{:} = NaN;
        profile_psal_qc{ifloat}{:} = "-";
        psal_adjusted_bad{ifloat}{:} = NaN;
        n_cycle(ifloat) = 0;
        continue
    end  
    % not including descent profiles
    prof_list(contains(prof_list,'D.nc')) = [];
    % not including deep files
    prof_list(contains(prof_list,'SR')) = [];
    prof_list(contains(prof_list,'MR')) = [];
    prof_list(contains(prof_list,'BR')) = [];
    prof_list = prof_list(3:end);
    n_cycle(ifloat) = length(prof_list);
      
        
    for icycle = 1:n_cycle(ifloat) % files loop (cycle loop)
        
        % get cycle number from file name
        cy = strsplit(char(prof_list(icycle)),'_');
        cycle_number{ifloat}{icycle} = str2double(erase(char(cy(2)),'.nc'));
        
        % DM done or not
        DM_done{ifloat}{icycle} = double(contains(prof_list(icycle),'D'));
        
        % file path
        file_path = [float_dir,char(prof_list(icycle))];
        % get varibles psal and psal_adjusted and QC
        psal = ncread(file_path,'PSAL');
        psal_adjusted = ncread(file_path,'PSAL_ADJUSTED');
        profile_qc = ncread(file_path,'PROFILE_PSAL_QC');
        psal_qc = ncread(file_path,'PSAL_ADJUSTED_QC');
        % if PSAL_ADJUSTED_QC is fill value, we use PSAL_QC
        if isspace(psal_qc)
            psal_qc = ncread(file_path,'PSAL_QC');
        end
        % delate second column
        [~,col_qc] = size(psal_qc);
        if col_qc > 1
            psal_qc(:,2) = [];
        end
        psal_qc = str2num(psal_qc);
        
        % calculate diference
        psal_corr = mean(psal_adjusted - psal,'omitnan');
        psal_corr = mean(psal_corr,'omitnan');
        
        % struc
        psal_correction{ifloat}{icycle} = psal_corr;
        profile_psal_qc{ifloat}{icycle} = profile_qc(1);
        
        % the most frequently occurring value in psal_qc
        if isempty(psal_qc)
            psal_qc_mode{ifloat}{icycle} = NaN;
        else
            psal_qc_mode{ifloat}{icycle} =  mode(psal_qc);
        end
        
    end
    
    % psal_correction(:,ifloat)'
    % psal_adjusted_qc(:,ifloat)'
    
end


% Formatting results: convert cells to matrix struc
% fill with NaN or "-" when diferent number of cycles
max_cycle = max(n_cycle);
DATA.cycle_number = NaN(max_cycle,n_floats);
DATA.DM_done = NaN(max_cycle,n_floats);
DATA.psal_correction = NaN(max_cycle,n_floats);
DATA.profile_psal_qc = repmat("-",max_cycle,n_floats); % string array
DATA.psal_qc_mode = NaN(max_cycle,n_floats); % percentage of 4 (bad data) in water column
for ifloat = 1:n_floats
       DATA.cycle_number(1:n_cycle(ifloat),ifloat) = cat(1,cycle_number{ifloat}{:});
       DATA.DM_done(1:n_cycle(ifloat),ifloat) = cat(1,DM_done{ifloat}{:});
       DATA.psal_correction(1:n_cycle(ifloat),ifloat) = cat(1,psal_correction{ifloat}{:});
       DATA.profile_psal_qc(1:n_cycle(ifloat),ifloat) = cat(1,profile_psal_qc{ifloat}{:}); % string array
       DATA.psal_qc_mode(1:n_cycle(ifloat),ifloat) = cat(1,psal_qc_mode{ifloat}{:});
end
% TODO : dimensions struc


%% Screen outputs
% get non cero adjustment floats reference
[~,c] = find(~isnan(DATA.psal_correction) & DATA.psal_correction ~= 0);
f_nonzero = DATA.WMO(c,:);
[~,irows] = unique(f_nonzero,'rows','stable');
f_nonzero = f_nonzero(irows,:);
disp(' ')
disp('Floats with non zero PSAL correction:')
if isempty(f_nonzero)
    disp('(none)')
else
    disp(f_nonzero)
end

% floats with F cycles
[f,c] = find(contains(DATA.profile_psal_qc, 'F'));
f_Fcycles = DATA.WMO(c,:);
[~,irows] = unique(f_Fcycles,'rows');
f_Fcycles = f_Fcycles(irows,:);
disp(' ')
disp('Floats with bad cycles (F in profile_psal_qc):')
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

% Fig.1: psal_correction diferent from 0 or not
% non zero correction
[~,ifloatN_1] = find(~isnan(DATA.psal_correction) & DATA.psal_correction ~= 0);
icycleN_1 = cat(1,DATA.cycle_number(~isnan(DATA.psal_correction) & DATA.psal_correction ~= 0));
sumN_1 = length(icycleN_1);
% zero correction
[~,ifloatZ_1] = find(DATA.psal_correction == 0);
icycleZ_1 = cat(1,DATA.cycle_number(DATA.psal_correction == 0));
sumZ_1 = length(icycleZ_1);

% Fig.3 : profile_psal_qc
[~,ifloatA] = find(contains(DATA.profile_psal_qc,'A')); % 100%
icycleA = cat(1,DATA.cycle_number(contains(DATA.profile_psal_qc,'A')));
sumA = sum(sum(contains(DATA.profile_psal_qc,'A')));
[~,ifloatB] = find(contains(DATA.profile_psal_qc,'B')); % 100% - 75%
icycleB = cat(1,DATA.cycle_number(contains(DATA.profile_psal_qc,'B')));
sumB = sum(sum(contains(DATA.profile_psal_qc,'B')));
[~,ifloatC] = find(contains(DATA.profile_psal_qc,'C')); % 75% - 50%
icycleC = cat(1,DATA.cycle_number(contains(DATA.profile_psal_qc,'C')));
sumC = sum(sum(contains(DATA.profile_psal_qc,'C')));
[~,ifloatD] = find(contains(DATA.profile_psal_qc,'D')); % 50% - 25%
icycleD = cat(1,DATA.cycle_number(contains(DATA.profile_psal_qc,'D')));
sumD = sum(sum(contains(DATA.profile_psal_qc,'D')));
[~,ifloatE] = find(contains(DATA.profile_psal_qc,'E')); % 25% - 0%
icycleE = cat(1,DATA.cycle_number(contains(DATA.profile_psal_qc,'E')));
sumE = sum(sum(contains(DATA.profile_psal_qc,'E')));
[~,ifloatF] = find(contains(DATA.profile_psal_qc,'F')); % no good data
icycleF = cat(1,DATA.cycle_number(contains(DATA.profile_psal_qc,'F')));
sumF = sum(sum(contains(DATA.profile_psal_qc,'F')));

% Fig.4 DMQC done
[~,ifloatY_4] = find(DATA.DM_done == 1); % DM done
icycleY_4 = cat(1,DATA.cycle_number(DATA.DM_done == 1));
sumY_4 = length(icycleY_4);
[~,ifloatN_4] = find(DATA.DM_done == 0); % DM not done
icycleN_4 = cat(1,DATA.cycle_number(DATA.DM_done == 0));
sumN_4 = length(icycleN_4);

% Fig.5 : Most frequently occuring QC flac in PSAL
[~,ifloat0] = find(DATA.psal_qc_mode == 0); % 0: No QC was performed
icycle0 = cat(1,DATA.cycle_number(DATA.psal_qc_mode == 0));
sum0 = sum(sum(DATA.psal_qc_mode == 0));
[~,ifloat1] = find(DATA.psal_qc_mode == 1); % 1: Good data
icycle1 = cat(1,DATA.cycle_number(DATA.psal_qc_mode == 1));
sum1 = sum(sum(DATA.psal_qc_mode == 1));
[~,ifloat2] = find(DATA.psal_qc_mode == 2); % 2: Probably good data
icycle2 = cat(1,DATA.cycle_number(DATA.psal_qc_mode == 2));
sum2 = sum(sum(DATA.psal_qc_mode == 2));
[~,ifloat3] = find(DATA.psal_qc_mode == 3); % 3: Bad data that are potentially correctable
icycle3 = cat(1,DATA.cycle_number(DATA.psal_qc_mode == 3));
sum3 = sum(sum(DATA.psal_qc_mode == 3));
[~,ifloat4] = find(DATA.psal_qc_mode == 4); % 4: Bad data
icycle4 = cat(1,DATA.cycle_number(DATA.psal_qc_mode == 4));
sum4 = sum(sum(DATA.psal_qc_mode == 4));
[~,ifloat5] = find(DATA.psal_qc_mode == 5); % 5: Value changed
icycle5 = cat(1,DATA.cycle_number(DATA.psal_qc_mode == 5));
sum5 = sum(sum(DATA.psal_qc_mode == 5));
[~,ifloat8] = find(DATA.psal_qc_mode == 8); % 8: Estimated value
icycle8 = cat(1,DATA.cycle_number(DATA.psal_qc_mode == 8));
sum8 = sum(sum(DATA.psal_qc_mode == 8));
[~,ifloat9] = find(DATA.psal_qc_mode == 9); % 9: Missing value
icycle9 = cat(1,DATA.cycle_number(DATA.psal_qc_mode == 9));
sum9 = sum(sum(DATA.psal_qc_mode == 9));



%% Fig.1 : psal_correction diferent from 0 or not

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
    set(gcf,'Name',['Non zero PSAL correction [Lot ' num2str(ifig) ']'])
    % axis labels size
    set(gca, 'FontSize', 12)
    %title
    if ifig == n_fig % last lot (diferent number of floats)
        title(['\fontsize{20}Non zero PSAL correction (updated ' update_date ')' newline ...
            '\rm\fontsize{18}[Lot ' num2str(ifig) ': ' num2str(n_floats - (ifig-1)*floats_per_fig) ' floats]'])
    else
        title(['\fontsize{20}Non zero PSAL correction (updated ' update_date ')' newline ...
            '\rm\fontsize{18}[Lot ' num2str(ifig) ': ' num2str(floats_per_fig) ' floats]'])
    end
     % colorbar
     lh = legend([hdl{ifig,5},hdl{ifig,7},hdl{ifig,3},hdl{ifig,1}],['PSAL correction = 0' newline '(',num2str(sumZ_1),' profiles)'],...
        ['PSAL correction \neq 0' newline '(',num2str(sumN_1),' profiles)'],...
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
     out_name = [pwd '/' working_date '/' project_name '_psalCorrectionStatus_lot' num2str(ifig) '_' working_date '.png'];
     export_fig(out_name)
     
end  % figures loop



%% Fig.2 : psal_correction for corrected floats

if ~isempty(f_nonzero) % if there is not nonzero correction cycles, this figure is not plotted
    %figure
    figure
    % full screen
    set(gcf, 'Position', get(0, 'Screensize'));
    % figure name
    set(gcf,'Name','Mean PSAL correction per cycle')
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
            %DATA.psal_correction(:,i)'
            h2 = scatter(DATA.cycle_number(:,i),cnt*ones(1,max_cycle), 15, DATA.psal_correction(:,i),'filled');
            %pause()
            % real time cycles
            index =ismember(ifloatN_4,i); 
            h3 = scatter(icycleN_4(index),cnt*ones(1,sum(index)),11,'o','k');
            
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
    title(sprintf('Mean PSAL correction per cycle (updated %s, only corrected floats)',update_date),'FontSize', 20)
    % axis labels size
    set(gca, 'FontSize', 10)
    % colorbar
    h = colorbar;
    ylabel(h, 'PSAL correction (PSU)','FontSize', 12)
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
    'title', '\bfMost frequent PSAL QC flag \rm(2nd line)');
    [hl(2).leg, hl(2).obj, hl(2).hout, hl(2).mout] = ...
    legendflex([h2,h3,h1], {'PSAL correction','Real time cycle','DM cycle not corrected'}, ...
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
    out_name = [pwd '/' working_date '/' project_name '_psalCorrectionValues_' working_date '.png'];
    export_fig(out_name,'-painters')

end



%% Fig.3 : profile_psal_qc

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
    set(gcf,'Name',['PSAL quality control flag [Lot '  num2str(lot) ']'])
    % axis labels size
    set(gca, 'FontSize', 12)
    %title
    if ifig == n_fig % last lot (diferent number of floats)
        title(['\fontsize{20}PSAL quality control flag (updated ' update_date ')' newline ...
            '\rm\fontsize{18}[Lot ' num2str(lot) ': ' num2str(n_floats - (lot-1)*floats_per_fig) ' floats]'])
    else
        title(['\fontsize{20}PSAL quality control flag (updated ' update_date ')' newline ...
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
     out_name = [pwd '/' working_date '/' project_name '_psalQCscatter_lot' num2str(lot) '_' working_date '.png'];
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



%% Fig.5 Most frequently occuring QC flag in PSAL

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
        title(['\fontsize{20}Most frequently occurring QC flag in PSAL (updated ' update_date ')' newline ...
            '\rm\fontsize{18}[Lot ' num2str(lot) ': ' num2str(n_floats - (lot-1)*floats_per_fig) ' floats]'])
    else
        title(['\fontsize{20}Most frequently occurring QC flag in PSAL (updated ' update_date ')' newline ...
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
     out_name = [pwd '/' working_date '/' project_name '_mostfreqpsalQC_lot' num2str(lot) '_' working_date '.png'];
     export_fig(out_name)
     
end  % figures loop

fprintf('Figures saved in %s folder\n',working_date)
