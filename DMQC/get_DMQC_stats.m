% get_DMQC_stats
% gets DMQC statistics for given floats
%
% Input (below)
% - argo_profile_detailled_index.txt
% - ar_greylist.txt
% - floats list: csv file, with in header WMO, RT (dac in charge of real time
% processing) and DM (institution in charge of delayed mode processing), 
% separated by ';'. Exemple: MOCCA_rt_dm.csv
% Paths to these files should be modified below
% - Dprof: if 1 descent profiles are consider, if 0 descent profiles are not
% consider.
% - sage: threshold (in days) for floats and observations statistics.
%
% Output
% -Figures saved in (update date) folder (ex: 20181210) 
% -DMQCstatus_{update_date}.txt containing all calculations
%
% Auxiliary functions needed:
%    read_csv
%    get_data_from_index_detailled
%    export_fig package
%    plotBarStackGroups
%
% NOTE 
% (1) Reading index file may take 2 or 3 minutes
% (2) It is better to use only ascent profiles (Dprof = 0) for coherence 
% with figures from get_DMQC_adjustment.m script 
%
% Modified on 2019/11/26



clear variables

% add paths (packages and auxiliary functions)
aux_functions_path = [pwd '/aux_functions'];
% aux_functions_path = '/home1/datahome/co_arg/agarciaj/DMQC_status/aux_functions';
addpath(genpath(aux_functions_path))


% INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
index_file_dir = '/home/ref-argo/gdac/etc/argo_profile_detailled_index.txt'
greylist_dir = '/home/ref-argo/gdac/ar_greylist.txt'
list_file_dir = '/home1/datahome/co_arg/agarciaj/DMQC_status/floats_noe_DAC.txt'
%list_file_dir = '/home1/datahome/co_arg/rcancoue/ROMAIN/ANDREA/DMQC_status/trueMOCCA_rt_dm.csv'
% file separated with ;
project_name = ' ';
Dprof= 0; % 0 not including D prof, 1 including D prof
sage = 365; % threshold (in days). More than 1 year (floats and obs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TODO: percentage of DMQC for each float (txt output) and warning when 0%
% and >1 year
% TODO: struc dimensions
% TODO: (warning, sametimes there are intercalate profiles with no DM done)
% in txt
% TODO: initialisations in loops
% TODO: delimiter in list inputs


%% Read input files
% index files path
% read csv with MOCCA floats
[Floats_list] = read_csv(list_file_dir,',');
disp(' ')
% read csv with grey list
[grey_list] = read_csv(greylist_dir,',');
disp(' ')
% read index file and get data for MOCCA floats
[IndexData] = get_data_from_index(index_file_dir, cellstr(Floats_list.WMO)',Dprof);

n_floats = size(Floats_list.WMO,1);
max_cycles = size(IndexData.latitude.data,1);


%% Output data calculation 
% Creating struct Output with one field per calculated variable.
% Variables dimension is 1*n_floats


disp(' ')
disp('Calculating statistics ...') 

Output.WMO = Floats_list.WMO; % floats WMO
Output.first_date = squeeze(IndexData.date.data(1,:,:))'; % date of first profile (it is not launch date)
Output.RT = Floats_list.RT; % real time responsible
Output.DM = Floats_list.DM; % delayed mode responsible
Output.index_update = IndexData.index_update.data; % index file update date



%%%%%%%%%% DMQC status %%%%%%%%%%

% find D files
match = contains(string(IndexData.file.data),'/D'); % n_cycles x n_floats

% initialisations
obs_age = NaN(max_cycles,n_floats);
Output.DM_done = NaN(n_floats,1);
Output.number_DMprof = NaN(n_floats,1);
Output.obs_more1year = NaN(n_floats,1);
Output.obs_more1year_noDM = NaN(n_floats,1);
Output.obs_more1year_DMQC = NaN(n_floats,1);

for i = 1:n_floats
    
    sum_DM = sum(match(:,i));
    
    % -- DMQC done or not --
    Output.DM_done(i) = (sum_DM ~= 0); % DM done = 1, DM not done = 0 (at least 1 D file)
    Output.number_DMprof(i) = sum_DM; % number of D profiles
    
    % -- Observations > 1 year --
    out_range = contains(string(IndexData.date.data(:,:,i)),'-');
    date_num = datenum(IndexData.date.data(:,:,i),'yyyymmddHHMMSS');
    date_num(out_range) = NaN;
    obs_age(:,i) = datenum(Output.index_update,'yyyymmddHHMMSS') - date_num;
    Output.obs_more1year(i) = sum(obs_age(:,i) > sage); % obs age > 1year = 1, else = 0
    
    % -- Observations > 1 year and DMQC --
    Output.obs_more1year_noDM(i) = sum((obs_age(:,i) > sage).*(~match(:,i))); % > 1 year * not DM done
    Output.obs_more1year_DMQC(i) = sum((obs_age(:,i) > sage).*(match(:,i))); % > 1 year * DMQC done
    
    % -- Year of observation --
    % get all different years
    obs_year = IndexData.date.data(:,1:4,i);
    obs_year(out_range,:) = NaN;
    obs_year = str2num(obs_year);
    year_names = unique(obs_year);
    % str to compare
    obs_year_str = cellstr(IndexData.date.data(:,1:4,i));
    % count obs for each year
    for iyear = 1 : length(year_names)
        Output.(['obs_' num2str(year_names(iyear))])(i) = sum(contains(obs_year_str, num2str(year_names(iyear))));
        Output.(['obs_' num2str(year_names(iyear)) '_DMQC'])(i) = sum(contains(obs_year_str, num2str(year_names(iyear))).*(match(:,i)));
    end
    
    % -- Last DM profile --
    % (warning, sametimes there are intercalate profiles with no DM done)
    last_date = string(IndexData.date.data(:,:,i));
    last_date = last_date(match(:,i));
    if isempty(last_date)
        Output.date_last_DMprof(i) = NaN;
    else
        Output.date_last_DMprof(i) = last_date(end);
    end
    
end

% -- Formating --
% problems with missing string when writting txt
imis=ismissing(Output.date_last_DMprof);
Output.date_last_DMprof(imis) = '--------------';
Output.date_last_DMprof = char(Output.date_last_DMprof);
Output.date_last_DMprof = squeeze(Output.date_last_DMprof)';
% formating obs per year
fields_numbers = regexp(fieldnames(Output),'\d*','Match');
year_names = unique(str2double(cellfun(@(x) char(x), fields_numbers,'UniformOutput',false)));
year_names(isnan(year_names)) = [];
year_names(year_names<2) = [];
for iyear = 1: length(year_names)
    if length(Output.(['obs_' num2str(year_names(iyear))])) < n_floats
        Output.(['obs_' num2str(year_names(iyear))])(end+1:n_floats) = zeros(n_floats - length(Output.(['obs_' num2str(year_names(iyear))])),1);
    end
     if length(Output.(['obs_' num2str(year_names(iyear)) '_DMQC'])) < n_floats
        Output.(['obs_' num2str(year_names(iyear)) '_DMQC'])(end+1:n_floats) = zeros(n_floats - length(Output.(['obs_' num2str(year_names(iyear)) '_DMQC'])),1);
    end
end


% -- Floats age --
out_range = contains(string(Output.first_date),'-');
Output.float_age = datenum(Output.index_update,'yyyymmddHHMMSS') - datenum(Output.first_date,'yyyymmddHHMMSS'); % float age in days (from index update date)
Output.float_age(out_range) = NaN;
% -- Floats > 1 year --
Output.float_more1year = (Output.float_age > sage); % float age > 1year = 1, else = 0
% -- Number of observations --
Output.obs_number = IndexData.n_obs.data'; % number of observations per float
% -- Percentage of DMQC --
Output.percen_DMQC = Output.number_DMprof./Output.obs_number*100; 


%%%%%%%%%% Data Quality %%%%%%%%%%

% -- psal qc --
for ifloat = 1:n_floats
    % number of profiles with qc = 'A' for each float
    Output.psal_qcA_prof(ifloat) = sum(contains(string(IndexData.profile_psal_qc.data(:,:,ifloat)),'A'));
    Output.psal_qcB_prof(ifloat) = sum(contains(string(IndexData.profile_psal_qc.data(:,:,ifloat)),'B')); % qc = 'B'
    Output.psal_qcC_prof(ifloat) = sum(contains(string(IndexData.profile_psal_qc.data(:,:,ifloat)),'C')); % qc = 'C'
    Output.psal_qcD_prof(ifloat) = sum(contains(string(IndexData.profile_psal_qc.data(:,:,ifloat)),'D')); % qc = 'D'
    Output.psal_qcE_prof(ifloat) = sum(contains(string(IndexData.profile_psal_qc.data(:,:,ifloat)),'E')); % qc = 'E'
    Output.psal_qcF_prof(ifloat) = sum(contains(string(IndexData.profile_psal_qc.data(:,:,ifloat)),'F')); % qc = 'F'
end


%%%%%%%%%% Grey list %%%%%%%%%%
for ifloat = 1:n_floats
    in_float = contains(string(grey_list.PLATFORMCODE),cellstr(Floats_list.WMO(ifloat,:)));
    Output.greylist(ifloat) = ~all(in_float == 0);
end


%% Calculations for graphs

% -- Data per DM operator --
% get all DM operators
[~,iuni] = unique(lower(cellstr(Output.DM))); % case insensitive
DM_operators = cellstr(Output.DM(iuni,:));
n_DM = length(DM_operators);

% initialisations
number_floats = NaN(1,n_DM);
n_floats_more1year = NaN(1,n_DM);
n_floats_DMQC = NaN(1,n_DM);
obs_number = NaN(1,n_DM);
n_obs_DMQC = NaN(1,n_DM);
n_obs_more1year = NaN(1,n_DM);
n_obs_more1year_DMQC = NaN(1,n_DM);

for idm = 1 : n_DM
    Output.(DM_operators{idm}) = contains(string(Output.DM),DM_operators{idm},'IgnoreCase',true);
    
    %%%%%%%%%%%%%%  Fig 1 : bar per float %%%%%%%%%%%%%% 
    number_floats(idm) = sum(Output.(DM_operators{idm}));
    n_floats_more1year(idm) = sum(Output.(DM_operators{idm}).*Output.float_more1year);
    n_floats_DMQC(idm) = sum(Output.(DM_operators{idm}).*Output.DM_done);
    
    %%%%%%%%%%%%%% Fig 2 : bar per observation %%%%%%%%%%%%%%
    obs_number(idm) = sum(Output.(DM_operators{idm}).*Output.obs_number);
    n_obs_DMQC(idm) = sum(Output.(DM_operators{idm}).*Output.number_DMprof);
    n_obs_more1year(idm) = sum(Output.(DM_operators{idm}).*Output.obs_more1year);
    n_obs_more1year_DMQC(idm) = sum(Output.(DM_operators{idm}).*Output.obs_more1year_DMQC);
end


%%%%%%%%%%%%%% Fig 3 : Data Quality %%%%%%%%%%%%%%
% per observation (all obs, not only DMQC done)
obs_psal_qc(1) = sum(Output.psal_qcA_prof); % A
obs_psal_qc(2) = sum(Output.psal_qcB_prof); % B
obs_psal_qc(3) = sum(Output.psal_qcC_prof); % C
obs_psal_qc(4) = sum(Output.psal_qcD_prof); % D
obs_psal_qc(5) = sum(Output.psal_qcE_prof); % E
obs_psal_qc(6) = sum(Output.psal_qcF_prof); % F

%%%%%%%%%%%%%% Fig 4 : Grey list  %%%%%%%%%%%%%%
total_DMfloats = sum(Output.DM_done);
total_greylist = sum(Output.greylist);
total_more1year = sum(Output.float_more1year);

%%%%%%%%%%%%%% Fig 5 : observation DMQC per year %%%%%%%%%%%%%%
n_obs_peryear = NaN(1,length(year_names));
n_obs_peryear_DMQC = NaN(1,length(year_names));
for iyear = 1: length(year_names)
    n_obs_peryear(iyear) = sum(Output.(['obs_' num2str(year_names(iyear))]));
    n_obs_peryear_DMQC(iyear) = sum(Output.(['obs_' num2str(year_names(iyear)) '_DMQC']));
end

%%%%%%%%%%%%%% Fig 6 : observation DMQC age histogram %%%%%%%%%%%%%%
obs_age_notDM = reshape(obs_age.*(~match),[],1);
obs_age_notDM(isnan(obs_age_notDM)) = [];
obs_age_notDM(obs_age_notDM == 0) = [];


%% Make graphics
close all

disp(' ')
disp('Making graphics...')

if Dprof== 0
    prof_included = '  [only ascent profiles]';
else
    prof_included = '  [including D profiles]';
end

% check if output folder exits
working_date = Output.index_update(1:8);
if ~exist(working_date, 'dir')
    mkdir(working_date)
end

% colors
bars_colors = [0.2422    0.1504    0.6603; ... % total obs/floats (dark blue)
               0.9769    0.9839    0.0805; ... % >1year obs/floats (yellow)
               0.1953    0.8008    0.1953; ... % DMQC done (green)
               0.5938    0.9833    0.5938; ... % >1year obs/floats + DMQC done (pale green)
               0.5430    0         0.5430; ... % psal qc A
               0.9375    0.5000    0.5000; ... % psal qc B
               1.0000    0.5469         0; ... % psal qc C
               0.8516    0.6445    0.1250; ... % psal qc D
               0.1250    0.6953    0.6641; ... % psal qc E
               0.5273    0.8047    0.9180; ... % psal qc F
               0.4102    0.4102    0.4102];    % grey list

%% %%%%%%%%%%%%%%  Fig 1 : bar per float %%%%%%%%%%%%%% 

figure(1)
% bigger figure
set(gcf, 'Position', [200, 200, 700, 600])
% figure name
set(gcf,'Name','DMQC per float')

if n_DM == 1
    hndl = bar(0:n_DM+1,[0 0 0; number_floats' n_floats_more1year' n_floats_DMQC'; 0 0 0], 1);
else
    hndl = bar(1:n_DM, [number_floats' n_floats_more1year' n_floats_DMQC'], 1);
end

% FIGURE FORMAT
% colors
set(hndl(1),'facecolor',bars_colors(1,:))
set(hndl(2),'facecolor',bars_colors(2,:))
set(hndl(3),'facecolor',bars_colors(3,:))
% xlabels
set(gca,'xtick',[1:n_DM],'xticklabel',DM_operators)
% title with update date
title(['Floats DMQC status (updated ',datestr(datenum(Output.index_update,'yyyymmddHHMMSS'),'yyyy-mm-dd'),')'])
ylabel('Number of floats')
% legend with total number
legend(['Number of floats (Total :',num2str(sum(number_floats)),')'],...
    ['Floats > 1 year (Total :',num2str(sum(n_floats_more1year)),')'],...
    ['Floats DMQC done (Total :',num2str(sum(n_floats_DMQC)),')', newline, '(at least 1 file)'])
% background color
set(gcf,'color','w');
% grid in y axis
ax = gca;
ax.YGrid = 'on';

% annotation: percentage of floats > 1 year with no DMQC done
floats_tobedone = sum(Output.float_more1year.*(~Output.DM_done))/sum(n_floats_more1year)*100;
H = figure(1);
set(H,'units','pix')
annotation('textbox', [0.8, 0.01, .1, .1], 'string', ['Floats (> 1 year) to be done: ',num2str(round(floats_tobedone)),'%'],...
    'FitBoxToText','on','verticalalignment', 'bottom','HorizontalAlignment', 'right','FontWeight','bold')

% save figure
out_name = [pwd '/' working_date '/' project_name '_float_DMQCstatus_' working_date '.png'];
export_fig(out_name)


%% %%%%%%%%%%%%%% Fig 2 : bar per observation %%%%%%%%%%%%%%

figure(2)
% bigger figure
set(gcf, 'Position', [200, 200, 700, 600])
% figure name
set(gcf,'Name','DMQC per observation')

if n_DM == 1
    stackData = cat(3,[0 0; n_obs_DMQC' n_obs_more1year_DMQC'; 0 0],[0 0; (obs_number-n_obs_DMQC)' (n_obs_more1year-n_obs_more1year_DMQC)'; 0 0]);
    hndl = plotBarStackGroups(stackData, [ {' '} DM_operators {' '}]);
else
    stackData = cat(3,[n_obs_DMQC' n_obs_more1year_DMQC'],[(obs_number-n_obs_DMQC)' (n_obs_more1year-n_obs_more1year_DMQC)']);
    hndl = plotBarStackGroups(stackData, DM_operators); 
end


% FIGURE FORMAT
hold on
hp = plot(1,0,'w'); % for comment in legend
% bar colors 
set(hndl(1,1),'facecolor',bars_colors(3,:))
set(hndl(1,2),'facecolor',bars_colors(1,:))
set(hndl(2,1),'facecolor',bars_colors(4,:))
set(hndl(2,2),'facecolor',bars_colors(2,:))
% xlabels
if n_DM == 1
    set(gca,'xtick',2,'xticklabel', DM_operators)
else
    set(gca,'xtick',[1:n_DM],'xticklabel', DM_operators)
end
% title with update date
title(['Observations DMQC status (updated ',datestr(datenum(Output.index_update,'yyyymmddHHMMSS'),'yyyy-mm-dd'),')'])
ylabel('Number of observations')
% legend with total number
legend([hndl(1,2), hndl(1,1), hndl(2,2), hndl(2,1), hp], ...
    {['Number of observations (Total: ',num2str(sum(obs_number)),')'],...
    ['Observations DMQC done (Total: ',num2str(sum(n_obs_DMQC)),')'], ...
    ['Observations > 1 year (Total: ',num2str(sum(n_obs_more1year)),')'],...
    ['Observations > 1 year' newline '  + DMQC done (Total: ',num2str(sum(n_obs_more1year_DMQC)),')'], ...
    prof_included})
% background color
set(gcf,'color','w');
% grid in y axis
ax = gca;
ax.YGrid = 'on';

% annotation: percentage of observations > 1 year with no DMQC done
obs_tobedone = sum(Output.obs_more1year_noDM)/sum(n_obs_more1year)*100;
H = figure(2);
set(H,'units','pix')
annotation('textbox', [0.8, 0.01, .1, .1], 'string', ['Observations (> 1 year) to be done: ',num2str(round(obs_tobedone)),'%'],...
    'FitBoxToText','on','verticalalignment', 'bottom','HorizontalAlignment', 'right','FontWeight','bold')

% save figure
out_name = [pwd '/' working_date '/' project_name '_obs_DMQCstatus_' working_date '.png'];
export_fig(out_name)


%% %%%%%%%%%%%%%% Fig 3 : Data Quality %%%%%%%%%%%%%%
    
figure(3)
% bigger figure
set(gcf, 'Position', [200, 200, 700, 600])
% figure name
set(gcf,'Name','Data Quality')

colormap(parula)
hold on
for i =1:6
   bar(i, obs_psal_qc(i), 'FaceColor',bars_colors(4+i,:))
end

% FIGURE FORMAT
hold on
plot(1,0,'w'); % for comment in legend
% xlabels
set(gca,'XTick',[])
% title with update date
title(['PSAL data quality control (updated ',datestr(datenum(Output.index_update,'yyyymmddHHMMSS'),'yyyy-mm-dd'),')'])
ylabel('Number of observations')
% legend with total number
legend(['A: 100% good data' newline '(' num2str(obs_psal_qc(1)) ' profiles)'],...
    ['B: 100% - 75% good data' newline '(' num2str(obs_psal_qc(2)) ' profiles)'],...
    ['C: 75% - 50% good data' newline '(' num2str(obs_psal_qc(3)) ' profiles)'],...
    ['D: 50% - 25% good data' newline '(' num2str(obs_psal_qc(4)) ' profiles)'],...
    ['E: 25% - 0% good data' newline '(' num2str(obs_psal_qc(5)) ' profiles)'],...
    ['F: no good data' newline '(' num2str(obs_psal_qc(6)) ' profiles)'],...
     prof_included)
% background color
set(gcf,'color','w');
% grid in y axis
ax = gca;
ax.YGrid = 'on';
box on

% save figure
out_name = [pwd '/' working_date '/' project_name '_psal_QCstats_' working_date '.png'];
export_fig(out_name)


%% %%%%%%%%%%%%%% Fig 4 : Grey list  %%%%%%%%%%%%%%

figure(4)
% bigger figure
set(gcf, 'Position', [200, 200, 700, 600])
% figure name
set(gcf,'Name','Grey list')

bar(1, n_floats, 'FaceColor', bars_colors(1,:))
hold on
bar(2, total_more1year, 'FaceColor', bars_colors(2,:))
bar(3, total_DMfloats, 'FaceColor', bars_colors(3,:))
bar(4, total_greylist, 'FaceColor', bars_colors(11,:))

% FIGURE FORMAT
% xlabels
set(gca,'XTick',[])
% title with update date
title(['Floats status (updated ',datestr(datenum(Output.index_update,'yyyymmddHHMMSS'),'yyyy-mm-dd'),')'])
ylabel('Number of floats')
% legend with total number
legend(['All floats (' num2str(n_floats) ')'],...
    ['Floats > 1 year (' num2str(total_more1year) ')'],...
    ['Floats with DMQC (' num2str(total_DMfloats) ')'],...
    ['Floats in grey list (' num2str(total_greylist) ')'],...
    'Location','southoutside')

% background color
set(gcf,'color','w');
% grid in y axis
ax = gca;
ax.YGrid = 'on';

% add WMO number for the greylist
Greylist_label = Output.WMO(find(Output.greylist'),:);
text(4, total_greylist, Greylist_label, 'HorizontalAlignment','center', 'VerticalAlignment','bottom')
% text(4-0.25, total_greylist-30, Greylist_label(1:22,:), 'HorizontalAlignment','center', 'VerticalAlignment','bottom')
% text(4+0.25, total_greylist-30, Greylist_label(23:end,:), 'HorizontalAlignment','center', 'VerticalAlignment','bottom')

% save figure
out_name = [pwd '/' working_date '/' project_name '_float_totals_' working_date '.png'];
export_fig(out_name)


%% %%%%%%%%%%%%%% Fig 5 : observation DMQC per year %%%%%%%%%%%%%%

figure(5)
% bigger figure
set(gcf, 'Position', [200, 200, 700, 600])
% figure name
set(gcf,'Name','Observations DMQC per year')

hndl = bar(1:length(year_names), [n_obs_peryear' n_obs_peryear_DMQC'], 1);

% FIGURE FORMAT
hold on
plot(1,0,'w'); % for comment in legend
% bars colors
set(hndl(1),'facecolor',bars_colors(1,:))
set(hndl(2),'facecolor',bars_colors(3,:))
% xlabels
set(gca,'xtick',[1:length(year_names)],'xticklabel',year_names)
% title with update date
title(['Observations DMQC status per year (updated ',datestr(datenum(Output.index_update,'yyyymmddHHMMSS'),'yyyy-mm-dd'),')'])
ylabel('Number of observations')
% legend with total number
legend(['Number of observations (Total: ',num2str(sum(n_obs_peryear)),')'],...
    ['Observations DMQC done (Total: ',num2str(sum(n_obs_peryear_DMQC)),')'], ...
    prof_included,'Location','southoutside')
% background color
set(gcf,'color','w');
% grid in y axis
ax = gca;
ax.YGrid = 'on';

% figure labels
for iyear = 1: length(year_names)
    text(iyear + 0.02, n_obs_peryear_DMQC(iyear) + 110,[num2str(round(n_obs_peryear_DMQC(iyear)/n_obs_peryear(iyear)*100,1)) ' %']);
end

% save figure
out_name = [pwd '/' working_date '/' project_name '_obs_DMQCstatus_byyear_' working_date '.png'];
export_fig(out_name)

%% %%%%%%%%%%%%%% Fig 6 : observation DMQC age histogram %%%%%%%%%%%%%%

figure(6)
% bigger figure
set(gcf, 'Position', [200, 200, 700, 600])
% figure name
set(gcf,'Name','Observations DMQC age histogram')

hndl = histogram(obs_age_notDM/365);

% FIGURE FORMAT
hold on
plot(ones(1,hndl.NumBins),hndl.Values,'--k','LineWidth',2.5); % for comment in legend
% xlabels
xlabel('Age of observations (years)')
% title with update date
title(['Age of observations with no DMQC (updated ',datestr(datenum(Output.index_update,'yyyymmddHHMMSS'),'yyyy-mm-dd'),')'])
ylabel('Number of observations')
% background color
set(gcf,'color','w');
% grid in y axis
ax = gca;
ax.YGrid = 'on';

% save figure
out_name = [pwd '/' working_date '/' project_name '_obs_DMQCstatus_agehist_' working_date '.png'];
export_fig(out_name)

%% Warnings

% floats with F cycles
Ffiles_WMO = cellstr(Output.WMO(Output.psal_qcF_prof > 0,:));
% floats > 1 year and never DMQC
neverDM_WMO = cellstr(Output.WMO(Output.percen_DMQC == 0 & Output.float_age > sage,:));
% floats grey list
greylist_WMO = cellstr(Output.WMO(Output.greylist,:));

% list all warnings
allwarning_WMO = unique(cellstr([Ffiles_WMO; neverDM_WMO; greylist_WMO]));
n_warfloats = length(allwarning_WMO);

if n_warfloats > 0
    [wmo_sorted, isort] = sort(cellstr(Output.WMO));
    dm_sorted = cellstr(Output.DM(isort,:));
%     warning_dm = Output.DM(contains(cellstr(Output.WMO), allwarning_WMO),:);
    warning_dm = dm_sorted(contains(wmo_sorted, allwarning_WMO),:);
    warning1_str = repmat(' ',n_warfloats,1);
    warning1_str(contains(allwarning_WMO, Ffiles_WMO)) = 'X';
    warning2_str = repmat(' ',n_warfloats,1);
    warning2_str(contains(allwarning_WMO, neverDM_WMO)) = 'X';
    warning3_str = repmat(' ',n_warfloats,1);
    warning3_str(contains(allwarning_WMO, greylist_WMO)) = 'X';

    % warning table
    warning_table = [allwarning_WMO, cellstr(warning_dm), cellstr(warning1_str), cellstr(warning2_str),cellstr(warning3_str)]';

    % header and disp (is table is not empty)
    fprintf(2,'\n\n--------------------------------- WARNINGS ---------------------------------\n')
    disp(' ')
    header_war = 'WMO    DMQC operator    At least 1cycle psal_QC= F     Never DMQC + >1 year    In grey list';
    fprintf('  %s \n', header_war);
    fprintf('%s%15s%15s%30s%20s\n',warning_table{:})
end
  

%% write text file with outputs
% TODO explain variables in headed and resume at the end with statistics in
% graphics (explain variables in pdf document)

outfile = [working_date '/' 'DMQCstatus_'  Output.index_update(1:8) '.txt'];
fprintf('\nSaving results in %s ...\n',outfile)

fid=fopen(outfile,'w');

% title and information

fprintf(fid, '# \n');
fprintf(fid, '# DMQC status and Data quality control statistics\n');

fprintf(fid, '# \n');

fprintf(fid, '# Project : %s\n', project_name);
fprintf(fid, '# Update date : %s\n', Output.index_update);
fprintf(fid, '# \n');
fprintf(fid, '# \n');

% Table 1: Statistics per float
fprintf(fid, '# Statistics per float\n');
% TODO explain variables
% header
header1 = ['WMO,' 'RT,' 'DM,' 'first_cycle_date,' 'DM_done,' 'float_age,'...
    'float_more1year,' 'greylist,' 'obs_number,' 'obs_more1year,'...
    'number_DMprof,' 'obs_more1year_noDM,' 'date_last_DMprof,' 'percentage_DMobs,'...
    'psal_qcA_prof,' 'psal_qcB_prof,' 'psal_qcC_prof,' 'psal_qcD_prof,'...
    'psal_qcE_prof,' 'psal_qcF_prof'];
fprintf(fid, '%s \n', header1);
% data
table1 = [cellstr(Output.WMO), cellstr(Output.RT), cellstr(Output.DM),...
          cellstr(Output.first_date), num2cell(Output.DM_done),...
          num2cell(Output.float_age), num2cell(Output.float_more1year),...
          num2cell(Output.greylist'), num2cell(Output.obs_number),...
          num2cell(Output.obs_more1year), num2cell(Output.number_DMprof),...
          num2cell(Output.obs_more1year_noDM), cellstr(Output.date_last_DMprof),...
          num2cell(Output.percen_DMQC), num2cell(Output.psal_qcA_prof'),...
          num2cell(Output.psal_qcB_prof'), num2cell(Output.psal_qcC_prof'),...
          num2cell(Output.psal_qcD_prof'), num2cell(Output.psal_qcE_prof'),...
          num2cell(Output.psal_qcF_prof')]';
fprintf(fid, '%s,%s,%s,%s,%d,%f,%d,%d,%d,%d,%d,%d,%s,%f,%d,%d,%d,%d,%d,%d\n',table1{:});


% Table 2: Statistics per DM responsible
fprintf(fid, '# \n');
fprintf(fid, '# Statistics per DM operator\n');
% TODO explain variables
% header
header2 = ['DM,' 'n_floats,' 'float_more1year,' 'float_DMdone,' 'n_obs,' 'obs_more1year,' ...
    'obs_DMdone,' 'obs_more1year_DMQC,'];
fprintf(fid, '%s \n', header2);
% data
table2 = [DM_operators, num2cell(number_floats'), num2cell(n_floats_more1year'),...
    num2cell(n_floats_DMQC'), num2cell(obs_number'), num2cell(n_obs_more1year'),...
    num2cell(n_obs_DMQC'), num2cell(n_obs_more1year_DMQC)']';
fprintf(fid, '%s,%d,%d,%d,%d,%d,%d,%d\n',table2{:});


% Table 3: Totals
fprintf(fid, '# \n');
fprintf(fid, '# Totals\n');
% TODO explain variables
% header
header3 = ['n_floats,' 'float_more1year,' 'float_DMdone,' 'floats_greylist,'...
    'n_obs,' 'obs_more1year,' 'obs_DMdone,' 'obs_more1year_noDM'];
fprintf(fid, '%s \n', header3);
% data
table3 = [num2cell(n_floats), num2cell(total_more1year), num2cell(total_DMfloats),...
    num2cell(total_greylist), num2cell(sum(obs_number)), num2cell(sum(n_obs_more1year)),...
    num2cell(sum(n_obs_DMQC)),num2cell(sum(Output.obs_more1year_noDM))]';
fprintf(fid, '%d,%d,%d,%d,%d,%d,%d,%d\n',table3{:});


% Table 4: Warnings
if n_warfloats > 0
    fprintf(fid, '# \n');
    fprintf(fid, '# Warnings\n');
    fprintf(fid, '  %s \n', header_war);
    fprintf(fid, '%s%15s%15s%30s%20s\n',warning_table{:});
end

fclose(fid);