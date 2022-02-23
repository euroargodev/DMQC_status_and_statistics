function [hdl] = split_figure(floats_per_fig, WMO, DATA)
% split figure data in lots and plot them in different figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% only inputs that came from find (matrix)
% the limit is 6 plots (6 colors)
%
% Modified on 20190528

%number of floats
n_floats = length(WMO);
% number of figures to be done
n_figs = floor(n_floats/floats_per_fig) + (mod(n_floats,floats_per_fig)>0);

% number of inputs ( plots per figure  x 2 with different color)
n_plots = length(DATA);

% TODO: if more than 6 plots
%colors
if n_plots == 8 % fig 1: psal_correction diferent from 0 or not
    % colors = {[1.0000 0.8516 0.7227],[1.0000 0.8516 0.7227],[0.8242 0.8242 0.8242],[0.8242 0.8242 0.8242],'b','b','r','r'};
    colors = {[0.7945 1.0000 0.7945],[0.3945 1.0000 0.3945],... % DM cycle not corrected
              [0.8242 0.8242 0.8242],[0.8242 0.8242 0.8242],... % RT cycle not corrected
              'b','b',... % Psal correction = 0
              'r','r'}; % Non zero psal correction
elseif n_plots == 12 % fig 3: profile_psal_qc
    colors = {'b','b',... % A
              'c','c',... % B
              'g','g',... % C
              'y','y',... % D
              'm','m',... % E
              'r','r'};   % F
elseif n_plots == 4 % fig 4: DMQC done
    colors = {'b','b',... % Real time cycle
              'g','g'}; % DMQC cycle
else % fig 5: Most frequently occuring QC flag
    colors = {[0.8242 0.8242 0.8242],[0.8242 0.8242 0.8242],... % 0
              'b','b',... % 1
              'g','g',... % 2
              [1.0000 0.5469 0],[1.0000 0.5469 0],... % 3
              'r','r',... % 4
              [0.5625 0.9297 0.5625],[0.5625 0.9297 0.5625],... % 5
              [0.4961 1.0000 0.8281],[0.4961 1.0000 0.8281],... % 8
              'k','k'}; % 9
end

ifloat = 0;
for ifig = 1 : n_figs % figures loop
    
    % getting last float of lot
    if (ifloat + floats_per_fig) > n_floats
        end_float = n_floats;
    else
        end_float = ifloat + floats_per_fig;
    end
    
    figure

    hold on
    for iplot = 1 :2: n_plots
        
        xData = DATA{iplot};
        yData = DATA{iplot+1};

        [fl,xi,~] = unique(yData);
        yindex = find(ismember(fl, ifloat+1:end_float));
        if isempty(yindex) % no data for this lot
            hdl{ifig,iplot} = scatter(-1,-1, 10, colors{iplot},'filled'); % for legend
            continue
        end
        
        float_start = xi(yindex(1));
        if xi(yindex(end)) == xi(end) % when last float in data array
            float_end = length(yData);
        else
            float_end = xi(yindex(end)+1)-1;
        end
        
        % plot
        hdl{ifig,iplot} = scatter(xData(float_start:float_end),yData(float_start:float_end),10,colors{iplot},'filled');
    end

     % y axis format
     yticks(ifloat+1:end_float)
     yticklabels(cellstr(WMO(ifloat+1:end_float,:)))
     ylim([ifloat,end_float+0.5])
     
      ifloat = ifloat + floats_per_fig;
end  % figures loop
