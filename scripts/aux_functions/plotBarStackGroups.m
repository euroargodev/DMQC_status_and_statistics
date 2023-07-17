function [h] = plotBarStackGroups(stackData, groupLabels)
%% Plot a set of stacked bars, but group them according to labels provided.
%
% Params: 
%      stackData is a 3D matrix (i.e., stackData(i, j, k) => (Group, Stack, StackElement)) 
%      groupLabels is a CELL type (i.e., { 'a', 1 , 20, 'because' };)
%
% Version: V2.01
%
% Copyright 2011 Evan Bollig (bollig at scs DOT fsu ANOTHERDOT edu
% V2.01 : 2023/07/07 update by Delphine Dobler to tackle issue with NumGroupsPerAxis == 1.
%
%% 
NumGroupsPerAxis = size(stackData, 1);
NumStacksPerGroup = size(stackData, 2);


% Count off the number of bins
groupBins = 1:NumGroupsPerAxis;
MaxGroupWidth = 0.65; % Fraction of 1. If 1, then we have all bars in groups touching
groupOffset = MaxGroupWidth/NumStacksPerGroup;
%figure
    hold on; 
for i=1:NumStacksPerGroup

    Y = squeeze(stackData(:,i,:));
    
    % Center the bars:
    
    internalPosCount = i - ((NumStacksPerGroup+1) / 2);
    
    % Offset the group draw positions:
    groupDrawPos = (internalPosCount)* groupOffset + groupBins;
    
    % D.D: There is an issue when NumGroupsPerAxis = 1
    % in 2019 and higher Matlav version, one can use
    % if NumGroupsPerAxis == 1 
    %    h(i,:) = bar(1, Y, 'stacked');
    % end
    % for lower release: use of a workaround:
    if NumGroupsPerAxis > 1 
        h(i,:) = bar(Y, 'stacked');
    else
        h(i,:) = bar([1;nan], [ Y'; nan(size(Y'))], 'stacked');
    end
    set(h(i,:),'BarWidth',groupOffset);
    set(h(i,:),'XData',groupDrawPos);
    
end
hold off;
set(gca,'XTickMode','manual');
set(gca,'XTick',1:NumGroupsPerAxis);
set(gca,'XTickLabelMode','manual');
set(gca,'XTickLabel',groupLabels);
end 
