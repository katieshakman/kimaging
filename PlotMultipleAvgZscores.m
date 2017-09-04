% Plot Multiple AvgZscores

clear; close all; 
%% Set labels
xLab = 'Seconds'; 
yLab = 'Z-Score'; 

%% Load and process data and make first plot

files = dir('2*score*'); 

figure; 
timeLen = 0; % initialize
longestTimes = cell(1,1); 
traces = cell(length(files),1); 
times = cell(length(files),1); 
hold on; 
for filei = 1:length(files)
    load(files(filei).name); % loads and timePoints and avgZscores
    if length(timePoints)>timeLen
        timeLen = length(timePoints); 
        longestTimes = timePoints;
    end
    traces{filei,1} = avgZscores; 
    times{filei,1} = timePoints; 
    plot(timePoints,avgZscores,'LineWidth',5); 
end

tracesNew = nan(length(traces),length(longestTimes)); 
for filei = 1:length(files)
    tracesNew(filei,:) = interp1(times{filei,1},traces{filei,1},longestTimes); 
end

avgZscoreManyFlies = mean(tracesNew,1); 
plot(longestTimes,avgZscoreManyFlies,'LineWidth',10,'Color','w'); 
xlabel(xLab);  ylabel(yLab); 
title('Trace For Each Fly and Average Across Flies'); 
lgd = legend('1','2','3','4','Average');
lgd.Title.String = 'Flies';

%% Make a figure with the mean and SEM error bars for all the traces in this group/folder: 
ntrials = size(tracesNew,1); 
SDdfEachTimepoint = std(tracesNew,1);
SEMdfEachTimepoint = SDdfEachTimepoint./sqrt(ntrials); 
figure; 
plot(longestTimes,avgZscoreManyFlies,'Color','w','LineWidth',3);
hold on
errorbar(longestTimes,avgZscoreManyFlies,SEMdfEachTimepoint,'Color','w','LineWidth',3); 
hold off
xlabel(xLab); ylabel(yLab);
title('Average and SEM Across Flies'); 