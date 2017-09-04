% interpolate_and_plot 
% Katie Shakman - 3/03/2017
%% Setup
clear; close all; 
addpath(genpath('/Users/katieshak/code/')); % on macbook pro
addpath(genpath('/Users/katherineshakman/code')); % on imac
%% Loading Data and Interpolation
dfFiles = uigetfile('*avg_dF*.mat','MultiSelect', 'on'); 
timeFiles = uigetfile('*tPts*.mat','MultiSelect', 'on');
interpTimepts = 0:0.05:10; % interpolates to get 200 pts from 1 to 10
for fileIdx = 1:length(dfFiles)
    try
        load(dfFiles{fileIdx}); % Loads as avg_dFoF  
        load(timeFiles{fileIdx}); % Loads as tPts, size trialsxtimes.
    catch % If only one of each type was selected.
        load(dfFiles); % Loads as avg_dFoF.
        load(timeFiles); % Loads as tPts, size trialsxtimes.
    end
    df{:,fileIdx} = avg_dFoF;
    times{:,fileIdx} = tPts(1,:);
    assert(length(df{:,fileIdx}) == length(times{:,fileIdx})); 
    
    dfInterp(:,fileIdx) = interp1(...
        times{:,fileIdx}, df{:,fileIdx},interpTimepts); % Interpolate.
end
dfAvgTseries = mean(dfInterp,2); 
%% make figure(s)
figH = figure; 
hold on
for fileIdx = 1:length(dfFiles) 
    plot( times{:,fileIdx}, df{:,fileIdx},'b--','LineWidth',3 );
end
plot(interpTimepts, dfAvgTseries,'b-','LineWidth',5); 
hold off

%save_figure_perhaps()
% Maybe the user wants to save the figure. 
% titleStr = input('Type a name for this figure','s');
% title(titleStr);
% savefig(titleStr); 