% Make empty figure to hold grid of data for comparison across odors and
% flies.

% clean up workspace
clear all; close all;
% addpath(genpath('/Users/katieshak/Desktop/Data_Analysis_Temp'));
% addpath(genpath('/Volumes/KATIELAB1/Data_Analysis')); 
%% adjustable parameters:
% odorDirs = {'1-oil','2-Oct','4-Vin', '5-OctVin', '6-H2O'};
% odorNames = {'Oil', 'Oct', 'Vin', 'OctVin', 'H2O'};
% odorDirs = {'1-oil','2-Oct','4-Vin', '5-OctVin', '6-H2O'};
% odorNames = {'Oil', 'Oct', 'Vin','OctVin', 'H2O'};
% odorDirs = {'1-oil', '2-Oct', '4-HiOVin','6-Aceto', '7-OilVin'};
% odorNames = {'oil', 'Oct', 'HiOVin','Aceto', 'OilVin'};
% odorDirs = {'1-oil', '2-Oct', '6-Aceto', '7-OilVin'};
% odorDirs = {'1-oil', '2-Oct', '4-HiWVin','5-H2O', '6-Aceto','7-WatVin'};
% odorNames = {'oil', 'Oct', 'HiWatVin','H2O','Aceto', 'WatVin'};
% odorDirs = {'1-oil', '2-Oct',  '6-Aceto'};
% odorNames = {'oil', 'Oct','Aceto'};
% odorDirs = {'1-oil', '2-Oct', '3-MCH','4-Far', '5-IpA', '6-HAc','7-Ben'};
% odorNames = {'oil', 'Oct', 'MCH','Far','IPA', 'HAC', 'Ben'};
% odorDirs = {'2-OctV1', '3-MCHV1','4-FarV1', '5-IpAV1', '6-HAcV1','7-BenV1'};
<<<<<<< HEAD
% odorDirs = {'2-OctMV1', '3-MCHMV1','4-FarMV1', '5-IpAMV1', '6-HAcMV1','7-BenMV1'};
% odorNames = {'Oct', 'MCH','Far','IPA', 'HAC', 'Ben'};
% odorDirs = {'2-OctMV1','4-FarMV1', '5-IpAMV1', '6-HAcMV1','7-BenMV1'};
odorDirs = {'2-Oct','4-Far', '5-IpA', '6-HAc','7-Ben'}; % {'2-OctV1','4-FarV1', '5-IpAV1', '6-HAcV1','7-BenV1'};
% odorDirs2 = {'2-Oct','4-Far', '5-IpA', '6-HAc','7-Ben'};
odorNames = {'Oct','Far','IPA', 'HAC', 'Ben'};
% odorDirs = {'1-oil','5-IpA'};
% odorNames = {'Oil','IPA'};
=======
odorDirs = {'2-Octv1','4-Farv1', '5-IpAv1', '6-HAcv1','7-Benv1'};
odorDirs2 = {'2-Oct','4-Far', '5-IpA', '6-HAc','7-Ben'};
odorNames = {'Oct','Far','IPA', 'HAC', 'Ben'};
>>>>>>> a19bfc0f32871c96b0b452be9bf63eda21c2bc04

numOdors = length(odorNames);
StartDir = pwd;
% addpath(genpath(StartDir));
yaxmin = -0.1; % -0.5 normally
<<<<<<< HEAD
yaxmax = 0.55; % 1.25 normally
=======
yaxmax = 0.9; % 1.25 normally
>>>>>>> a19bfc0f32871c96b0b452be9bf63eda21c2bc04
xaxmin = 0; % 0;
xaxmax = 10; % 8;
odorOn = 3;
odorOff = 5;
% note: user should also specify a file with the directories for each
% fly/hemi in that fly from which the user wants data imported.
<<<<<<< HEAD
% cd('/Users/katieshak/Desktop/Data_Analysis_Temp');
cd('/Volumes/KATIELAB1/FliesToAnalyze');% user specifies this folder
% groupFileString = 'DirsForImport_MV1.txt';
groupFileString = 'DirsForImport_KC-TeT_71D08L-6f.txt';
% groupFileString = 'MBONs_71D08LexA.txt' ;% user specifies this file 
% groupFileString = 'OilVSiPA.txt';
=======
cd('/Users/katieshak/Desktop/Data_Analysis_Temp');
% cd('/Volumes/KATIELAB1/FliesToAnalyze');% user specifies this folder
% groupFileString = 'DirsForImport_sorted_V1_vinegar_oct.txt';
groupFileString = 'DirsForImport_V1.txt' ;% user specifies this file 
>>>>>>> a19bfc0f32871c96b0b452be9bf63eda21c2bc04

fileID = fopen(groupFileString);
flyDirList = textscan(fileID,'%s');
fclose(fileID);
flyDirList;
cd(StartDir)
numFlies = length(flyDirList{1,1});
% add-ons for plotting etc.:
% addpath('/Users/katieshak/Documents/MATLAB/Extensions'); % contains subtitle.m

%% make empty figure:
gfig = figure;
for flyIdx = 1:numFlies
    for odorIdx = 1:numOdors*flyIdx
        subplot(numOdors,numFlies,odorIdx)
%         subplot(numOdors,numFlies,odorIdx);
    end
end

%% fill in figure
splotNum = 0; % initialize subplot counter
% % set up to fill in subplots by column (row is default):
plotindex = reshape(1:(numFlies*numOdors), numFlies, numOdors).';

% import data and add plots for each fly:

for flyIdx = 1:numFlies
    % get data directory for given fly:
    %     tempFlyDir = uigetdir(); % for manually choosing fly dirs
    tempFlyDir = flyDirList{1,1}{flyIdx,1}(1:end);
    
    cd(tempFlyDir);
    
    % fill in figure:
    for odorIdx = 1:numOdors
        splotNum = splotNum + 1; % iterate subplot #
        % plot all traces for given odor for this fly, plus average
        %         subplot(numFlies,numOdors,splotNum)
        subplot(numOdors, numFlies, plotindex(splotNum))
        hold on
        % grab the correct data
        %         odorDir = uigetdir(); % for manually choosing odor dirs
        odorDir = odorDirs{odorIdx};
        try
            try
                cd(odorDir);
            catch
                display('Odor directory not found, trying shorter version of folder name.');
                odorDir = odorDirs2{odorIdx};
                cd(odorDir)
            end
            % get time points
            tPtsFiles = dir('*timePoints.mat');
            tPtsDat = load(tPtsFiles(1).name);
            tPtsDat = tPtsDat.timePoints(1,:); % all in tPts should be the same.
            % get trace files
            tempTraceFiles = dir('*delFoverF.mat');
            tempTraces = nan(length(tempTraceFiles),length(tPtsDat)); % init
            tempAvgFile = dir('*avg_dFoF.mat');
            for traceIdx = 1:length(tempTraceFiles)
                tempTrace = load(tempTraceFiles(traceIdx).name);
                tempTrace = tempTrace.delFoverF;
                plot(tPtsDat,tempTrace);
                tempTraces(traceIdx,:) = tempTrace;
            end
            % add average
            avgTrace = mean(tempTraces,1);
            %         avgTrace = load(tempAvgFile.name);
            %         avgTrace = avgTrace.avg_dFoF;
            plot(tPtsDat,avgTrace,'Color','k','LineWidth',1.5);
            ylim([yaxmin, yaxmax]);
            xlim([xaxmin, xaxmax]);
            ylabel(odorNames{odorIdx});
            xlabel('Seconds');
            legend();
            % add odor on and off markers:
            lineOn = plot([odorOn,odorOn],[yaxmin,yaxmax]);
            lineOff = plot([odorOff,odorOff], [yaxmin,yaxmax]);
            % Make the new vertical lines into dash-dotted black lines.
            set(lineOn,'Color',[0.5 0.5 0.5],'LineWidth', 1, 'LineStyle','-.')
            set(lineOff,'Color',[0.5 0.5 0.5],'LineWidth', 1, 'LineStyle','-.')
            subtitle(groupFileString(15:end-4));
        catch
            display('Odor directory not found.');
            plot(NaN);
            %             tempTraces = nan(size(tempTraces)); % use length of last tempTrace
            
        end
        display('movin back up');
        cd('..'); % back up one level

    end
    display ('fly done:'); display(flyIdx);
    display('on to the next fly');
end % end for loop thru flies


cd(StartDir);
