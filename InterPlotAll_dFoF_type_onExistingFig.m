% InterplPlotAll_dFoF_ofType_onExistingFig

% Plot from previously stored deltaF/F data and timepoints. 
clear all; close all;

StartDir = pwd; 
useFigureFile = uigetfile('*.fig'); 
useFigH = openfig(useFigureFile); 

%% Load data files from several user-selected directories. 
imDir = uigetdir(); % Allow user to select directory. 
cd(imDir);
currentDirectory = pwd; % For use in save paths. 


all_delFoverF_files = dir('*delFoverF.mat'); % Nx1 matrix
all_timePoints_files = dir('*timePoints.mat'); % Nx1 matrix
all_baseline_files = dir('*baselineF0.mat'); % Nx1 matrix
if length(all_baseline_files) 
    mode = 1;
else    
    all_dFoF_files = dir('*dFoF.mat'); 
    all_tPts_files = dir('*tPts.mat'); 
    mode = 2;
    clear all_delFoverF_files all_timePoints_files all_baseline_files;
end

addpath(pwd); 
% Ask if user wants to add files from another folder:
openMore = 1; % Initialize openMore = 1; 
while openMore
    button = questdlg('Would you like to add files from another folder?',...
        'title','Yes','No','No');
    if isequal(button, 'No')
        openMore = 0;
    else
        imDir = uigetdir(); % Allow user to select directory.
        cd(imDir);
        currentDirectory = pwd; % For use in save paths.
        if mode == 1
            all_delFoverF_files = [all_delFoverF_files; dir('*delFoverF.mat')]; % Nx1 matrix
            all_timePoints_files = [all_timePoints_files; dir('*timePoints.mat')]; % Nx1 matrix
            all_baseline_files = [all_baseline_files; dir('*baselineF0.mat')]; % Nx1 matrix
        elseif mode == 2
            all_dFoF_files = [all_dFoF_files; dir('*dFoF.mat')];
            all_tPts_files = [all_tPts_files; dir('*tPts.mat')];
        end
    end
    addpath(pwd)
end

pattern = '._' ;
xq = 0:0.1:10;
yq = []; % initialize
if mode == 1
    % dFoF = cell(size(all_delFoverF_files));
    dFoF = [];
    tPts = [];
    F0s = [];
    for idx = 1:max(size(all_timePoints_files))
        notTrueMAT = strfind(all_timePoints_files(idx,1).name, pattern);
        if notTrueMAT
            ; % do nothing for now.  maybe delete these files in future.
        else
            load(all_timePoints_files(idx,1).name);
            tPts = [xq; tPts];
            %         tPts = [timePoints; tPts];
        end
    end
    clear idx;
    for idx = 1:max(size(all_delFoverF_files))
        %     new_delFoverF = load(all_delFoverF_files(idx,1).name);
        notTrueMAT = strfind(all_delFoverF_files(idx,1).name, pattern);
        if notTrueMAT
            display('passed, not true .mat'); % maybe delete these files in future.
        else
            load(all_timePoints_files(idx,1).name); % to get timePoints
            load(all_delFoverF_files(idx,1).name);
            display('loaded .mat file');
            %         for sub_idx = 1:min(size(tPts))
            tempY = interp1(timePoints,delFoverF,xq,'linear');
            yq = [tempY; yq];
            %         end
            %         dFoF(idx,1) = [delFoverF; dFoF];
        end
    end
    clear idx;
    for idx = 1:max(size(all_baseline_files)) % also get baselines
        notTrueMAT = strfind(all_baseline_files(idx,1).name, pattern);
        if notTrueMAT
            ; % do nothing for now.  maybe delete these files in future.
        else
            load(all_baseline_files(idx,1).name);
            F0s = [F_0; F0s];
        end
    end
    clear idx;
    clear pattern;
    
elseif mode == 2
    dFoF_all = [];
    tPts_all = [];
    F0s_all = [];
%      for idx = 1:max(size(all_tPts_files))
%         notTrueMAT = strfind(all_tPts_files(idx,1).name, pattern);
%         if notTrueMAT
%             ; % do nothing for now.  maybe delete these files in future.
%         else
%             load(all_tPts_files(idx,1).name);
%             tPts_all = [xq; tPts_all];
%             %         tPts = [timePoints; tPts];
%         end
%     end
    clear idx;
    for idx = 1:max(size(all_dFoF_files))
        %     new_delFoverF = load(all_delFoverF_files(idx,1).name);
        notTrueMAT = strfind(all_dFoF_files(idx,1).name, pattern);
        if notTrueMAT
            display('passed, not true .mat'); % maybe delete these files in future.
        else
            load(all_tPts_files(idx,1).name); % loads as tPts, which is M(trials)xT(timepoints)
            load(all_dFoF_files(idx,1).name); % loads as avg_dFoF
            display('loaded .mat file');
            %         for sub_idx = 1:min(size(tPts))
            tempY = interp1(tPts(1,:),avg_dFoF,xq,'linear');
            yq = [tempY; yq];
            %         end
            %         dFoF(idx,1) = [delFoverF; dFoF];
        end
    end
    clear idx;
    clear pattern;
end
%% Plot GCaMP (for each trial and on average) versus time.
% fig = figure;
xvalsOn = [3,3];
xvalsOff = [5,5];
% subplot(1,2,1); % trial-by-trial plotting.
% hold all
% ylim([-.5 1.5])
% evokedTrialSums = zeros(max(size(yq)),1); % initialize
% if exist('yq', 'var')
%     plotTimes = xq;
%     
%     for idx = 1:min(size(yq))
%         plot(tPts(1,:), yq(idx,:));
%         %
%         ind = 0;
%         for i = 1:length(plotTimes)
%             if (plotTimes(i) > 3) && (plotTimes(i) < 5)
%                 ind = ind +1;
%                 evokedTrialSums(idx) = evokedTrialSums(idx) + yq(idx, i);
%                 denom = ind;
%             end
%         end
%     end
%     evokedTrialAvgs = evokedTrialSums/denom;
%     %
%     legend('show')
%     avg_dFoF = mean(yq,1);
%     max_dFoF = max(yq);
%     min_dFoF = min(yq);
% end
% title('deltaF/F vs Time')
% % Annotate the point (3, delFoverF(3))
% % text(3, delFoverF(3),'\leftarrow (delFoverF(3sec))',...
% %      'HorizontalAlignment','left')
% % title('GCaMP Average vs Time')
% % % Add a line at timepoint = 3 sec
%     yvals = get(gca, 'ylim');
%     max_y_vals = yvals;
%     lineOn = plot(xvalsOn, yvals);
%     lineOff = plot(xvalsOff, yvals);
%     % Make the new vertical line at 3s a dash-dotted black line. 
%     set(lineOn,'Color','k','LineWidth', 1, 'LineStyle','-.')
%     set(lineOff,'Color','k','LineWidth', 1, 'LineStyle','-.')
% hold off
% 
% % % Also plot the mean of all traces:
% subplot(1,2,2);
avgTrace = mean(yq); 
% avg_dFoF_plot = plot(xq, avgTrace);
% set(avg_dFoF_plot, 'Color', 'blue', 'LineWidth', 1)
% 
% hold on
% ylim([-.5 1.5])
% title('Average deltaF/F vs Time')
% 
%     set(gca, 'ylim', max_y_vals);
%     yvals = get(gca, 'ylim');
% %     xvalsOn = [3,3];
%     lineOn = plot(xvalsOn, yvals);
% %     xvalsOff = [4,4];
%     lineOff = plot(xvalsOff, yvals);
%     % Make the new vertical line at 3s a dash-dotted black line. 
%     set(lineOn,'Color','k','LineWidth', 1, 'LineStyle','-.')
%      % Make the new vertical line at 4s a dash-dotted black line.
%     set(lineOff,'Color','k','LineWidth', 1, 'LineStyle','-.')
%     legend
% hold off

% % Replot the mean of all traces as a new figure: 
useFigH; 
hold on; 
avg_dFoF_plot = plot(xq, avgTrace);
set(avg_dFoF_plot, 'Color', 'r', 'LineWidth', 5) % [0 0 1] for blue
% hold on
title('Average deltaF/F vs Time')
yvals = get(gca,'ylim'); 
    lineOn = plot(xvalsOn, yvals);
    lineOff = plot(xvalsOff, yvals);
    % Make the new vertical line at 3s a dash-dotted black line. 
    set(lineOn,'Color','k','LineWidth', 1, 'LineStyle','-.')
     % Make the new vertical line at 4s a dash-dotted black line.
    set(lineOff,'Color','k','LineWidth', 1, 'LineStyle','-.')
    legend
hold off
% Optionally: load the avg_dFoF files in the folder and plot those as
% dashed lines on top of the average trace for this whole group: 
hold on 
flyAvgFiles = dir('*avg_dFoF.mat'); 
flyAvgFiles = {flyAvgFiles.name}'; 
flyAvgs = {}; 
flyTptsFiles = dir('*tPts.mat'); 
flyTptsFiles = {flyTptsFiles.name}';
for fly = 1:length(flyAvgFiles)
    load(flyAvgFiles{fly}); % loads as avg_dFoF
    flyAvgs{fly} = avg_dFoF; 
    load(flyTptsFiles{fly}); % loads as tPts
    tempTpts = tPts(1,:); 
    flyAvgsInterp{fly} = interp1(tempTpts,flyAvgs{fly},xq,'linear');
    plot(xq,flyAvgsInterp{fly},'--','Color','r'); 
end
hold off
whitebg('k');


%% Save figure files and evokedAvg value
% % % dateStart = strfind(imDir, '201'); 
% % % dateEnd = dateStart-1 + min(strfind(imDir(dateStart:end), '/'))-1; 
% % % date = imDir(dateStart:dateEnd);
% % % trialStart = max(strfind(imDir, '/'))+1; 
% % % trialID = imDir(trialStart:end); 
% % % figName = [date,'_', trialID];
% % % %savePath = currentDirectory(...
% % %             %strfind(currentDirectory, 'Raw_Un')+16:strfind(...
% % %             %currentDirectory, 'TSeries')-1);
% % % %figName = savePath; 
% % % saveas(fig, figName, 'fig'); 
% % % saveas(fig, figName, 'png'); 
% % % 
% % % ind = 0; 
% % % 
% % % for i = 1:length(timePoints)
% % %     if (timePoints(i) > 3) && (timePoints(i) < 5)
% % %         ind = ind +1;
% % %         evokedTimes(ind) = timePoints(i);
% % %         evokedResp(ind) = avg_dFoF(i);
% % %     end
% % % end
% % % 
% % % % convert currentDirectory to currentDirName
% % % currentDirName = strrep(currentDirectory, '/', '_')  % Print to console.
% % % 
% % % evokedAvg = mean(evokedResp); 
% % % evokedAvg_savefile = [figName, '_evokedAvg']; 
% % % save(evokedAvg_savefile, 'evokedAvg'); 
% % % load(evokedAvg_savefile); 
% % % %
% % % evokedPeak = max(evokedResp); 
% % % evokedPeak_savefile = [figName, '_evokedPeak']; 
% % % save(evokedPeak_savefile, 'evokedPeak'); 
% % % load(evokedPeak_savefile); 
% % % %
% % % avgBaseline = mean(F0s); 
% % % avgBaseline_savefile = [figName, '_avgBaseline'];
% % % save(avgBaseline_savefile, 'avgBaseline'); 
% % % load(avgBaseline_savefile); 
% % % 
% % % %display('evokedAvg ='); 
% % % evokedAvg
% % % %display('evokedPeak ='); 
% % % evokedPeak
% % % %
% % % avgBaseline
% % % %
% % % evokedTrialAvgs % display values
% % % 
% % % 
% % % avg_dFoF_savefile = [figName, '_avg_dFoF'];
% % % save(avg_dFoF_savefile, 'avg_dFoF'); 
% % % 
% % % % save a copy in the GroupedData folder)
% % % savepath = '/Volumes/KATIELAB1/Data_Analysis/GroupData/' ;
% % % cd(savepath);
% % % save(evokedAvg_savefile, 'evokedAvg');
% % % save(avg_dFoF_savefile, 'avg_dFoF', 'tPts'); 
% % % % save(fullfile(savepath,evokedAvg_savefile,'evokedAvg'));
% % % % save(fullfile(savepath,avg_dFoF_savefile,'avg_dFoF'));
% % % %save(avg_dFoF_savefile, 'avg_dFoF'); 
% % % 
% % % cd(StartDir); 