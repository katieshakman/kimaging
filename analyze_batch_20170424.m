%%% KB Shakman, 5/16/2016 %%%

% temp script, can modify for analysis of a given set of trials (organized
% in folders by odor).

% add necessary dirs to search path: 
% addpath('/Users/katieshak/code/kimaging/gcamp_analysis_copies')
% %forLaptop
% OR
% addpath(genpath('/Users/katherineshakman/code/kimaging/gcamp_analysis_copies'))
% % foriMacDesktop
% 
% May also need to add the folders where the data will be, e.g.:
% addpath(genpath('/Volumes/KATIELAB1/Raw_Unprocessed/'))
% addpath(genpath('/Volumes/KATIELAB1/Data_Analysis/'))

 clear all; close all; % if not using as a function

% Check that 'names' is the same in both cells below. 

%function [] = analyze_batch(raw_folder)

%% Inital analysis of raw data (navigate to desired parent folder first)

% To prespecify starting folder:
%raw_folder = '/Users/katieshak/Desktop/Data_Analysis_Temp/Raw_Unprocessed/2015_10_21_Live/br1_58B_6fA/p1_hemi1_left';

% % To let user load starting folder:
% raw_folder = uigetdir('','Select raw data folder');
% cd(raw_folder);
% OR
if ~exist('raw_folder','var')
    raw_folder = pwd; 
end

% Get subfolder names (one for each odor/stim, generally)
folders = dir(); 
folders = folders([folders.isdir]);
folders = struct2cell(folders(arrayfun(@(x) x.name(1), folders) ~= '.')); % Ref: Steve Lord via http://www.mathworks.com/matlabcentral/newsreader/view_thread/258220
names = folders(1,:); % keep only the folder names
% names = { '1-oil'   , '2-octCSplus'  ,  '6-aceto'};

% Choose ROI for this analysis batch if desired: 

chooseROI = 0; % Set to 0 for no ROI choice, or 1 to use the drawROI function. 

if chooseROI == 1
    [ROI, roiName] = drawROI();
else
    ROI = NaN; roiName = [];
end

% try
%     tempFolder = names{1,1}; % tempFolder is the stimulus type folder (e.g. 1-oil)
%     cd(tempFolder);
%     new_FluorAnalysis_batch_withPID_function(ROI, roiName, tempFolder); 
%     cd ..; % go back up a level   
% catch
%     display('First folder not analyzed.');
%     cd(raw_folder); % go back to starting folder
%     display('Moved back to batch folder.'); 
% end

for nameIdx = 2%:length(names)
    tempFolder = names{1,nameIdx}; % tempFolder is the stimulus type folder (e.g. 2-Oct)
    cd(raw_folder); % go into the experiment folder (e.g. p1_hemi1)
    cd(tempFolder); % go into the stimulus type folder (e.g. 2-Oct)
    new_FluorAnalysis_batch_withPID_function(ROI, roiName, tempFolder);
    cd ..; % back into experiment folder (e.g. p1_hemi1)
end

display('Preprocessing done.');

%% Plots by group, for pre-processed data

% findPattern = 'Raw_Unprocessed/'; newPattern = 'Data_Analysis/';
% % an_folder = strrep(raw_folder, findPattern, newPattern); 
% %To let user load an_folder (starting from expected an_folder):
% analysis_path = strrep(raw_folder, findPattern, newPattern); 
% % if ~exist(char(analysis_path),'dir')
% %     mkdir(analysis_path); 
% % end
% an_folder = uigetdir(analysis_path,'Select analysis folder');
% 
% cd(an_folder);
% % try
% %     cd(names{1,1});
% %     PlotAll_DelFoverF_ofType;
% % %     PlotAll_DelFoverF_ofType_newfunction();
% % %     PlotAll_DelFoverF_ofType_new_function();
% %     cd .. ;
% % catch
% %     an_folder = uigetdir('','Select analysis folder');
% %     cd(an_folder);
% % end
% 
% % Get subfolder names (one for each odor/stim, generally)
% folders = dir(); 
% folders = folders([folders.isdir]);
% folders = struct2cell(folders(arrayfun(@(x) x.name(1), folders) ~= '.')); % Ref: Steve Lord via http://www.mathworks.com/matlabcentral/newsreader/view_thread/258220
% names = folders(1,:); % keep only the folder names
% 
% for nameIdx = 1:length(names)
%     
%     try 
%         currDir = names{1,nameIdx}; 
%         cd(currDir);
%         display(currDir);
%         PlotAll_DelFoverF_ofType_func(); 
%         cd .. ;
%     catch
%        cd .. ; % skips over non-existent analysis dirs
%        display('PlotAll failed.  Will try with an argument.');
%     end
%     
%     try
%         currDir = names{1,nameIdx}; 
%         cd(currDir);
%         PlotAll_DelFoverF_ofType_func(currDir); 
%         cd .. ;
%     catch
%        cd .. ; % skips over non-existent analysis dirs
%        display('Directory skipped in analysis.');
%     end
% 
% end
% 
% display('Analysis done.');
% %end