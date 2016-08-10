%%% KB Shakman, 5/16/2016 %%%

% temp script, can modify for analysis of a given set of trials (organized
% in folders by odor).

clear all; close all; 

% add necessary dirs to search path: 
% addpath('/Users/katieshak/code/kimaging/gcamp_analysis_copies') 

% Check that 'names' is the same in both cells below. 

%% Inital analysis of raw data (navigate to desired parent folder first)

% To prespecify starting folder:
%raw_folder = '/Users/katieshak/Desktop/Data_Analysis_Temp/Raw_Unprocessed/2015_10_21_Live/br1_58B_6fA/p1_hemi1_left';

% To let user load starting folder:
raw_folder = uigetdir('','Select raw data folder');
cd(raw_folder);


% Get subfolder names (one for each odor/stim, generally)
folders = dir(); 
folders = folders([folders.isdir]);
folders = struct2cell(folders(arrayfun(@(x) x.name(1), folders) ~= '.')); % Ref: Steve Lord via http://www.mathworks.com/matlabcentral/newsreader/view_thread/258220
names = folders(1,:); % keep only the folder names

try
    tempFolder = names{1,1};
    cd(tempFolder);
    new_FluorAnalysis_batch_withPID_function(); 
    cd .. ; 
catch
    display('First folder not analyzed.');
    cd ..; % go back up a level
end

for nameIdx = 2:length(names)
    tempFolder = names{1,nameIdx};
    cd(tempFolder);
    new_FluorAnalysis_batch_withPID_function();
    cd ..;
end

display('Preprocessing done.');

%% Plots by group, for pre-processed data

findPattern = 'Raw_Unprocessed/'; newPattern = 'Data_Analysis/';
an_folder = strrep(raw_folder, findPattern, newPattern); 
% %To let user load an_folder:
% an_folder = uigetdir('','Select analysis folder');

cd(an_folder);
% try
%     cd(names{1,1});
%     PlotAll_DelFoverF_ofType;
% %     PlotAll_DelFoverF_ofType_newfunction();
% %     PlotAll_DelFoverF_ofType_new_function();
%     cd .. ;
% catch
%     an_folder = uigetdir('','Select analysis folder');
%     cd(an_folder);
% end

for nameIdx = 1:length(names)
    try
        currDir = names{1,nameIdx}; 
        cd(currDir);
        PlotAll_DelFoverF_ofType_func(currDir); 
        cd .. ;
    catch
       cd .. ; % skips over non-existent analysis dirs
    end

end


display('Analysis done.');
