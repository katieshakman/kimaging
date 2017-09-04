% listall_avgAOdortrapz_inFolder
% Lists all avgAOdortrapz in folder

clear; close all; 

%% Get the folder and get names of the appropriate files in it + values in files

myFolder = uigetdir(); 
dirResults = dir('*AOCOdor_trapz.mat'); 
numResults = max(size(dirResults)); 
fileNames = cell(numResults,1); % initialize
values = nan(numResults,1); % initialize 

for idx = 1:numResults
    fileNames{idx} = dirResults(idx).name; 
    load(fileNames{idx}); % loads as: avgAOdortrapz
    values(idx) = avgAOdortrapz; 
end

%% Display the values and write them to a csv file. 
disp(values); 
csvwrite('avgAOdortrapz_values_csv.csv',values); 