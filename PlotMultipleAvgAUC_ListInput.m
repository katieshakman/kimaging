% Plot Multiple AvgZscores

clear; 
close all; 
cd('/Volumes/KATIELAB1/Data_Analysis/B-1_Stim/MB543Bstim_B-1_Live')
%% Set labels and user parameters
xLab = ''; 
yLab = 'AUC'; 

%% Choose data to load

group = 'expt_post_Mch'; 
files = MB0543BcV1gLive_loadDataFor_PLotMultipleAvgAUC_ListInput(group);
textname = 'MB0543BcV1gLive'; 

%% Load data and process data and make first plot
figure; 
traces = cell(length(files),1); 
hold on; 
flyLabels = cell(1,length(files)+1); 
for filei = 1:length(files)
    filename = files{filei,1};
    load(filename); % loads and timePoints and avgZscores
    disp(strcat('Loaded:',filename)); 
    AUCs(filei,:) = Atrapz; 
    avgAUCs(filei,1) = mean(Atrapz); 
    flyLabels{filei} = ['Fly', num2str(filei)];
end
    scatter(1:length(avgAUCs),avgAUCs); 

%% Plot the AUCs for each fly and the average
figure; 
try
    scatter(1:length(avgAUCs), avgAUCs,'LineWidth',5,'Color','k'); 
    plot(mean(1:length(avgAUCs)), mean(avgAUCs), 'LineWidth',5,'Color','b'); 
    xlabel(xLab);  ylabel(yLab);
    title('AUC For Each Fly and Average Across Flies');
catch
    disp('Something went wonky!'); 
end
flyLabels{end} = 'Average'; 
lgdMany = legend(flyLabels);

%% Save avg and SEM and max for each fly (after shuttered values removed). 
saveName = [textname,'_',group,'_avgAUCs']; 
save(saveName, 'avgAUCs'); 
display(pwd)
