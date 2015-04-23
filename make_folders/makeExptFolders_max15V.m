% Make folders for imaging expts
BASIC_NAME = 'TRIAL';
NUM_TRIAL_TYPES = 6;
for idx = 1:45 %50
    
    if idx > 5 && idx < 36
      
            if mod(idx,NUM_TRIAL_TYPES) == 0
                trialType = '_PT5V';
            elseif mod(idx,NUM_TRIAL_TYPES) == 1
                trialType = '_15V';
            elseif mod(idx,NUM_TRIAL_TYPES) == 2
                trialType = '_1V';
            elseif mod(idx,NUM_TRIAL_TYPES) == 3
                trialType = '_8V';
            elseif mod(idx,NUM_TRIAL_TYPES) == 4
                trialType = '_2V';
            elseif mod(idx,NUM_TRIAL_TYPES) == 5
                trialType = '_4V';
            end
    else
        trialType = '_BASELINE';
%     if idx > 59
%         trialType = '_BASELINE';
%     end
    end
    newFolderName = [BASIC_NAME, num2str(idx), trialType];
    mkdir(newFolderName);
    disp(newFolderName)
    disp('---')
end
% Cycle folders
for i=1:7
    if i==1
        newFolderName = 'PRE-EXPT_BASELINE';
    elseif i ==2
        newFolderName = 'CYCLE1';
    elseif i == 3
        newFolderName = 'CYCLE2';
    elseif i == 4
        newFolderName = 'CYCLE3';
    elseif i == 5
        newFolderName = 'CYCLE4';
    elseif i == 6
        newFolderName = 'CYCLE5';
    elseif i == 7
        newFolderName = 'POST-EXPT_BASELINE';
    end
    mkdir(newFolderName);
    disp(newFolderName)
    disp('---')
end