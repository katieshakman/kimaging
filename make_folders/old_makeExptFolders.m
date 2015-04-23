% Make folders for imaging expts
BASIC_NAME = 'TRIAL';
NUM_TRIAL_TYPES = 6;
for idx = 1:40 %50
    
    for idx > 5 || idx < 36
        if mod(idx,NUM_TRIAL_TYPES) == 0
            trialType = 'PT5V';
        elseif mod(idx,NUM_TRIAL_TYPES) == 1
            trialType = '10V';
        elseif mod(idx,NUM_TRIAL_TYPES) == 2
            trialType = '1V';
        elseif mod(idx,NUM_TRIAL_TYPES) == 3
            trialType = '8V';
        elseif mod(idx,NUM_TRIAL_TYPES) == 4
            trialType = '2V';
        elseif mod(idx,NUM_TRIAL_TYPES) == 5
            trialType = '4V';
        end
        else
            trialType = 'BASELINE';
    end
    
    newFolderName = [BASIC_NAME, num2str(idx), trialType];
    mkdir(newFolderName);
    disp(newFolderName)
end