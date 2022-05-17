clear;
close all;

%% Start the Journey

% Generally, there will be a few sessions for the preprocessing:
% Step 1: Merge the beh data on to the raw eeg file(done on the CRC due to the size of raw eeg file)
% Step 2: Epoching the data file(done on the CRC due to the size of the raw eeg file)
% Step 3: Preprocessing the eeg files for the IEM anaylsis and ERP/ERSP analysis
% 
% After the preprocessing, we will go to the analysis

% This script is only for SUMO111 - 120
% For the rest of SUMO 100 - 109, we used the already merged eeg file for
% processing (/afs/crc.nd.edu/group/roselab/EEG_fixed_merged)

addpath(genpath('/afs/crc.nd.edu/group/roselab/vol2/zx/new_results/analysis_code'));
addpath(genpath('/afs/crc.nd.edu/group/roselab/vol2/zx/matlab_envi/'));
fpath = '/afs/crc.nd.edu/group/roselab/vol2/EEG_raw/'; %path for raw EEG data
fpath2 = '/afs/crc.nd.edu/group/roselab/zx/beh/'; %path for beh data for merging
fpath3 = '/afs/crc.nd.edu/group/roselab/vol2/zx/eeg_data/eeg_before_preproc/'; 
%path for saving preproc EEG data

%for SUMO_0100 - 109, there are issues about the event marker. We then used the fixed_merged_EEG files direclty
subject_behavioral = {'SUMO_111','SUMO_114', 'SUMO_117', 'SUMO_118', 'SUMO_119', 'SUMO_120'};
subject = {'SUMO_0111', 'SUMO_0114', 'SUMO_0117', 'SUMO_0118', 'SUMO_0119', 'SUMO_0120'};
subjectnames = {'Subject 111', 'Subject 114', 'Subject 117', 'Subject 118', 'Subject 119', 'Subject 120'};
blocks = {8,8,7,9,7,7};
etypes = {'stim_trig', 'cue1_trig', 'TMS1_trig', 'probe1_trig', 'cue2_trig', 'TMS2_trig', 'probe2_trig'};
etypes3 = {'S 48', 'S 49', 'S 49', 'T  1', 'T  1', 'S 51', 'S 51'};
etypes4 = {'stimuli', 'cue1', 'TMS1', 'probe1', 'cue2', 'TMS2', 'probe2'};
epochs = {'epoch1', 'epoch2'};

eeglab;

%% Merging Beh data on the raw EEG files(SUMO 111 - 120)
% Main scripts are from Morgan -zx
for l = 1 : 6
    %create array to store the merged behavioral files we compile later
    mergedbehave = [];
    
    %load behavioral data
    cd(fpath2);
    load(char(strcat('DataWM_TMSEEG_',subject_behavioral{l},'.mat')));
    %transfer behavioral data from Data to mergedbehave

    %rename allorients
    for blk = 1:6 %blocks{l}
        Data(blk).eventslabel{14} = 'leftori'; %Data is a structure with each block as a field; reference the field with ()
        Data(blk).eventslabel{15} = 'rightori'; %to reference the text inside eventslabel(14) need to use {} i.e. eventslabel{15}
        
        %append each block to the end of mergedbehave (check that mergedbehave is 560 x 15 after this step)
        mergedbehave = [mergedbehave; Data(blk).events]; %this will vertically append the full behavioral data file for each block into one "mergedbehave" array (length = #trials, width = 15)
    end
    
    %% LOAD IN EEG DATA AND MERGE
    cd(fpath);

    EEG = pop_loadbv(fpath, (strcat(subject{l}, '.vhdr'))); 
    eeglab redraw
    
    EEG.event(strcmp({EEG.event.type}, 'empty')) = []; %this gets rid of random EEG markers that are empty
    EEG.event(strcmp({EEG.event.type}, 'boundary')) = []; %this gets rid of boundary markers
            
    % Clean up data(Frome Morgan's script: mannually delete some event marker problem trials in EEG)
    if l == 2 % also delete the first trial of the 4th Data block and trials 1-12 of the 5th Data block
        EEG.event(1:2) = [];
        EEG.event(1009) = [];
        EEG.event(1261:1265) = [];
        EEG.event(1513:1519) = [];
    elseif l == 4
        EEG.event(253:287) = [];
    elseif l == 5
        EEG.event(2017:2018) = [];
        EEG.event(3025:3028) = [];
    elseif l == 6
        EEG.event(1009:1011) = [];
        EEG.event(1261:1264) = [];
        EEG.event(2017:2022) = [];
        EEG.event(3025:3028) = [];
        EEG.event(3277:3279) = [];
    end

    %Merge the EEG and behavioral data 
    %%questionable -> if there is a ... faster way? this part takes much
    %%longer time than expected
    start = 1;
    stop = 9;
    for j = 1:length(mergedbehave) %for all trials
        for i = 1:size(Data(1).events, 2) %for each column of Data(1).events
           EEG = pop_editeventfield(EEG,'indices',start:stop,Data(1).eventslabel{i},mergedbehave(j,i)); %adds the label for the field and all of the merged behavioral data to the EEG.event structure
        end
        start = start + 9;
        
        stop = stop + 9;
    end
    
    cd(fpath3);
    EEG = pop_saveset(EEG, char(strcat(subject(l), '_before_preproc')), fpath3, '7.3');
    fprintf("Done ! \n");
end 
fprintf("All Done! Now move forward to preprocessing!\n");
% Now let's move to the next script: eeg_preproc_part2_Epoching.m
