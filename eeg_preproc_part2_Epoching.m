clear;
close all;

%% Start the Journey

% Generally, there will be a few sessions for the preprocessing:
% Step 1: Merge the beh data on to the raw eeg file(done on the CRC due to the size of raw eeg file)
% Step 2: Epoching the data file(done on the CRC due to the size of the raw eeg file)
% Step 3: Preprocessing the eeg files for the IEM anaylsis and ERP/ERSP analysis
% 
% After the preprocessing, we will go to the analysis

% In this part, I will preprocess in 2 for loops for epoching on qstat tasks
% One for SUMO 111-120, the other for SUMO 100-109
% After this script, we will epoch data on different taskes

% This script requires EEGLAB v2020.0 in MATLAB

% The flow for this part:
%                               Loading EEG sets
%                                       |
%                                Load 63 channels
%                                       |
%                              Add channel locations
%                                       |
%                     Remove photosensor channel and HEOG/VEOG channels                                                   |
%                             Remove Early cue trials
%                                       |
%                             Epoch on different tasks
%                                       |
%                                    Save out

%% Start the Preprocessing
%% For SUMO 100 - 109

addpath(genpath('/afs/crc.nd.edu/group/roselab/vol2/zx/new_results/analysis_code'));
addpath(genpath('/afs/crc.nd.edu/group/roselab/vol2/zx/matlab_envi/'));
% folder to extract the raw EEG data
fpath = '/afs/crc.nd.edu/group/roselab/EEG_fixed_merged/'; %path for raw EEG data
% folder to save the epoched dataset for later preprocessing
fpath3 = '/afs/crc.nd.edu/group/roselab/vol2/zx/eeg_data/eeg_before_preproc/';

subject = {'SUMO_0100_fixedEEG_mergedF', ...
           'SUMO_0104_fixedEEG_mergedF', ...
           'SUMO_0106_fixedEEG_mergedF', ...
           'SUMO_0105_mergedEEG_mergedF', ...
           'SUMO_0102_mergedF', ...
           'SUMO_0103_mergedF', ...
           'SUMO_0108_fixedEEG_mergedF', ...
           'SUMO_0109_mergedF'};
       
subject2 = {'SUMO_0100', 'SUMO_0104', 'SUMO_0106', 'SUMO_0105', ...
            'SUMO_0102', 'SUMO_0103', 'SUMO_0108', 'SUMO_0109'};

%% Preprocessing part
for l = 1 : 8
    % load the merged eeg sets from /afs
    cd(fpath);
    EEG = pop_loadset('filename', strcat(subject{l},'.set'), 'filepath', fpath);
    
    % For some subjects, there is only T 1 labels rather than T 2 labels
    % Fixing the beh labels on EEG.event, change the T  1 to T  2 labels
    for x = 1:length(EEG.event)
        if strcmp(EEG.event(x).type, 'T  1') && strcmp(EEG.event(x-1).type,'S 52')
            EEG.event(x).type = 'T  2'; %TMS2 is now T  2
        end
    end
    
    % Save out
    cd(fpath3);
    EEG = pop_saveset(EEG, char(strcat(subject2{l}, '_before_preproc')));
    fprintf(strcat(subject2{l}, " Labels have been fixed! \n"));
    
%% Start the preprocessing
    %% Step 1:Add channel locations
    % We can use the old ced file for the electrodes locations or use the eeglab's one
    EEG = pop_chanedit(EEG, 'lookup','/afs/crc.nd.edu/group/roselab/zx/matlab_envi/bv_chanlocs_new.ced');
    EEG.allchan = EEG.chanlocs;
    eeglab redraw;
    
    %% Step 2: Remove bad/unused channels
    % This section, we will remove the photosensor channel and HEOG/VEOG
    EEG = pop_select( EEG,'nochannel',{'Photodiode' 'HEOG' 'VEOG'}); %e.g., 'FC1' 'TP9' 'TP9' 'FT8' 'T7'
    eeglab redraw
    
    %% Step 3: Fixing trigger issues
    % Epoching eeg datasets around TMS. One epoch will be one trial
    EEG = pop_epoch( EEG, { 'T  1' }, [-4 8]);

    % In some trials, 'S 49'(cue) may occur really fast after 'S 48'(stimulus)
    % Then we need to find those trials out and remove them
    % In this case, we delete trials where the latency difference between S49
    % and S48 is less than 200ms
    latency_diff = [];
    delete_epoch = [];
    j = 1;
    k = 1;
    for i = 1: length(EEG.event) - 1
        if strcmp(EEG.event(i).type, 'S 48')
          latency_diff(j,1) = EEG.event(i).epoch;
          latency_diff(j,2) = EEG.event(i+1).latency - EEG.event(i).latency;
          if latency_diff(j,2) < 200 % if the time is less than 200ms
              delete_epoch(1,k) = latency_diff(j,1);
              k = k + 1;
          end
          j = j + 1;
        end
    end
    % Visually inspect the trial again before rejecting
    EEG = pop_rejepoch (EEG, delete_epoch, 0);  
    
    fprintf(strcat(subject2{l}," Channels have been referenced! \n"));

%% Epoching data based on different tasks
    % We will save seven epoched data sets:
    % Stimulus, Cue 1, TMS 1, Probe 1, Cue 2, TMS 2, Probe 2 
    %%%% This part can be reconsidered. We can epoch a little bit lager epoch, e.g from -500ms to
    %%%% +1000ms for later preprocessing
    EEG_Stim = pop_epoch( EEG, { 'S 48' }, [-0.3 .9]);
    % In order to have a smaller size of dataset, we downsample the dataset to 1kHz except TMS
    EEG_Stim = pop_resample(EEG_Stim, 1000);
    EEG_Stim = pop_saveset(EEG_Stim, char(strcat(subject2{l}, '_after_epoch_stim')));
    
    EEG_Cue1 = pop_epoch( EEG, { 'S 49' }, [-0.3 .9]);
    EEG_Cue1 = pop_resample(EEG_Cue1, 1000);
    EEG_Cue1 = pop_saveset(EEG_Cue1, char(strcat(subject2{l}, '_after_epoch_cue1')));
    
    EEG_TMS1 = pop_epoch( EEG, { 'T  1' }, [-0.3 .9]);
    EEG_TMS1 = pop_saveset(EEG_TMS1, char(strcat(subject2{l}, '_after_epoch_tms1')));
    
    EEG_Probe1 = pop_epoch( EEG, { 'S 51' }, [-0.3 .9]);
    EEG_Probe1 = pop_resample(EEG_Probe1, 1000);
    EEG_Probe1 = pop_saveset(EEG_Probe1, char(strcat(subject2{l}, '_after_epoch_probe1')));
    
    EEG_Cue2 = pop_epoch( EEG, { 'S 52' }, [-0.3 .9]);
    EEG_Cue2 = pop_resample(EEG_Cue2, 1000);
    EEG_Cue2 = pop_saveset(EEG_Cue2, char(strcat(subject2{l}, '_after_epoch_cue2')));
    
    EEG_TMS2 = pop_epoch( EEG, { 'T  2' }, [-0.3 .9]);
    EEG_TMS2 = pop_saveset(EEG_TMS2, char(strcat(subject2{l}, '_after_epoch_tms2')));
   
    EEG_Probe2 = pop_epoch( EEG, { 'S 54' }, [-0.3 .9]);
    EEG_Probe2 = pop_resample(EEG_Probe2, 1000);
    EEG_Probe2 = pop_saveset(EEG_Probe2, char(strcat(subject2{l}, '_after_epoch_probe2')));
    
end    

%% for SUMO 111 - 120
clearvars -except keep fpath3;
subject = {'SUMO_0111', 'SUMO_0114', 'SUMO_0117', 'SUMO_0118', 'SUMO_0119', 'SUMO_0120'};

for l = 1 : 6
    % Load EEG data from SUMO 111 to SUMO 120
    cd(fpath3);
    EEG = pop_loadset('filename', strcat(subject{l},'_before_preproc.set'), 'filepath', fpath3);
    
    %% Step 1:Add channel locations
    EEG = pop_chanedit(EEG, 'lookup','/afs/crc.nd.edu/group/roselab/zx/matlab_envi/bv_chanlocs_new.ced');
    EEG.allchan = EEG.chanlocs;
    eeglab redraw;
    
    %% Step 2: Remove bad/unused channels
    % This section, we will remove the photosensor channel and VEOG
    EEG = pop_select( EEG,'nochannel',{'HEOG' 'VEOG' 'Photodiode'}); %e.g., 'FC1' 'TP9' 'TP9' 'FT8' 'T7'
    eeglab redraw
    
    %% Step 3: Fixing Trigger issues
    % Epoching eeg datasets around TMS 1. One epoch will be one trial
    % Splitting into EEG1: around TMS 1; 
    % Splitting into EEG2: around TMS 2;
    EEG1 = pop_epoch( EEG, {  'T  1'  }, [-4  8], 'epochinfo', 'yes');
    EEG2 = pop_epoch( EEG, {  'T  2'  }, [-2  2], 'epochinfo', 'yes');
    
    % In some trials, 'S 49'(cue) may occur really fast after 'S 48'(stimulus)
    % Then we need to find those trials out and remove them
    % In this case, we delete trials where the latency difference between S49
    % and S48 is less than 200ms
    latency_diff = [];
    delete_epoch = [];
    j = 1;
    k = 1;
    for i = 1: length(EEG1.event) - 1
        if strcmp(EEG1.event(i).type, 'S 48')
          latency_diff(j,1) = EEG1.event(i).epoch;
          latency_diff(j,2) = EEG1.event(i+1).latency - EEG1.event(i).latency;
          if latency_diff(j,2) < 200 % if the time is less than 200ms
              delete_epoch(1,k) = latency_diff(j,1);
              k = k + 1;
          end
          j = j + 1;
        end
    end
    % Visually inspect the trial again before rejecting
    EEG1 = pop_rejepoch (EEG1, delete_epoch, 0);  
    EEG2 = pop_rejepoch (EEG2, delete_epoch, 0);
    % Save out eeg files
    cd(fpath3);
    EEG1 = pop_saveset(EEG1, char(strcat(subject{l}, '_before_epoch_epoch1')));
    EEG2 = pop_saveset(EEG2, char(strcat(subject{l}, '_before_epoch_epoch2')));
    fprintf(strcat(subject{l}," Channels have been referenced! \n"));
    
%% Epoching data based on different tasks
    % We will save seven epoched data sets:
    % Stimulus, Cue 1, TMS 1, Probe 1, Cue 2, TMS 2, Probe 2 
    EEG_Stim = pop_epoch( EEG1, { 'S 48' }, [-0.3 .9]);
    EEG_Stim = pop_saveset(EEG_Stim, char(strcat(subject{l}, '_after_epoch_stim')));
    
    EEG_Cue1 = pop_epoch( EEG1, { 'S 49' }, [-0.3 .9]);
    EEG_Cue1 = pop_saveset(EEG_Cue1, char(strcat(subject{l}, '_after_epoch_cue1')));
    
    EEG_TMS1 = pop_epoch( EEG1, { 'T  1' }, [-0.3 .9]);
    EEG_TMS1 = pop_saveset(EEG_TMS1, char(strcat(subject{l}, '_after_epoch_tms1')));
    
    EEG_Probe1 = pop_epoch( EEG1, { 'S 51' }, [-0.3 .9]);
    EEG_Probe1 = pop_saveset(EEG_Probe1, char(strcat(subject{l}, '_after_epoch_probe1')));
    
    EEG_Cue2 = pop_epoch( EEG2, { 'S 52' }, [-0.3 .9]);
    EEG_Cue2 = pop_saveset(EEG_Cue2, char(strcat(subject{l}, '_after_epoch_cue2')));
    
    EEG_TMS2 = pop_epoch( EEG2, { 'T  2' }, [-0.3 .9]);
    EEG_TMS2 = pop_saveset(EEG_TMS2, char(strcat(subject{l}, '_after_epoch_tms2')));
    
    EEG_Probe2 = pop_epoch( EEG2, { 'S 54' }, [-0.3 .9]);
    EEG_Probe2 = pop_saveset(EEG_Probe2, char(strcat(subject{l}, '_after_epoch_probe2')));
    
end
