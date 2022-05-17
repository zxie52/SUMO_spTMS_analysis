close all;
% close the eeglab GUI and then clear the workspace, to avoid possible bug from eeglab
clear;

%% Start the Journey

% Generally, there will be a few sessions for the preprocessing:
% Step 1: Merge the beh data on to the raw eeg file(done on the CRC due to the size of raw eeg file)
% Step 2: Epoching the data file(done on the CRC due to the size of the raw eeg file)
% Step 3: Preprocessing the eeg files for the IEM anaylsis and ERP/ERSP analysis
% 
% After the preprocessing, we will go to the analysis
% This part is only for the TMS part of the EEG
% for other parts(i.e stim, cue, and probe, please use xxx_noTMS.m in this folder)
%
% This script requires EEGLAB v2020.0 in MATLAB

% The flow for this part:
%                               Loading EEG sets
%                                       |
%                    Remove pre-stimulus baseline to remove the DC offset
%                                       |
%                              Remove poor channels
%                                       |
%                     Downsample data (10,000 Hz to 1,000 Hz)
%                                       |                                       |
%                    Detrend the EEG.data for removing the remainng DC offset
%                                       |
%        High pass through 1 Hz before 1st round of ICA to remove the remaining DC offset
%                                       |
%                      Bandpass Filter(1-100Hz, excluding 55-65Hz)
%                                       |
%                   fastICA to remove other artifacts(e.g blinking)
%                                       |
%                Interpolate the missing data and removed poor channels
%                                       |
%                               Re-reference to median
%                                       |
%                                    Save out
%                                       
% Be careful: This script is only for preprocessing the no-TMS intervals
% Also, because the ICA needs interaction between frontend and user, it is highly recommended to run
% this preprocessing scipt locally

%% Loading basics

% before the start of the preprocessing, download the eeg files from CRC to local machine
fpath3 = 'E:\SUMO_further_data_pack_zx\N2pc_IEM\new_results\eeg_before_preproc';
fpath7 = 'E:\SUMO_further_data_pack_zx\N2pc_IEM\new_results\prac';
fpath8 = 'E:\SUMO_further_data_pack_zx\N2pc_IEM\new_results\prac';
% fpath7 = 'E:\SUMO_further_data_pack_zx\N2pc_IEM\new_results\eeg_before_2ndICA';
% fpath8 = 'E:\SUMO_further_data_pack_zx\N2pc_IEM\new_results\eeg_before_IEM';

subject = {'SUMO_0102', 'SUMO_0104', 'SUMO_0105', 'SUMO_0106',  ...
           'SUMO_0108', 'SUMO_0111', 'SUMO_0114', 'SUMO_0120',...
           'SUMO_3001', 'SUMO_3017', 'SUMO_3015'};
label = {'stim', 'cue1', 'probe1', 'cue2', 'probe2'};

% trigger label from EEG vmrk files
trigger = {'S 48', 'S 49', 'S 51', 'S 52', 'S 54'};

% choose the condition, cue1, cue2, probe1, probe2 or stimulus
warning("It is highly recommended to preprocess in order: stim-cue-probe epochs, to check if subject blinks during the onset of the cue.");
type  = input("What type of condition do you want to preprocess? stim, cue1, probe1, cue2 or probe2\n", 's');
while ~strcmp(type, label)
    warning("For tms intervals, please use the script xxx_spTMS.m");
    fprintf("Please type in the valid type: stim, cue1, probe1, cue2 or probe2 \n");
    type  = input("What type of condition do you want to preprocess? stim, cue1, probe1, cue2 or probe2 \n", 's');
end
%% Saving baseline for later tesa
for l = 1 : length(subject)
    %% Step 1: Load eeg files which were epoched around TMS 1
    cd(fpath3);
    EEG = pop_loadset('filename', strcat(subject{l}, '_after_epoch_', type, '.set'), 'filepath', fpath3);
    eeglab redraw;
    fprintf(strcat("Preprocessing for Subject: ", subject(l), ", has been started! \n"));
        
    %% Step 2: Downsample the Dataset to 1000Hz
    EEG = pop_resample(EEG, 1000);
    
    %% Step 3: Remove DC offset
    % still questionable to remove the DC offset
    EEG = pop_rmbase( EEG, [-200, -5]);
    eeglab redraw;

    %% Step 4: Remove Bad Channels
    EEG = pop_select( EEG,'nochannel',{'HEOG' 'VEOG' 'Photodiode'});
    % save the rest of the channel location 
    EEG.allchan = EEG.chanlocs;

    if strcmp(subject{l}, 'SUMO_0102')
        % has huge TMS artifcat on FT9
        EEG = pop_select( EEG,'nochannel',{'TP9'}); %e.g., 'FC1' 'TP9' 'TP9' 'FT8' 'T7'
    end
    if strcmp(subject{l}, 'SUMO_0105')
        % has huge TMS artifcat on FT9
        EEG = pop_select( EEG,'nochannel',{'P3'}); %e.g., 'FC1' 'TP9' 'TP9' 'FT8' 'T7'
    end
    if strcmp(subject{l}, 'SUMO_0108')
        % has huge TMS artifcat on FT9
        EEG = pop_select( EEG,'nochannel',{'TP9'}); %e.g., 'FC1' 'TP9' 'TP9' 'FT8' 'T7'
    end
    if strcmp(subject{l}, 'SUMO_0120')
        % has huge TMS artifcat on FT9
        EEG = pop_select( EEG,'nochannel',{'FT9'}); %e.g., 'FC1' 'TP9' 'TP9' 'FT8' 'T7'
    end
    if strcmp(subject{l}, 'SUMO_3001')
        % has huge TMS artifcat on F1 and F3
        EEG = pop_select( EEG,'nochannel',{'F3' 'F1'}); %e.g., 'FC1' 'TP9' 'TP9' 'FT8' 'T7'
    end

    %% Step 5: Filter data - from  TESA
    % Detrend the EEG.data for removing the remaining DC offset drifting
    for i = 1:EEG.trials, EEG.data(:,:,i) = detrend(EEG.data(:,:,i)')'; end
    eeglab redraw
    
    % Bandpass filter from 1 to 100Hz and exclude the 55-65Hz
    EEG = pop_tesa_filtbutter( EEG, 1, 100, 4, 'bandpass' ); %bandpass 1-100Hz like in Julkunen paper
    EEG = pop_tesa_filtbutter( EEG, 55, 65, 4, 'bandstop' ); %bandstop 55-65Hz like in Julkunen paper
    eeglab redraw;
    
    % Remove first 100ms and last 100ms to avoid the border effect after filtering
    EEG = pop_epoch( EEG, trigger(strcmp(type, label)), [-0.2 .8]);
    eeglab redraw;
    
    %% Extra Step:  Remove bad trials - from TESA (run whole chunk at a time)
    % meaning of values: (use data stored in EEG, use electrode data, look at all electrodes, locthresh > 5 SD, globalthresh > 5SD, do not reject but store marks instead)
    EEG = pop_jointprob(EEG, 1, 1:size(EEG.data,1),5,5,0,0);
    pop_rejmenu(EEG,1);
    pause_script = input('Highlight bad trials, update marks and then press enter');
    
    EEG = rmBadTr(EEG, type, fpath8, subject{l});
    eeglab redraw;
    
    % save out before the ICA
    cd(fpath7);
    EEG = pop_saveset(EEG, strcat(subject{l}, '_before_2ndICA_', type, '.set'));
    
    %% Step 7: Remove all other artifacts (using FastICA and auto component selection) - from  TESA
    EEG = pop_tesa_fastica( EEG, 'approach', 'symm', 'g', 'tanh', 'stabilization', 'off' ); %default should be off; if it faiils to converge turn it on
    EEG = pop_tesa_compselect( EEG,'compCheck','on','remove','on','saveWeights','off','figSize','medium',...
                               'plotTimeX',[-50 400],'plotFreqX',[1 100],'freqScale','log',...
                               'tmsMuscle','off','tmsMuscleThresh',8,'tmsMuscleWin',[-20 30],'tmsMuscleFeedback','off',...
                               'blink','on','blinkThresh',2.5,'blinkElecs',{'Fp1' 'Fp2'},'blinkFeedback','off',...
                               'move','on','moveThresh',2,'moveElecs',{'F7', 'F8'},'moveFeedback','off',...
                               'muscle','on','muscleThresh',-0.31,'muscleFreqIn',[7 70],'muscleFreqEx',[58 62],'muscleFeedback','off',...
                               'elecNoise','on','elecNoiseThresh',4,'elecNoiseFeedback','off' );
    %ICA_reject_count = sum(EEG.icaCompClass.TESA2.compClass ~= 1)

    %% Step 8: Interpolate missing channels
    EEG = pop_interp(EEG, EEG.allchan, 'spherical');

    %% Step 9: Re-reference to median
    EEG.data=bsxfun(@minus,EEG.data,median(EEG.data,1)); %median re-reference - better than average b/c average sensitive to outliers
    eeglab redraw
    
    %% Step 13: Save out the eeg dataset for later IEM analysis
    cd(fpath8);
    % save the EEG dataset for the later data analysis
    EEG = pop_saveset(EEG, strcat(subject{l}, '_before_iem_', type, '.set')); 
    fprintf(strcat("Preprocessing for Subject: ", subject(l), ", has been finished! \n"));
    eeglab redraw;
end 

fprintf(strcat("All preprocessings have been finished! \n"));

