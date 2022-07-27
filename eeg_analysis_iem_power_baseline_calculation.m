% Calculate the pre-stimulus baseline for alpha and beta power
addpath(genpath('/afs/crc.nd.edu/group/roselab/vol2/zx/matlab_envi/'));
fpath8 = 'E:\SUMO_further_data_pack_zx\N2pc_IEM\new_results\eeg_before_IEM';
% There will be 11 subjects available for the SUMO spTMS
subject = {'SUMO_0102', 'SUMO_0104', 'SUMO_0105', 'SUMO_0106',  ...
           'SUMO_0108', 'SUMO_0111', 'SUMO_0114', 'SUMO_0120', ...
           'SUMO_3001', 'SUMO_3017', 'SUMO_3015'};
        
type = {'stim', 'cue1', 'tms1', 'probe1', 'cue2', 'tms2', 'probe2'};

for l = 1 : length(subject)
    for t = 1
        cd(fpath8);
        EEG = pop_loadset('filename', strcat(subject{l}, '_before_iem_', type{t}, '.set'), 'filepath', fpath8);
     
        % Calculate the pre-stimulus baseline for theta power
        % see the dothewave.m/ or dothewave_broadband.m
        data = double(EEG.data);
        [pow, ~, ~, dstimes, freqs] = dothewave_broadband(data, 1000, 4, 1, [], EEG.times);
        % Only calcaulate the pre-stimulus baseline from -200ms to -2ms
        % Average across timepoints and trials
        pow_base = mean(mean(pow(:,:,find(EEG.times == -200):find(EEG .times == -2),:),3), 4);
        % Saveout
        cd(fpath8);
        save(strcat("Pre_stimulus_broadband_power_baseline", subject(l)), 'pow_base', 'dstimes');
    end
end
