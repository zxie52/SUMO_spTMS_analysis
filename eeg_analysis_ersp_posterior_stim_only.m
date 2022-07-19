%% Load basics before the ersp preprocessing
close all
clear;

% using CRC
addpath(genpath('/afs/crc.nd.edu/group/roselab/vol2/zx/new_results/analysis_code'));
addpath(genpath('/afs/crc.nd.edu/group/roselab/vol2/zx/matlab_envi/'));
fpath8 = '/afs/crc.nd.edu/group/roselab/vol2/zx/new_results/eeg_before_IEM';
% this code is only for the stim epoch
fpath12 = '/afs/crc.nd.edu/group/roselab/vol2/zx/new_results/ERSP_results_stim_epoch';

subject = {'SUMO_0102', 'SUMO_0104', 'SUMO_0105', 'SUMO_0106',  ...
           'SUMO_0108', 'SUMO_0111', 'SUMO_0114', 'SUMO_0120',...
           'SUMO_3001', 'SUMO_3017', 'SUMO_3015'};
type = {'stim', 'cue1', 'tms1', 'probe1', 'cue2', 'tms2', 'probe2'};

ROI = {'Pz', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', ...
       'POz', 'PO3', 'PO4', 'PO7', 'PO8', ...
       'Oz', 'O1', 'O2'};

%% Load subjects' eeg for preprocessing and FFT
for t = 1
    for s = 1 : length(subject)
        % load eeg data
        cd(fpath8);
        EEG = pop_loadset('filename', strcat(subject{s}, '_before_iem_', type{t}, '.set'), 'filepath', fpath8);
        
        % get the meta data for the eeg
        p = struct2table(EEG.event); % change the EEG.event to a cell array and then categorize into different bins 
 
        % for stimulus epoch, select all the trials
        % get left trial and right trial numebers 
        trial = p.epoch;

        % only pick up midline frontal channels
        impchan = find(ismember({EEG.chanlocs.labels}, ROI)); 
        data = EEG.data(impchan,:,:);
        
        % FFT, from Jason's code
        [pow, ~, ~, dstimes, freqs] = dothewave(data, 1000, [2 50], 49, 4, (1), [-200 -20], EEG.times);

        % calculate the mean across the trials, then put into the 4-d tensor
        pow = squeeze(mean(pow(:,:,:,trial),4));
        pow_group(:,:,:,s) = pow; 
        
        %%%% pow_group: channels x nfreqs x time x subject
    end
    
    %% Preare for the visualization    
    % shorten the time range from -150ms to 700ms to avoid the border effect from FFT
    dstimes2 = find(dstimes==-150):find(dstimes==700);
    dstimes3 = dstimes(dstimes2);
    
    %% Plotting the ersp for midline frontal for all trials
    cd(fpath12);
    % calculate the average across all channels and all trials
    ersp_group_mean = squeeze(mean(mean(pow_group,1),4));
    imagesc(dstimes3, freqs, ersp_group_mean(:,dstimes2));
    ersp_fig_setup();
    title(strcat("Posterior Channels on All Trials on ", type{t}), 'Fontsize', 25);
    savefig(strcat("allsubjects_ersp_", type{t},"_heatmap"));
    saveas(gcf, strcat("allsubjects_ersp_", type{t},"_heatmap.png"));
    saveas(gcf, strcat("allsubjects_ersp_", type{t},"_heatmap.svg"));
    
    % cluster permutation with FDR corrected t-test between channel response and the baseline(0)
    test = permute(pow_group, [2 3 1 4]);
    baseline = zeros(size(test));
    
    % do the permutation with FDR-corrected t-test
    pvals = std_stat({ test baseline }', 'method', 'permutation', 'condstats', 'on', 'correctm', 'fdr');  
    tmpersp = squeeze(mean(test,4)); % average ERSP for all subjects
    tmpersp(pvals{1} > 0.05) = 0; % zero out non-significant values
    
    imagesc(dstimes3, freqs, ersp_group_mean(:,dstimes2));
    ersp_fig_setup();
    hold on;
    % mark the significant regions
    contour(dstimes3, 2:50, squeeze(mean(tmpersp(:,dstimes2,:),3)), [1,1], '--k', 'LineWidth', 3)
    hold off;
    title(strcat("Posterior Channels on All Trials on ", type{t}, " (FDR-corrected)"), 'Fontsize', 25);
    savefig(strcat("allsubjects_ersp_", type{t},"_heatmap_fdr_corrected"));
    saveas(gcf, strcat("allsubjects_ersp_", type{t},"_heatmap_fdr_corrected.png"));
    saveas(gcf, strcat("allsubjects_ersp_", type{t},"_heatmap_fdr_corrected.svg"));
end