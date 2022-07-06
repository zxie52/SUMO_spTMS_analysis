%% Load basics before the ersp preprocessing
close all
clear;

fpath8 = 'E:\SUMO_further_data_pack_zx\N2pc_IEM\new_results\eeg_before_IEM';
fpath12 = 'E:\SUMO_further_data_pack_zx\N2pc_IEM\new_results\ERSP_results_po3_po4';

subject = {'SUMO_0102', 'SUMO_0104', 'SUMO_0105', 'SUMO_0106',  ...
           'SUMO_0108', 'SUMO_0111', 'SUMO_0114', 'SUMO_0120',...
           'SUMO_3001', 'SUMO_3017', 'SUMO_3015'};
type = {'stim', 'cue1', 'tms1', 'probe1', 'cue2', 'tms2', 'probe2'};
ROI = {'PO3', 'PO4'};

%% Load subjects' eeg for preprocessing and FFT
for t = 2 %1:length(type)
    for s = 1 : length(subject)
        % load eeg data
        cd(fpath8);
        EEG = pop_loadset('filename', strcat(subject{s}, '_before_iem_', type{t}, '.set'), 'filepath', fpath8);
        
        % get the meta data for the eeg
        p = struct2table(EEG.event); % change the EEG.event to a cell array and then categorize into different bins 
 
        switch t
            case 1
                % for stimulus epoch, select all the trials
                % preallocate the space
                leftTrial = zeros(height(p), 1);
                rightTrial = zeros(height(p), 1);
                % get left trial and right trial numebers 
                leftTrial = p.epoch(p.targetlocation == 1);
                rightTrial = p.epoch(p.targetlocation == 2);
            otherwise
                % for other epochs, only pick up the correct trials
                % filter only correct trials in the probe1
                q = p((p.response == 1), :);
                % preallocate the space
                leftTrial = zeros(height(q), 1);
                rightTrial = zeros(height(q), 1);
                % get left trial and right trial numebers 
                leftTrial = q.epoch(q.targetlocation == 1);
                rightTrial = q.epoch(q.targetlocation == 2);
        end
        % only pick up po3 and po4 channels
        impchan = find(ismember({EEG.chanlocs.labels}, ROI)); 
        data = EEG.data(impchan,:,:);
        
        % FFT, from Jason's code
        [pow, ~, ~, dstimes, freqs] = dothewave(data, 1000, [2 50], 49, 4, (1), [-200 -20], EEG.times);

        % calculate the mean across the trials, then put into the 4-d tensor
        pow_left = squeeze(mean(pow(:,:,:,leftTrial),4));
        pow_left_group(:,:,:,s) = pow_left; 
        pow_right = squeeze(mean(pow(:,:,:,rightTrial),4));
        pow_right_group(:,:,:,s) = pow_right;
        
        % pow_left_group: channels(po3&po4) x nfreqs x time x subject
    end

    %% Preare for the visualization
    % calculate the average of ersp across subjects
    pow_left_group_mean = squeeze(mean(pow_left_group,4));
    pow_right_group_mean = squeeze(mean(pow_right_group,4));

    % shorten the time range from -150ms to 700ms to avoid the border effect from FFT
    dstimes2 = find(dstimes==-150):find(dstimes==700);
    dstimes3 = dstimes(dstimes2);

    %% Plotting the ersp for PO3, PO4 for left trials and right trials
    cd(fpath12);
    % plot the left trials
    % 1 for PO3 and 2 for PO4
    for i = 1:2
        imagesc(dstimes3, freqs, squeeze(pow_left_group_mean(i,:,dstimes2)));
        ersp_fig_setup();
        title(strcat(ROI{i}, " on Left Trials on ", type{t}), 'Fontsize', 25);
        savefig(strcat("allsubjects_ersp_left_", type{t},"_", ROI{i},"_heatmap"));
        saveas(gcf, strcat("allsubjects_ersp_left_", type{t},"_", ROI{i},"_heatmap.png"));
        saveas(gcf, strcat("allsubjects_ersp_left_", type{t},"_", ROI{i},"_heatmap.svg"));
    end
    
    % contra - ipsi = PO4 - PO3
    contra_ipsi = squeeze(pow_left_group_mean(2,:,dstimes2)) - squeeze(pow_left_group_mean(1,:,dstimes2));
    imagesc(dstimes3, freqs, contra_ipsi);
    ersp_fig_setup();
    title(strcat("Contra minus Ipsi on Left Trials on ", type{t}), 'Fontsize', 25);
    savefig(strcat("allsubjects_ersp_contra_ipsi_left_", type{t},"_heatmap"));
    saveas(gcf, strcat("allsubjects_ersp_contra_ipsi_left_", type{t},"_heatmap.png"));
    saveas(gcf, strcat("allsubjects_ersp_contra_ipsi_left_", type{t},"_heatmap.svg"));

    % plot the right trials
    for i = 1:2
        imagesc(dstimes3, freqs, squeeze(pow_right_group_mean(i,:,dstimes2)));
        ersp_fig_setup();
        title(strcat(ROI{i}, " on Right Trials on ", type{t}), 'Fontsize', 25);
        savefig(strcat("allsubjects_ersp_right_", type{t},"_", ROI{i},"_heatmap"));
        saveas(gcf, strcat("allsubjects_ersp_right_", type{t},"_", ROI{i},"_heatmap.png"));
        saveas(gcf, strcat("allsubjects_ersp_right_", type{t},"_", ROI{i},"_heatmap.svg"));
    end
    
    % contra - ipsi = PO3 - PO4
    contra_ipsi = squeeze(pow_right_group_mean(1,:,dstimes2)) - squeeze(pow_right_group_mean(2,:,dstimes2));
    imagesc(dstimes3, freqs, contra_ipsi);
    ersp_fig_setup();
    title(strcat("Contra minus Ipsi on Right Trials on ", type{t}), 'Fontsize', 25);
    savefig(strcat("allsubjects_ersp_contra_ipsi_right_", type{t},"_heatmap"));
    saveas(gcf, strcat("allsubjects_ersp_contra_ipsi_right_", type{t},"_heatmap.png"));
    saveas(gcf, strcat("allsubjects_ersp_contra_ipsi_right_", type{t},"_heatmap.svg"));
end