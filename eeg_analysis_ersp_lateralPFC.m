%% Load basics before the ersp preprocessing
close all
clear;

fpath8 = 'E:\SUMO_further_data_pack_zx\N2pc_IEM\new_results\eeg_before_IEM';
fpath12 = 'E:\SUMO_further_data_pack_zx\N2pc_IEM\new_results\ERSP_results_lateralPFC';

subject = {'SUMO_0102', 'SUMO_0104', 'SUMO_0105', 'SUMO_0106',  ...
           'SUMO_0108', 'SUMO_0111', 'SUMO_0114', 'SUMO_0120',...
           'SUMO_3001', 'SUMO_3017', 'SUMO_3015'};
type = {'stim', 'cue1', 'tms1', 'probe1', 'cue2', 'tms2', 'probe2'};

left_ROI = {'F5', 'F3', 'F1', 'FC5', 'FC3', 'FC1'};
right_ROI = {'F6', 'F4', 'F2', 'FC6', 'FC4', 'FC2'};
ROI = [left_ROI, right_ROI];

% choose the type between posterior lateral ersp(po3, po4) and frontal midline cluster ersp or dlPFC
condition = input("What clusetr of electrodes you want to pick, lateral posterior(lp) or frontal midline(fm) or lateral PFC(dlpfc)?\n", 's');
while ~ismember(condition, {'dlpfc'})
    warning("For frontal midline, please check the MATLAB code: eeg_analysis_ersp_frontal_midline.m")
    warning("For lateral posterior, please check the MATLAB code: eeg_analysis_ersp_po3_po4.m")
    warning("The input condition is not available, please re-type the condition: dlpfc");
    condition = input("What clusetr of electrodes you want to pick, lateral posterior(dlpfc)?\n", 's');
end

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
        % only pick up two clusters of DLPFC channels
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

    %% Plotting the ersp for DLPFC for left trials and right trials
    % (can be modulized in a few helper functions)
    cd(fpath12);
    % plot the left trials
    % 1 for leftDLPFC and 2 for rightDLPFC
    cond = ["Left", "Right"];
    for i = 1:2
        if i == 1, g = 1:6; else, g = 7:12; end
        imagesc(dstimes3, freqs, squeeze(mean(pow_left_group_mean(g,:,dstimes2),1)));
        ersp_fig_setup();
        title(strcat(cond(i)," DLPFC on Left Trials on ", type{t}), 'Fontsize', 25);
        savefig(strcat("allsubjects_ersp_left_", type{t},"_", cond(i),"_heatmap"));
        saveas(gcf, strcat("allsubjects_ersp_left_", type{t},"_", cond(i),"_heatmap.png"));
        saveas(gcf, strcat("allsubjects_ersp_left_", type{t},"_", cond(i),"_heatmap.svg"));
        
        % cluster permutation with FDR corrected t-test between channel response and the baseline(0)
        % (can be done in a function)
        % nfreqs x dstimes x nchan x subjects
        test = permute(pow_left_group,[2,3,1,4]);
        test = test(:,:,g,:);
        baseline = zeros(size(test));
        
        % do the permutation with FDR-corrected t-test
        pvals = std_stat({ test baseline }', 'method', 'permutation', 'condstats', 'on', 'correctm', 'fdr');  
        tmpersp = squeeze(mean(test,4)); % average ERSP for all subjects
        tmpersp(pvals{1} > 0.05) = 0; % zero out non-significant values
        
        % plot the ersp(FDR-corrected)
        imagesc(dstimes3, freqs, squeeze(mean(pow_left_group_mean(g,:,dstimes2),1)));
        ersp_fig_setup();
        hold on;
        % mark the significant regions
        contour(dstimes3, 2:50, squeeze(mean(tmpersp(:,dstimes2,:),3)), [1,1], '--k', 'LineWidth', 3)
        hold off;
        title(strcat(cond(i)," DLPFC on Left Trials on ", type{t}, " (FDR-corrected)"), 'Fontsize', 25);
        savefig(strcat("allsubjects_ersp_left_", type{t},"_", cond(i),"_heatmap_fdr_corrected"));
        saveas(gcf, strcat("allsubjects_ersp_left_", type{t},"_", cond(i),"_heatmap_fdr_corrected.png"));
        saveas(gcf, strcat("allsubjects_ersp_left_", type{t},"_", cond(i),"_heatmap_fdr_corrected.svg"));
    end
    
    % contra - ipsi = rightDLPFC - leftDLPFC
    contra_ipsi = mean(pow_left_group(7:12,:,dstimes2,:),1) - mean(pow_left_group(1:6,:,dstimes2,:),1);
    imagesc(dstimes3, freqs, squeeze(mean(contra_ipsi,4)));
    ersp_fig_setup();
    title(strcat("Contra minus Ipsi on Left Trials on ", type{t}), 'Fontsize', 25);
    savefig(strcat("allsubjects_ersp_contra_ipsi_left_", type{t},"_heatmap"));
    saveas(gcf, strcat("allsubjects_ersp_contra_ipsi_left_", type{t},"_heatmap.png"));
    saveas(gcf, strcat("allsubjects_ersp_contra_ipsi_left_", type{t},"_heatmap.svg"));
    
    % cluster permutation with FDR corrected t-test between channel response and the baseline(0)
    % (can be done in a function)
    % nfreqs x dstimes x nchan x subjects
    test = permute(contra_ipsi, [2 3 1 4]);
    baseline = zeros(size(test));
    
    % permutation statistics with FDR correction
    pvals = std_stat({ test baseline }', 'method', 'permutation', 'condstats', 'on', 'correctm', 'fdr');  
    tmpersp = mean(test,4); % average ERSP for all subjects
    tmpersp(pvals{1} > 0.05) = 0; % zero out non-significant values
    
    imagesc(dstimes3, freqs, squeeze(mean(contra_ipsi,4)));
    ersp_fig_setup();
    hold on;
    % mark the significant regions
    contour(dstimes3, 2:50, tmpersp, [1,1], '--k', 'LineWidth', 3)
    hold off;
    title(strcat("Contra minus Ipsi on Left Trials on ", type{t}, " (FDR-corrected)"), 'Fontsize', 25);
    savefig(strcat("allsubjects_ersp_contra_ipsi_left_", type{t},"_heatmap_fdr_corrected"));
    saveas(gcf, strcat("allsubjects_ersp_contra_ipsi_left_", type{t},"_heatmap_fdr_corrected.png"));
    saveas(gcf, strcat("allsubjects_ersp_contra_ipsi_left_", type{t},"_heatmap_fdr_corrected.svg"));
    
    % plot the right trials
    for i = 1:2
        if i == 1, g = 1:6; else, g = 7:12; end
        imagesc(dstimes3, freqs, squeeze(mean(pow_right_group_mean(g,:,dstimes2),1)));
        ersp_fig_setup();
        title(strcat(cond(i)," DLPFC on Right Trials on ", type{t}), 'Fontsize', 25);
        savefig(strcat("allsubjects_ersp_right_", type{t},"_", cond(i),"_heatmap"));
        saveas(gcf, strcat("allsubjects_ersp_right_", type{t},"_", cond(i),"_heatmap.png"));
        saveas(gcf, strcat("allsubjects_ersp_right_", type{t},"_", cond(i),"_heatmap.svg"));
        
        % cluster permutation with FDR corrected t-test between channel response and the baseline(0)
        % (can be done in a function)
        test = permute(pow_right_group,[2,3,1,4]);
        test = test(:,:,g,:);
        baseline = zeros(size(test));
        
        % do the permutation with FDR-corrected t-test
        pvals = std_stat({ test baseline }', 'method', 'permutation', 'condstats', 'on', 'correctm', 'fdr');  
        tmpersp = squeeze(mean(test,4)); % average ERSP for all subjects
        tmpersp(pvals{1} > 0.05) = 0; % zero out non-significant values
        
        % plot the ersp(FDR-corrected)
        imagesc(dstimes3, freqs, squeeze(mean(pow_right_group_mean(g,:,dstimes2),1)));
        ersp_fig_setup();
        hold on;
        % mark the significant regions
        contour(dstimes3, 2:50, squeeze(mean(tmpersp(:,dstimes2,:),3)), [1,1], '--k', 'LineWidth', 3)
        hold off;
        title(strcat(cond(i)," DLPFC on Right Trials on ", type{t}, " (FDR-corrected)"), 'Fontsize', 25);
        savefig(strcat("allsubjects_ersp_right_", type{t},"_", cond(i),"_heatmap_fdr_corrected"));
        saveas(gcf, strcat("allsubjects_ersp_right_", type{t},"_", cond(i),"_heatmap_fdr_corrected.png"));
        saveas(gcf, strcat("allsubjects_ersp_right_", type{t},"_", cond(i),"_heatmap_fdr_corrected.svg"));
    end
    
    % contra - ipsi = leftDLPFC - rightDLPFC
    contra_ipsi = mean(pow_right_group(1:6,:,dstimes2,:),1) - mean(pow_right_group_mean(7:12,:,dstimes2,:),1);
    imagesc(dstimes3, freqs, squeeze(mean(contra_ipsi,4)));
    ersp_fig_setup();
    title(strcat("Contra minus Ipsi on Right Trials on ", type{t}), 'Fontsize', 25);
    savefig(strcat("allsubjects_ersp_contra_ipsi_right_", type{t},"_heatmap"));
    saveas(gcf, strcat("allsubjects_ersp_contra_ipsi_right_", type{t},"_heatmap.png"));
    saveas(gcf, strcat("allsubjects_ersp_contra_ipsi_right_", type{t},"_heatmap.svg"));
    
    % cluster permutation with FDR corrected t-test between channel response and the baseline(0)
    % (can be done in a function)
    % nfreqs x dstimes x nchan x subjects
    test = permute(contra_ipsi, [2 3 1 4]);
    baseline = zeros(size(test));
    
    % permutation statistics with FDR correction
    pvals = std_stat({ test baseline }', 'method', 'permutation', 'condstats', 'on', 'correctm', 'fdr');  
    tmpersp = mean(test,4); % average ERSP for all subjects
    tmpersp(pvals{1} > 0.05) = 0; % zero out non-significant values
    
    imagesc(dstimes3, freqs, squeeze(mean(contra_ipsi,4)));
    ersp_fig_setup();
    hold on;
    % mark the significant regions
    contour(dstimes3, 2:50, tmpersp, [1,1], '--k', 'LineWidth', 3)
    hold off;
    title(strcat("Contra minus Ipsi on Right Trials on ", type{t}, " (FDR-corrected)"), 'Fontsize', 25);
    savefig(strcat("allsubjects_ersp_contra_ipsi_right_", type{t},"_heatmap_fdr_corrected"));
    saveas(gcf, strcat("allsubjects_ersp_contra_ipsi_right_", type{t},"_heatmap_fdr_corrected.png"));
    saveas(gcf, strcat("allsubjects_ersp_contra_ipsi_right_", type{t},"_heatmap_fdr_corrected.svg"));
end