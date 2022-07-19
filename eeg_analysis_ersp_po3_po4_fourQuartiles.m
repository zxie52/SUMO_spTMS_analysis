%% Load basics before the ersp preprocessing
close all
clear;

fpath8 = 'E:\SUMO_further_data_pack_zx\N2pc_IEM\new_results\eeg_before_IEM';
% pathway to save out the results
fpath12 = 'E:\SUMO_further_data_pack_zx\N2pc_IEM\new_results\ERSP_results_4_quartiles';

subject = {'SUMO_0102', 'SUMO_0104', 'SUMO_0105', 'SUMO_0106',  ...
           'SUMO_0108', 'SUMO_0111', 'SUMO_0114', 'SUMO_0120',...
           'SUMO_3001', 'SUMO_3017', 'SUMO_3015'};
type = {'stim', 'cue1', 'tms1', 'probe1', 'cue2', 'tms2', 'probe2'};

ROI = {'PO3', 'PO4'};

%% Doing FFTs
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
            error("Wrong code, using 'eeg_analysis_ersp_posterior_stim_only.m' instead");
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
    % have the 4 quartiles to split the trials for both left and right trials
    lf = fix(linspace(1,length(leftTrial),5));
    rf = fix(linspace(1,length(rightTrial),5));
    
    % FFT, from Jason's code
    [pow, ~, ~, dstimes, freqs] = dothewave(data, 1000, [2 50], 49, 4, (1), [-200 -20], EEG.times);
    
    % pre-allocate the space for cells
    pow_left = {};
    pow_right = {};
    
    for i = 1 : 4
        % get the pow for specific trials
        pow_left{i} = squeeze(mean(pow(:,:,:,leftTrial(lf(i):lf(i+1))),4));
        pow_right{i} = squeeze(mean(pow(:,:,:,rightTrial(rf(i):rf(i+1))),4));
    end   
    
    % calculate the mean across the trials, then put into the 4-d tensor(save in a cell array structure)
    pow_left_group(:,s) = pow_left; 
    pow_right_group(:,s) = pow_right;
    % pow_left_group: channels(po3&po4) x nfreqs x time x subject
end

%% Data Analysis & Visualization
pow_left_group = squeeze(pow_left_group);
pow_right_group = squeeze(pow_right_group);

% Names for later plotting
tp = ["1st", "2nd", "3rd", "4th"];

% shorten the time range from -150ms to 700ms to avoid the border effect from FFT
dstimes2 = find(dstimes==-150):find(dstimes==700);
dstimes3 = dstimes(dstimes2);

cd(fpath12);

% for four quartiles in total
for i = 1 : 4
    pl = [];
    pr = [];
    for s = 1 : length(subject)
        % first dimension: 1:PO3, 2:PO4
        pl(:,:,:,s) = pow_left_group{i,s};
        pr(:,:,:,s) = pow_right_group{i,s};
    end
    
    %% calculate the contra Minus ipsi
    %% left trials: PO4 - PO3
    contra_ipsi = pl(2,:,dstimes2,:) - pl(1,:,dstimes2,:);
    
    % cluster permutation with FDR corrected t-test between channel response and the baseline(0)
    % (can be done in a function)
    % nfreqs x dstimes x nchan x subjects
    test = permute(contra_ipsi, [2 3 1 4]);
    baseline = zeros(size(test));
    
    % permutation statistics with FDR correction
    pvals = std_stat({ test baseline }', 'method', 'permutation', 'condstats', 'on', 'correctm', 'fdr');  
    tmpersp = mean(test,4); % average ERSP for all subjects
    tmpersp(pvals{1} > 0.05) = 0; % zero out non-significant values
    
    % Data visulization
    imagesc(dstimes3, freqs, squeeze(mean(contra_ipsi,4)));
    ersp_fig_setup();
    hold on;
    % mark the significant regions
    contour(dstimes3, 2:50, tmpersp, [1,1], '--k', 'LineWidth', 3)
    hold off;
    title(strcat("Contra minus Ipsi on Left Trials(", tp(i), " Quartile) on ", type{t}, " (FDR-corrected)"), 'Fontsize', 25);
    savefig(strcat("allsubjects_ersp_contra_ipsi_left_", type{t},"_heatmap_fdr_corrected_", tp(i), "_Quartile"));
    saveas(gcf, strcat("allsubjects_ersp_contra_ipsi_left_", type{t},"_heatmap_fdr_corrected_", tp(i), "_Quartile.png"));
    saveas(gcf, strcat("allsubjects_ersp_contra_ipsi_left_", type{t},"_heatmap_fdr_corrected_", tp(i), "_Quartile.svg"));
    
    %% right trials: PO3 - PO4
    clear contra_ipsi;
    contra_ipsi = pr(1,:,dstimes2,:) - pr(2,:,dstimes2,:);
    
    % cluster permutation with FDR corrected t-test between channel response and the baseline(0)
    % (can be done in a function)
    % nfreqs x dstimes x nchan x subjects
    test = permute(contra_ipsi, [2 3 1 4]);
    baseline = zeros(size(test));
    
    % permutation statistics with FDR correction
    pvals = std_stat({ test baseline }', 'method', 'permutation', 'condstats', 'on', 'correctm', 'fdr');  
    tmpersp = mean(test,4); % average ERSP for all subjects
    tmpersp(pvals{1} > 0.05) = 0; % zero out non-significant values
    
    % Data visulization
    imagesc(dstimes3, freqs, squeeze(mean(contra_ipsi,4)));
    ersp_fig_setup();
    hold on;
    % mark the significant regions
    contour(dstimes3, 2:50, tmpersp, [1,1], '--k', 'LineWidth', 3)
    hold off;
    title(strcat("Contra minus Ipsi on Right Trials(", tp(i), " Quartile) on ", type{t}, " (FDR-corrected)"), 'Fontsize', 25);
    savefig(strcat("allsubjects_ersp_contra_ipsi_right_", type{t},"_heatmap_fdr_corrected_", tp(i), "_Quartile"));
    saveas(gcf, strcat("allsubjects_ersp_contra_ipsi_right_", type{t},"_heatmap_fdr_corrected_", tp(i), "_Quartile.png"));
    saveas(gcf, strcat("allsubjects_ersp_contra_ipsi_right_", type{t},"_heatmap_fdr_corrected_", tp(i), "_Quartile.svg"));    
end
