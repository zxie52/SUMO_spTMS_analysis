% % Smoothing data per channel per trial, accross the 1200ms interval
for i = 1 : size(chanrespsR_allsub, 1) % # of subjects
    for j = 1 : size(chanrespsR_allsub, 2) % # of bins
        % The channel response is 1000ms (250 time points)
        % Use 12 as the time window to get 50ms per sample
        cond1(i,j,:) = smoothdata(chanrespsR_allsub(i,j,:), 'gaussian', 12);
    end
end
% 
% % only pick up the center channel's channel response
chanRespR = squeeze(cond1(:,4,:))';
permutRespR = squeeze(chanresp_permsR_allsub(:,4,:))';
% 
% 
for i = 1 : length(chanRespR)
    [~, p(1,i)] = ttest(chanRespR(i,:),permutRespR(i,:),.05,'right');
end
% 
% plot the uncorrected p-value
figure;
hold on;
xlim([dstimes(1), dstimes(end)]);
ylim([0 1]);
yline(.05, 'lineWidth', 4);
scatter(dstimes, p, 'filled');
% 
% % plot the FDR-corrected p-value(q-value from Storey, 2002)
% figure;
% [fdr, q, pi0] = mafdr(p);
% hold on;
% xlim([dstimes(1), dstimes(end)]);
% ylim([0 1]);
% yline(.05, 'lineWidth', 4);
% scatter(dstimes(1:3:end), q, 'filled');
% 
% % plot the FDR-corrected p-value (using LSU Benjamini and Hochberg (1995))
% % finally we will pick the BHFDR
figure;
fdr2 = mafdr(p, 'BHFDR', true);
hold on;
xlim([dstimes(1), dstimes(end)]);
ylim([0 1]);
yline(.05, 'lineWidth', 4);
scatter(dstimes, fdr2, 'filled');
% 
% figure;
% fdr2 = mafdr(p, 'BHFDR', true, 'showplot', true);
% 
% plot Bonferonni-Holme Correction
% figure;
% cor_p = bonf_holm(p, .05);
% hold on;
% xlim([dstimes(1), dstimes(end)]);
% ylim([0 1]);
% yline(.05, 'lineWidth', 4);
% scatter(dstimes, cor_p, 'filled');
% 
% % same using fdr_bh fucntion to run the BHFDR
% [~, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p,.05);
% figure;
% hold on;
% xlim([dstimes(1), dstimes(end)]);
% ylim([0 1]);
% yline(.05, 'lineWidth', 4);
% scatter(dstimes(1:3:end), adj_p, 'filled');
% 
% % using the cluster permutation to correct the p-value
% threshSize = clustthresh1D(chanRespR, permutRespR, 10000, 'right');
% S=regionprops(p<.05,'PIxelIdxList','Area');
% idx = {};
% j = 1;
% for i = 1: length(S)
%     if S(i).Area >= threshSize
%         idx{j} = S(i).PixelIdxList';
%         j = j + 1;
%     end
% end
% %% Plotting the heatmap
% 
% nbins = 7;
% binedges = linspace(1,181,nbins+1);
% bincent = round(mean([binedges(1:end-1);binedges(2:end)]));
% 
% % Load Parameters
% centerind = 4;
% nchan = 7;
% 
% % Set up title font size
% titleFontSize = 20;
% axisFontSize = 14;
% textFontSize = 12;
% 
% % the process for generating the heat map
% orioffset = linspace(-90, 90, nchan); %for plotting
% figure; contourf(dstimes, orioffset, squeeze(mean(cond1, 1)),30,'linec','none'); %make a contour plot
% xlim([min(dstimes) max(dstimes)]); %can trim the time range of the plot here if you want. Right now setting to min and max does nothing.
% h = colorbar;
% caxis([-0.1 0.5]); %scale the colorbar axis to just region of interest (-0.1 to 0.6)
% 
% % Opening the figure in the full size
% set(gcf, 'Position', get(0, 'Screensize'));
% 
% ylabel(h, 'Channel Response (µV^2)')
% title(strcat("Reconstruction from Voltage in Contra-Posterior Electrodes (TMS2 ", group_name{i}, ")(q-value)"), 'Fontsize', titleFontSize)
% xlabel('Time from TMS onset(msec)','Fontsize', axisFontSize)
% ylabel('Centered Orientation Channel (°)','Fontsize', axisFontSize)
% 
% hold on;
% if ~isempty(t(fdr2<.05))
%     plot([t(fdr2<.05),t(fdr2<.05)+1,t(fdr2<.05)+2], 75, 'r.', 'MarkerSize', 25);
% end
% 
% hold on;
% if ~isempty(t(p<.05))
%     plot([t(p<.05),t(p<.05)+1,t(p<.05)+2], 80, 'g.', 'MarkerSize', 25);
% end

%% Refined code for Cluster-Permutation
% Smoothing data per channel per trial, accross the 1200ms interval
for i = 1 : size(chanrespsR_allsub, 1) % # of subjects
    for j = 1 : size(chanrespsR_allsub, 2) % # of bins
        % The channel response is 1000ms (250 time points)
        % Use 12 as the time window to get 50ms per sample
        cond1(i,j,:) = smoothdata(chanrespsR_allsub(i,j,:), 'gaussian', 12);
    end
end

% only pick up the center channel's channel response
chanRespR = squeeze(cond1(:,4,:))';
permutRespR = squeeze(chanresp_permsR_allsub(:,4,:))';

% run the t-test on each timepoint(IEM channel response vs IEM permutation)
for i = 1 : length(chanRespR)
    [~, p(1,i)] = ttest(chanRespR(i,:),permutRespR(i,:),.05,'right');
end

% using the cluster permutation to correct the p-value
threshSize = clustthresh1D(chanRespR, permutRespR, 10000, 'right');
S=regionprops(p<.05,'PIxelIdxList','Area');
idx = {};
j = 1;
for i = 1: length(S)
    if S(i).Area >= threshSize
        idx{j} = S(i).PixelIdxList';
        j = j + 1;
    end
end

%% Plotting the heatmap
% Load Parameters
nbins = 7;
binedges = linspace(1,181,nbins+1);
bincent = round(mean([binedges(1:end-1);binedges(2:end)]));

centerind = 4;
nchan = 7;

% Set up title font size
titleFontSize = 20;
axisFontSize = 14;
textFontSize = 12;

% the process for generating the heat map
orioffset = linspace(-90, 90, nchan); %for plotting
figure; contourf(dstimes, orioffset, squeeze(mean(cond1, 1)),30,'linec','none'); %make a contour plot
xlim([min(dstimes) max(dstimes)]); %can trim the time range of the plot here if you want. Right now setting to min and max does nothing.
h = colorbar;
caxis([-0.1 0.5]); %scale the colorbar axis to just region of interest (-0.1 to 0.6)

% Opening the figure in the full size
set(gcf, 'Position', get(0, 'Screensize'));

% Set the labels
% ylabel(h, 'Channel Response (µV^2)')
% title(strcat("Reconstruction from Voltage in Contra-Posterior Electrodes (TMS2 ", group_name{i}, ")(q-value)"), 'Fontsize', titleFontSize)
% xlabel('Time from TMS onset(msec)','Fontsize', axisFontSize)
% ylabel('Centered Orientation Channel (°)','Fontsize', axisFontSize)

% green for uncorrected p-values that are less than .05
% red for corrected p-values that are less than .05

% Plotting the uncorrected p-value on the heatmap
hold on;
if ~isempty(p(p<.05))
    plot(dstimes(p<.05), 80, 'g.', 'MarkerSize', 25);
end

% Plotting the clusters on the heatmap
hold on;
if ~isempty(cell2mat(idx))
    plot(dstimes(cell2mat(idx)), 75, 'r.', 'MarkerSize', 25);
else
    warning("There is no significant p-value after clusted-based permutation correction");
    warning("The largest cluster size is, %d", max([S.Area]));
end
