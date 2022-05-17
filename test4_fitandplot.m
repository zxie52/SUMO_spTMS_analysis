%% Yhought from linear fitting
% % Smoothing data per channel per trial, accross the 1200ms interval
for i = 1 : size(chanrespsR_allsub, 1) % # of subjects
    for j = 1 : size(chanrespsR_allsub, 2) % # of bins
        % The channel response is 1000ms (250 time points)
        % Use 12 as the time window to get 50ms per sample
        cond1(i,j,:) = smoothdata(chanrespsR_allsub(i,j,:), 'gaussian', 12);
    end
end


% have the center channel response for the left
chanRespR = squeeze(mean(cond1(:,4,:),1));

% have the permutation for the center channel
permutChanRespR = squeeze(mean(chanresp_permsR_allsub(:,4,:),1));

% plot the response one the base of permutation
figure;
hold on;
scatter(permutChanRespR, chanRespR, 25, 'b', '*');

% linear fit the channel response and permutation results
c = polyfit(permutChanRespR, chanRespR, 1);

% Evaluate fit equation using polyval
y_est = polyval(c,permutChanRespR);

% Add trend line to plot
hold on
plot(permutChanRespR,y_est,'r--','LineWidth',2)
hold off

% x = 1:10; 
% y1 = x + randn(1,10); 
% scatter(x,y1,25,'b','*') 

%% From fit and plot script
% % Smoothing data per channel per trial, accross the 1200ms interval
for i = 1 : size(chanrespsR_allsub, 1) % # of subjects
    for j = 1 : size(chanrespsR_allsub, 2) % # of bins
        % The channel response is 1000ms (250 time points)
        % Use 12 as the time window to get 50ms per sample
        cond1(i,j,:) = smoothdata(chanrespsR_allsub(i,j,:), 'gaussian', 12);
    end
end

% for IEM channel response 
for f = 1:size(chanrespsR_allsub, 1)

    r = squeeze(cond1(f,:,:)); %MLW: 6 orientation bins over 188 time points (average over 1D = subjects)
    centerind = 4; %MLW: center of 6 bins is technically 3

    %loop over timepoints and fit the slope ('selectivity') of the iem at each timepoint
    for t = 1:size(r,2) 

        %average the data over equidistant channels
        % why flip the array? -z.x
        avgdat = mean([r(:,t),flipud(r(:,t))],2);
        % why only pick the first half of the channels(from 1 to the center)? -z.x
        avgdat = avgdat(1:centerind);

        %fit a linear model and extract slope and save for each timepoint
        m = polyfit(1:centerind, avgdat', 1);
        slope1(f,t) = m(1); 
    end
end

% calculate the mean across all subjects
figure; plot(dstimes, squeeze(mean(slope1,1)));
xlim([dstimes(1), dstimes(end)]);
ylim([-.25, .25]);
xlabel('time from stimulus onset');
ylabel('slope');
xlim([dstimes(1) dstimes(end)]);

% for permuation
for f = 1:size(chanresp_permsR_allsub, 1)

    r = squeeze(chanresp_permsR_allsub(f,:,:)); %MLW: 6 orientation bins over 188 time points (average over 1D = subjects)
    centerind = 4; %MLW: center of 6 bins is technically 3

    %loop over timepoints and fit the slope ('selectivity') of the iem at each timepoint
    for t = 1:size(r,2) 

        %average the data over equidistant channels
        % why flip the array? -z.x
        avgdat = mean([r(:,t),flipud(r(:,t))],2);
        % why only pick the first half of the channels(from 1 to the center)? -z.x
        avgdat = avgdat(1:centerind);

        %fit a linear model and extract slope and save for each timepoint
        m = polyfit(1:centerind, avgdat', 1);
        slope2(f,t) = m(1); 
    end
end

% calculate the mean across all subjects
hold on;
plot(dstimes, squeeze(mean(slope2,1)));
hold off;


% do the t-test, right tailed, .05
[~, p] = ttest(slope1,slope2,.05,'right');

% plot the uncorrected p-value
figure;
hold on;
xlim([dstimes(1), dstimes(end)]);
ylim([0 1]);
yline(.05, 'lineWidth', 4);
scatter(dstimes, p, 'filled');

% Use cluster based permutation to corrrect the p-value
threshSize = clustthresh1D(slope1, slope2, 10000, 'right');
S=regionprops(p<.05,'PIxelIdxList','Area');
idx = {};
j = 1;
for i = 1: length(S)
    if S(i).Area >= threshSize
        idx{j} = S(i).PixelIdxList';
        j = j + 1;
    end
end
% plot the clusetered corrected p-value
scatter(dstimes(cell2mat(idx)), p(cell2mat(idx)), 'filled', 'r');
hold off;

% do the BH-FDR for correcting the p-value
figure;
fdr2 = mafdr(p, 'BHFDR', true);
hold on;
xlim([dstimes(1), dstimes(end)]);
ylim([0 1]);
yline(.05, 'lineWidth', 4);
scatter(dstimes, fdr2, 'filled');
hold off;
