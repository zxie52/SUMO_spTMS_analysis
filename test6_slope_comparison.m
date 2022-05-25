% do the voltage cue1, correlate the left AMI channel response with right UMI channel response from IEM reconstruction
% do the voltage cue1, correlate the right AMI channel response with left UMI channel response from IEM reconstruction

fpath13 = 'E:\SUMO_further_data_pack_zx\N2pc_IEM\new_results\IEM_results_voltage_cue1';
subject = {'SUMO_0102', 'SUMO_0104', 'SUMO_0105', 'SUMO_0106',  ...
           'SUMO_0108', 'SUMO_0111', 'SUMO_0114', 'SUMO_0120',...
           'SUMO_3001', 'SUMO_3017', 'SUMO_3015'};

group_name = {'left AMI', 'left UMI', 'right AMI', 'right UMI'};

for i = 1 : 2 % two comparisons in cue1, voltage
    for l = 1 : length(subject)
        switch i
            case 1 % for comparing left AMI and right UMI
                load(strcat('IEM_Exp1_contraposterior', 'SUMO_', subject{l}, '_epoch1_left AMI.mat'));
                load(strcat('IEM_Exp1_contraposterior', 'SUMO_', subject{l}, '_epoch1_right UMI.mat'));
            case 2 % for comparing right AMI and left UMI
                load(strcat('IEM_Exp1_contraposterior', 'SUMO_', subject{l}, '_epoch1_left UMI.mat'));
                load(strcat('IEM_Exp1_contraposterior', 'SUMO_', subject{l}, '_epoch1_right AMI.mat'));
            otherwise
                error("there should be only two types of comparisons")
        end
        % for IEM channel response 
        r = squeeze(chanrespsL); %MLW: 6 orientation bins over 188 time points (average over 1D = subjects)
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
            slope1(t) = m(1); 
        end

        % calculate the mean across all subjects
        slope3 = smoothdata(slope1, 'gaussian', 16);

        figure; plot(dstimes, slope3);
        xlim([dstimes(1), dstimes(end)]);
        ylim([-1, 1]);
        xlabel('time from stimulus onset');
        ylabel('slope');
        xlim([dstimes(1) dstimes(end)]);

        % for permuation
        %%
        r = squeeze(chanrespsR); %MLW: 6 orientation bins over 188 time points (average over 1D = subjects)
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
            slope2(t) = m(1); 
        end

        slope4 = smoothdata(slope2, 'gaussian', 16);
        % calculate the mean across all subjects
    %     hold on;
    %     plot(dstimes, slope4);
    %     hold off;
        switch i
            case 1
                a(1,l) = corr(slope1', slope2');
            case 2
                a(2,l) = corr(slope1', slope2');
        end
    end
end


% % do the t-test, right tailed, .05
% [~, p] = ttest(slope1,slope2,.05,'right');
% 
% % plot the uncorrected p-value
% figure;
% hold on;
% xlim([dstimes(1), dstimes(end)]);
% ylim([0 1]);
% yline(.05, 'lineWidth', 4);
% scatter(dstimes, p, 'filled');
% 
% % Use cluster based permutation to corrrect the p-value
% threshSize = clustthresh1D(slope1, slope2, 10000, 'right');
% S=regionprops(p<.05,'PIxelIdxList','Area');
% idx = {};
% j = 1;
% for i = 1: length(S)
%     if S(i).Area >= threshSize
%         idx{j} = S(i).PixelIdxList';
%         j = j + 1;
%     end
% end
% % plot the clusetered corrected p-value
% scatter(dstimes(cell2mat(idx)), p(cell2mat(idx)), 'filled', 'r');
% hold off;
% 
% % do the BH-FDR for correcting the p-value
% figure;
% fdr2 = mafdr(p, 'BHFDR', true);
% hold on;
% xlim([dstimes(1), dstimes(end)]);
% ylim([0 1]);
% yline(.05, 'lineWidth', 4);
% scatter(dstimes, fdr2, 'filled');
% hold off;
