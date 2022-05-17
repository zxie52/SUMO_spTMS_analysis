%% Step 1: filter trials into four groups: left_AMI, left_UMI, right_AMI, right_UMI
p = struct2table(EEG.event); % change the EEG.event to a cell array and then categorize into different bins

% we will take 7 bins for the IEM
nbins = 7;
% calculate the edge between each bin
binedges = linspace(1,181,nbins+1);
% calculate the center of each bin
bincent = round(mean([binedges(1:end-1);binedges(2:end)]));

% filter only correct trials in the probe1
q = p((p.response == 1), :);

%Filter trials into left and right group(tms1)
[left_bin_AMI, left_bin_UMI, right_bin_AMI, right_bin_UMI] = filter_bins_for_iem_tms(q, binedges, bincent);

%% Step 2: Preprocess data befor calculating Mahalanobis Distance
groups = {left_bin_AMI, left_bin_UMI, right_bin_AMI, right_bin_UMI};

% have the EEG data and stimlabels for the iem function
% the output h: epoch# * stimlabels
h = nonzeros(groups{3});
[~, colIdcs] = find(groups{3} ~= 0);
h(:,2) = bincent(colIdcs);
% sort the h to have the ordered stimlabels
h = sortrows(h);

stimlabels = h(:,2);
super_charge = EEG.data(:,:,h(:,1));

%% Calculate the Mahalanobis Distance
data = single(permute(super_charge, [3 1 2]));
theta = deg2rad(stimlabels * 2);
n_folds = 5;

% Ester
[distance_cos,distances] = mahal_func_theta_kfold(data,theta,n_folds);

% Wolff, 2017
%[cos_amp, d_tune] = mahalTune_func(data,theta,deg2rad(bincent*2),2*pi/7);

%% Plotting
figure;
hold on;
ylim([-0.02, 0.02]);
scatter(EEG.times, squeeze(mean(distance_cos,1)), 'filled');
plot(EEG.times, squeeze(mean(distance_cos,1)));

figure;
hold on;
ylim([-0.02, 0.02]);
scatter(downsample(EEG.times, 4), downsample(squeeze(mean(distance_cos,1)), 4));

figure;
hold on;
ylim([-0.02, 0.02]);
scatter(downsample(EEG.times, 4), smoothdata(downsample(squeeze(mean(distance_cos,1)), 4), 'gaussian', 12));

title("Cosine Weighted Average of Mahalanobis Distances across timepoints on Probe1(right AMI)(testorient)", 'Fontsize', 20);

% % figure;
% % imagesc(EEG.times, bincent, squeeze(mean(distances, 2)));
% % set(gca, 'YDir', 'normal');
% 
% A = squeeze(mean(distances,2));
% 
% figure;
% for i = 1:size(A,1)
%     hold on;
%     scatter(EEG.times, A(i,:), [], rand(1,3), 'filled');
% end
% 
% figure;
% imagesc(EEG.times, bincent, A);
% set(gca, 'YDir', 'normal');
% 
% % orioffset = [1:7];
% % contourf(EEG.times, orioffset, squeeze(mean(distances, 1)),30,'linec','none');
% for i = 1 : size(A,1)
%     B(i,:) = smoothdata(A(i,:), 'gaussian', 12);
% end
% figure;
% for i = 1:size(B,1)
%     hold on;
%     scatter(EEG.times, B(i,:), [], rand(1,3), 'filled');
% end
% figure;
% for i = 1:size(B,1)
%     hold on;
%     plot(EEG.times, B(i,:), 'Color', rand(1,3));
% end
% 
% 
% figure;
% imagesc(EEG.times, bincent, B);
% set(gca, 'YDir', 'normal');
