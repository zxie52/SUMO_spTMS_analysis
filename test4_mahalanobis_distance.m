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

%Filter trials into left and right group
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
