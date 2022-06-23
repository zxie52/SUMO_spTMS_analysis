%% plot the distribution of the difference
figure;
hold on;
C = CircHist(q.orientdiff);
hold off;

diary("log2.txt");
for i = 2:20
    cid = circ_clust(q.orientdiff, i);
    fprintf(strcat("the number of cluster:", num2str(i), " and the average is:", num2str(mean(cid)), "\n"));
    for j = 1:i
%         fprintf(strcat(num2str(j),":", num2str(sum(cid == j)/length(cid)), "\n"));
        fprintf(strcat(num2str(j),":", num2str(sum(cid == j)), "\n"));
    end
end
diary('off');

%% have 3 clusters in total, by using the circular clustering

p = struct2table(EEG.event); % change the EEG.event to a cell array and then categorize into different bins
q = p((p.response == 1), :);
q.orientdiff = q.leftori - q.rightori;

clust_size = 3;
cid = circ_clust(q.orientdiff, 3, 1);
super_charge = EEG.data(:,:,q.epoch);
data = single(permute(super_charge, [3 1 2]));

% change the second column as cluster's average degree
for i = 1:length(cid)
    switch cid(i)
        case 1
            cid(i,2) = deg2rad(mean(q.orientdiff(cid==1)));
        case 2
            cid(i,2) = deg2rad(mean(q.orientdiff(cid==2)));
        case 3
            cid(i,2) = deg2rad(mean(q.orientdiff(cid==3)));
        otherwise
            error(strcat("we only pick ", num2str(clust_size), " clusters"));
    end
end
% theta = cid;
theta = cid(:,2);

j = 1;
for i = 1 : length(tmp)
    if tmp(i,2) ~= 0
        tmp2(j,1) = tmp(i,1);
        tmp2(j,2) = tmp(i,2);
        j = j + 1;
    end
end
tmp2(:,3) = q.orientdiff(tmp2(:,1));

ROI = {'CP2', 'CP4', 'CP6', 'TP8', 'TP10', ...
       'P2', 'P4', 'P6', 'P8', ...
       'PO4', 'PO8', 'O2',...
       'CP1', 'CP3', 'CP5', 'TP7', 'TP9', ...
       'P1', 'P3', 'P5', 'P7', ...
       'PO3', 'PO7', 'O1'};
impchan = find(ismember({EEG.chanlocs.labels}, ROI)); %channels in R hem

super_charge = EEG.data(impchan,:,tmp2(:,1));
theta = tmp2(:,2);
%theta = deg2rad(tmp2(:,4));
data = single(permute(super_charge, [3 1 2]));
n_folds = 10;

% Ester
[distance_cos,distances] = mahal_func_theta_kfold(data,theta,n_folds);
% distance_cos_group(l,:,i) = squeeze(mean(distance_cos,1));

figure; hold on;plot(EEG.times, mean(distance_cos,1));
figure; hold on;scatter(tmp2(:,3), tmp2(:,1), 20, tmp2(:,2), 'filled')

%% IEM
stimlabels = tmp2(:,4)';

[chanresp, ~, dstimes] = iemori(super_charge, stimlabels,4,EEG.times);
nperm = 100;
parfor p = 1:nperm
    disp(p)
    [tmp, ~, dstimes] = iemori(super_charge,stimlabels(randperm(length(stimlabels))),4,EEG.times); 
    %the 4,EEG.times is the vector of times with a downsampling factor of 4
    chanresp_perm(:,:,p) = mean(tmp,3);
end

%% Using linear K-means clustering
figure; hold on; scatter(1:length(q.orientdiff),q.orientdiff, 'filled');

clust = zeros(size(q.orientdiff,1), 20);
for K = 1 : 20
    clust(:,K) = kmeans(q.orientdiff, K);
end

eva = evalclusters(q.orientdiff, clust, 'CalinskiHarabasz');
figure; hold on; plot(eva);

% do clustering
classes = kmeans(q.orientdiff, 9);

figure; hold on; scatter(q.orientdiff, 100, classes, 'filled')

n_folds = 5;
% Ester
[distance_cos,distances] = mahal_func_theta_kfold(data,classes,n_folds);
figure; hold on;plot(EEG.times, mean(distance_cos,1));
