function  [distance_cos,distances] = mahal_func_theta_kfold_b(data,theta,n_folds)

%% input

% data format is trial by channel by time

% theta is vector with angles in radians for each trial (must comprise the
% whole circle, thus, for orientation data, which is only 180 degrees, make
% sure to multiply by 2)

% n_folds is the desired number of folds cross validation
%% output
% output is trial by time, to summarize average over trials

% distance_cos is a measure of decoding accuracy, cosine weighted distances
% of pattern-difference between trials of increasinglt dissimilar
% orientations

% distances is the ordered mean-centred distances
%%
% have the pairwise circular distance between stimlabel and 0
theta=circ_dist(theta,0);

u_theta=unique(theta);

train_partitions = cvpartition(theta,'KFold',n_folds); % split data n times using Kfold
distances=nan(length(u_theta),size(data,1),size(data,3)); % prepare for output

theta_dist=circ_dist2(u_theta',theta)';

for tst=1:n_folds % run for each fold
    trn_ind = training(train_partitions,tst); % get training trial rows
    tst_ind = test(train_partitions,tst); % get test trial rows
    trn_dat = data(trn_ind,:,:); % isolate training data
    tst_dat = data(tst_ind,:,:); % isolate test data
    trn_theta =theta(trn_ind);
    
    % subsample to equalize trial numbers of training set
    m=double(nan(length(unique(u_theta)),size(data,2),size(data,3)));
    % original: n_conds = [u_theta',histc(trn_theta,u_theta)']; -> z.x
    n_conds = [u_theta,histc(trn_theta,u_theta)];
%     n_conds2 = [u_theta',histcounts(trn_theta,u_theta)];

    for c=1:length(u_theta)
        temp1=trn_dat(trn_theta==u_theta(c),:,:);
        ind=randsample(1:size(temp1,1),min(n_conds(:,2)));
        m(c,:,:)=mean(temp1(ind,:,:),1);
    end
    m = single(m);
    for ti=1:size(data,3) % decode at each time-point
        % compute pair-wise mahalabonis distance between test-trials
        % and averaged training data, using the covariance matrix
        % computed from the training data
        temp=pdist2(squeeze(m(:,:,ti)), squeeze(tst_dat(:,:,ti)),'mahalanobis',covdiag(trn_dat(:,:,ti)));
        distances(:,tst_ind,ti)=temp;
    end
end
distance_cos=squeeze(-mean(bsxfun(@times,cos(theta_dist)',distances),1)); % take cosine-weigthed mean of distances
% reorder distances so that same condition distance is in the middle
for c=1:length(u_theta)
    temp=round(circ_dist(u_theta,u_theta(c)),4);
    temp(temp==round(pi,4))=round(-pi,4);
    [~,i]=sort(temp);
    distances(:,theta==u_theta(c),:)=distances(i,theta==u_theta(c),:);
end
distances=-bsxfun(@minus,distances,mean(distances,1)); % mean-centre distances for prettier visuals
%%
    function sigma=covdiag(x)
        
        % x (t*n): t iid observations on n random variables
        % sigma (n*n): invertible covariance matrix estimator
        %
        % Shrinks towards diagonal matrix
        % as described in Ledoit and Wolf, 2004
        
        % de-mean returns
        [t,n]=size(x);
        meanx=mean(x);
        x=x-meanx(ones(t,1),:);
        
        % compute sample covariance matrix
        sample=(1/t).*(x'*x);
        
        % compute prior
        prior=diag(diag(sample));
        
        % compute shrinkage parameters
        d=1/n*norm(sample-prior,'fro')^2;
        y=x.^2;
        r2=1/n/t^2*sum(sum(y'*y))-1/n/t*sum(sum(sample.^2));
        
        % compute the estimator
        shrinkage=max(0,min(1,r2/d));
        sigma=shrinkage*prior+(1-shrinkage)*sample;
    end
end