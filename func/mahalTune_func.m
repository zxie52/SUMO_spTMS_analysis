function  [cos_amp, d_tune] = mahalTune_func(data,theta,angspace,bin_width)

% computes the mahalanobis distances between the data of single test-trial of a particular orientation and the rest of the data,
% averaged into specific orientation bins relative to orientation of the test-trial.

%% input
% data is trial by channel by time

% theta is angles in radians (-pi to pi)

% angspace (in radians) is angular space of the tuning curve, where each
% number is the center of the angle bins relative to the orienation of the test-trial.

% bin_width (in radians) is the range of orientations that are included in
% each angle bin

%% output
% d_tune is the trial-wise distance tuning curves
% cos_amp is the cosine amplitude of the tuning curve, a good summary value
% for decoding accuracy (see Sprague, T. C., Ester, E. F. & Serences, J. T. Restoring Latent Visual Working Memory Representations in Human Cortex. Neuron 91, 694ñ707 (2016))
% This code was written by Michael Wolff and Mark Stokes (see Wolff et al., 2017, Nature Neuroscience)
%%
if size(angspace,2)==length(angspace)
    angspace=angspace';
end
d_tune=nan(size(data,1),length(angspace),size(data,3));
cos_amp=nan(size(data,1),size(data,3));
trl_ind=1:size(data,1);
reverseStr=''; %this is simply to present the percentage completed
data = single(data);
for trl=1:size(data,1)
    trn_dat = data(setdiff(trl_ind,trl),:,:);
    trn_angle = theta(setdiff(trl_ind,trl));
    m=nan(length(angspace),size(trn_dat,2),size(trn_dat,3));
    for b=1:length(angspace)
        % average the training data into orientation bins relative to the test-trial's orientation
        m(b,:,:)=mean(trn_dat(abs(angle(exp(1i*trn_angle)./exp(1i*(theta(trl)-angspace(b)))))<bin_width,:,:),1);
    end
    msg=sprintf('%d percent\n',round((trl/size(data,1))*100));
    fprintf([reverseStr,msg]);
    reverseStr=repmat(sprintf('\b'),1,length(msg));
    for ti=1:size(data,3)
        if ~isnan(trn_dat(:,:,ti))
            % the covariance matrix is computed for each time-point and excluding the test-trial
            sigma = covdiag(trn_dat(:,:,ti)); % here covariance is estimated using shrinkage (see below from Ledoit & Wolf, 2004).
            % calculates the distances between the trial and all angle bins
            d_tune(trl,:,ti) = pdist2(squeeze(m(:,:,ti)), squeeze(single(data(trl,:,ti))),'mahalanobis',sigma);
            % convolve cosine of angspace with the tuning curve.
            % Since the decoding tuning curve should resemble a reversed cosine (higher distance=higher value), the value is reversed for ease of interpretation, so that high=high decoding
            cos_amp(trl,ti)=-(mean(cos(angspace).*squeeze(d_tune(trl,:,ti))'));
        end
    end
end
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