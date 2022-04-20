
function [chanresp, weights, dstimes] = iemori(dat, stimlabels, dsfac, timevec)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   REQUIRED INPUTS
%
%   dat:        chans X timepoints X trials of features (e.g., voltage, power)
%               or, chans X trials in the case of non-timeresolved inputs
%   
%   stimlabels: stimulus feature on each trial (i.e., orientation *in degrees*)
%               center of the bin
%
%
%   OPTIONAL INPUTS
%
%   dsfac:      optional downsampling factor to speed up (e.g., dsfac = 3
%               cuts data in thirds along the temporal dimension
%   
%   times:      vector of time (e.g., in ms or s). Used if dsfac is given  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% to do: 
% - get weights
% - check that it works for different dimensions
% - think about implimenting train test on different data


%check dsfac
if ~exist('dsfac','var')
    dsfac = 1;
elseif exist('timevec','var')
    dstimes = timevec(1:dsfac:end);
    dat = dat(:,1:dsfac:end,:);
end


%check whether dat has temporal dimension and get dimension with trials
if size(dat) == 2
    nsamp = 1;
    trldim = 2;
else
    nsamp = size(dat,2);
    trldim = 3;
end

%assume number of channels equals number of unique stimulus orientations
n_ori_chans = length(unique(stimlabels));

%function handle to make sinusiod basis function raised to the nth power
make_basis_function = @(xx,mu) (cosd(xx-mu)).^(n_ori_chans-mod(n_ori_chans,2));

%create basis set - assuming stimuli presented across all ranges of orientations
xx = linspace(1,180,180);
basis_set = nan(180,n_ori_chans);
chan_center = unique(stimlabels);%linspace(180/n_ori_chans,180,n_ori_chans);

for cc = 1:n_ori_chans
    basis_set(:,cc) = make_basis_function(xx,chan_center(cc));
end

%recode stimlabels to be between 1 and 180 - e.g., assumes 270 deg == 90 deg
stimlabels = mod(stimlabels,180);
stimlabels(stimlabels==0) = 180;

%create C1 matrix of expected channel responses for each trial
stim_mask = zeros(length(stimlabels),length(xx));

for tt = 1:size(stim_mask,1)  % loop over trials
    stim_mask(tt,stimlabels(tt))=1; 
end

%make C1
C1 = stim_mask*basis_set;

%set up leave-one-trial-out cross validation structure
tpart = cvpartition(stimlabels,'Leaveout');
predresp = zeros(tpart.NumTestSets,n_ori_chans,nsamp);
weights = zeros(n_ori_chans,size(dat,1),tpart.NumTestSets,nsamp);

%start loop over timepoints
for samp = 1:nsamp   
    %get current sample
    sampdat = squeeze(dat(:,samp,:));
    sampdat = sampdat';
    
    %preallocate
    tmpweights = zeros(n_ori_chans, size(dat,1), size(dat,trldim));

    %start loop of folds
    for w = 1:tpart.NumTestSets            
        %estimate training weights for this fold
        tWeights = C1(tpart.training(w),:)\sampdat(tpart.training(w),:);  
        
        %get predicted response
        predresp(tpart.test(w),:,samp) = (tWeights'\sampdat(tpart.test(w),:)')';
        
        %save weights from different runs
        tmpweights(:,:,w) = tWeights;
    end
    
    weights(:,:,:,samp) = tmpweights;
end

%align to common center of 90 degs
cent_ori = chan_center(round(length(chan_center)/2));
cent_ori_idx = find(chan_center==cent_ori);

for tr = 1:size(predresp,1)
    thisidx = find(chan_center == stimlabels(tr));
    chanresp(:,:,tr) = circshift(predresp(tr,:,:),cent_ori_idx - thisidx,2);    
    C1_cent(tr,:) = circshift(C1(tr,:),cent_ori_idx - thisidx,2);
end
