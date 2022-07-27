function [pow, phase, comp, dstimes, freqs] = dothewave_broadband(timeseries, srate, cyclerange, dsfac, baseline, timevec)
% dsfac = 5;
% calculate the base at the first time for each subj
%  stim epoch      [pow, ~, ~, dstimes, freqs] = dothewave(EEG.data, 500, [8 13], 6, 4, [], [], EEG.times);
%  non-stim epoch  [pow, ~, ~, dstimes, freqs] = dothewave(EEG2.databs, 500, [8 13], 6, 4, [], [-800 -500], EEG.times_total)


% expects timeseries to be voltage data: CHAN x TIME x TRIALS (if input is a single channel, do not squeeze the data - leave as 1 X TIME X TRIALS)
% freqrange = range you want (for alpha: [8 13])
% nfreqs = # of frequencies you want to estimate in the range -- linearly spaced
% cyclerange = sets trade-off between temporal and frequency // how many cycles of oscillation to use to estimate power for each frequency
%       first value min. of 3 ... second value ~= 10
%       good value when just looking at alpha = 4
% dsfac = do we want to temporally down-sample the output? (makes output file smaller) dsfac = 3 means every third sample is used
% baseline = when you define this, it sets percent signal change relative to baseline i.e. [-400 -100]
% timevec = vector index that maps samples to ms ex: [EEG.times]
%
% returns power phase and complex values as CHAN x FREQ x TIME x TRIALS
%     can put '~' instead of 'phase' and 'comp' to not get that output because it'll be large
% pow will be a 4D output -- plot with dstimes on X axis, freqs on Y axis
%
% trials is optional if timeseries is 2d: then timeseries must be CHAN x TIME
% any undesired input must be set to [] 
%           jasonsamaha@gmail.com UW-Madison Psychology
%           equations based on Cohen MX (2014).

%default dsfac
if ~exist('dsfac','var')
    dsfac = 1;
end

if exist('timevec','var')
    dstimes = timevec(1:dsfac:end);
end

%length of final output
outlen = length(1:dsfac:size(timeseries,2));

%freq params
%%%this part can be changed (zx)
freqs = [5:7 8 10 12 14 18 22 26 34 42];
num_frex = length(freqs);


% other wavelet parameters
s = logspace(log10(cyclerange(1)),log10(cyclerange(end)),num_frex) ./ (2*pi*freqs);
wavtime = -2:1/srate:2;
half_wave = (length(wavtime)-1)/2;

% FFT parameters and preallocate
nWave = length(wavtime);
if length(size(timeseries)) == 3
    nData = size(timeseries,2) * size(timeseries,3);
    pow = zeros(size(timeseries,1), length(freqs),outlen,size(timeseries,3));
    phase = zeros(size(timeseries,1), length(freqs),outlen,size(timeseries,3));
    comp = zeros(size(timeseries,1), length(freqs),outlen,size(timeseries,3));
else
    nData = size(timeseries,2);
    pow = zeros(size(timeseries,1), length(freqs),outlen);
    phase = zeros(size(timeseries,1), length(freqs),outlen);
    comp = zeros(size(timeseries,1), length(freqs),outlen);
end
nConv = nWave + nData - 1;

%fprintf('\n performing wavelet convolution ') 
fprintf('\n ch ') 
parfor ch = 1:size(timeseries,1)
    fprintf('%d ',ch)

    % now compute the FFT of all trials concatenated
    if length(size(timeseries)) == 3
        chandata = reshape( timeseries(ch,:,:) ,1,[]);
        % initialize output time-frequency data
        tf = zeros(length(freqs),size(timeseries,2),size(timeseries,3));
    else
        chandata = timeseries(ch,:);
        tf = zeros(length(freqs),size(timeseries,2));
    end
    dataX = fft( chandata ,nConv );

    % now perform convolution
    % loop over frequencies
    for fi=1:length(freqs)

        % create wavelet and get its FFT
        wavelet  = exp(2*1i*pi*freqs(fi).*wavtime) .* exp(-wavtime.^2./(2*s(fi)^2));
        waveletX = fft(wavelet,nConv);
        waveletX = waveletX ./ max(waveletX); % normalize to put back in original recording units

        % now run convolution in one step
        ast = ifft(waveletX .* dataX);
        ast = ast(half_wave+1:end-half_wave);

        % and reshape back to time X trials
        if length(size(timeseries)) == 3
            ast = reshape( ast, size(timeseries,2), size(timeseries,3) );
        end

        % collect complex result over frequencies
        tf(fi,:,:) = ast;
        ast = [];
    end

    if length(size(timeseries)) == 3 && isempty(baseline)
        pow(ch,:,:,:) = abs(tf(:,1:dsfac:end,:)).^2;
        phase(ch,:,:,:) = angle(tf(:,1:dsfac:end,:));
        comp(ch,:,:,:) = tf(:,1:dsfac:end,:);
        %continue    
    end
%     elseif length(size(timeseries)) == 2 && isempty(baseline)
%         pow(ch,:,:) = abs(tf(:,1:dsfac:end)).^2;
%         phase(ch,:,:) = angle(tf(:,1:dsfac:end));
%         comp(ch,:,:) = tf(:,1:dsfac:end);
%         %continue
%     else
%         % percent signal change baseline using all data in timevec (i.e., possibly all conditions)
%         baselineidx = [];
%         baselineidx(1)=dsearchn(timevec',baseline(1));
%         baselineidx(2)=dsearchn(timevec',baseline(2));
%         if length(size(timeseries))==3
%             base = mean(mean(abs(tf(:,baselineidx(1):baselineidx(2),:)).^2,2),3);
%             pow(ch,:,:,:) = 100 * bsxfun(@rdivide, bsxfun(@minus,abs(tf(:,1:dsfac:end,:)).^2,base), base);
%             
%             % pre_stim_pow(ch,:,:,:) = 100 * bsxfun(@rdivide, bsxfun(@minus,pow, base), base);
%             phase(ch,:,:,:) = angle(tf(:,1:dsfac:end,:));
%             comp(ch,:,:,:) = tf(:,1:dsfac:end,:);
%         else
%         	base = mean(abs(tf(:,baselineidx(1):baselineidx(2))).^2,2);
%             pow(ch,:,:) = 100 * bsxfun(@rdivide, bsxfun(@minus,abs(tf(:,1:dsfac:end)).^2,base), base);
%             phase(ch,:,:) = angle(tf(:,1:dsfac:end));
%             comp(ch,:,:) = tf(:,1:dsfac:end);
%         end

    
end

fprintf(' done! \n')


