%%%%%
%%%%% This file is the test for the dothewave function:
%%%%% one is using the original dothewave fucntion
%%%%% The other is using the manual subtraction from the pre-stimulus baseline


%% prac1 for original way to calculate the alpha power
prac1 = EEG.data(1:3,100:200,1:5);
[pow, ~, ~, dstimes, freqs] = dothewave(prac1, 1000, [8 13], 6, 4, 1, [1 20], 1:101);
result1 = squeeze(mean(pow,2));

%% prac2 for manual way to calculate the alpha(baseline)
prac2 = EEG.data(1:3, 100:119, 1:5);
[pow_base, ~, ~, dstimes, freqs] = dothewave(prac2, 1000, [8 13], 6, 4, 1, [], 1:20);
new_pow_base = mean(mean(pow_base,2), 3);

prac3 = EEG.data(1:3,100:200,1:5);
[pow2, ~, ~, dstimes, freqs] = dothewave(prac3, 1000, [8 13], 6, 4, 1, [], 1:101);
tmp = 100 * bsxfun(@rdivide, bsxfun(@minus, pow2, new_pow_base), new_pow_base);

result2 = squeeze(mean(tmp,2));

%% Without comments
[pow, ~, ~, dstimes, freqs] = dothewave(EEG.data, 1000, [8 13], 6, 4, 1, [-200 -20], EEG.times);
result1 = squeeze(mean(pow,2));
p1 = result1(1,1,1);

[pow_base, ~, ~, dstimes, freqs] = dothewave(EEG.data, 1000, [8 13], 6, 4, 1, [], EEG.times);
A = pow_base(:,:,find(EEG.times == -200):find(EEG.times == -20),:);
new_pow_base = mean(mean(A,2), 3);
[pow2, ~, ~, dstimes, freqs] = dothewave(EEG.data, 1000, [8 13], 6, 4, 1, [], EEG.times);
tmp = 100 * bsxfun(@rdivide, bsxfun(@minus, pow2, new_pow_base), new_pow_base);
result2 = squeeze(mean(tmp,2));
p2 = result2(1,1,1);

A = EEG.data(:, find(EEG.times == -200):find(EEG.times == -20), :);
B = EEG.times(1, find(EEG.times == -200):find(EEG.times == -20));
[pow_base2, ~, ~, dstimes, freqs] = dothewave(A, 1000, [8 13], 6, 4, 1, [], B);
new_pow_base2 = mean(mean(pow_base2,2), 3);
[pow3, ~, ~, dstimes, freqs] = dothewave(EEG.data, 1000, [8 13], 6, 4, 1, [], EEG.times);
tmp = 100 * bsxfun(@rdivide, bsxfun(@minus, pow2, new_pow_base2), new_pow_base2);
result3 = squeeze(mean(tmp,2));
p3 = result3(1,1,1);


[pow_base, ~, ~, dstimes, freqs] = dothewave(EEG.data, 1000, [8 13], 6, 4, 1, [], EEG.times);
% average across timepoints
new_pow_base = mean(pow_base(:,:,find(EEG.times == -200):find(EEG.times == -20),:),3);
[pow4, ~, ~, dstimes, freqs] = dothewave(EEG.data, 1000, [8 13], 6, 4, 1, [], EEG.times);
tmp = 100 * bsxfun(@rdivide, bsxfun(@minus, pow4, new_pow_base), new_pow_base);
result4 = squeeze(mean(tmp,2));


[pow_base, ~, ~, dstimes, freqs] = dothewave(EEG.data, 1000, [8 13], 6, 4, 1, [], EEG.times);
% average across timepoints and acros trials
new_pow_base = mean(mean(pow_base(:,:,find(EEG.times == -200):find(EEG.times == -20),:),3), 4);
[pow5, ~, ~, dstimes, freqs] = dothewave(EEG.data, 1000, [8 13], 6, 4, 1, [], EEG.times);

tmp = 100 * bsxfun(@rdivide, bsxfun(@minus, pow5, new_pow_base), new_pow_base);
result5 = squeeze(mean(tmp,2));

randChan = randi(length(EEG.allchan));
randTrial = randi(length(EEG.epoch));
figure; plot(squeeze(result1(randChan,:,randTrial))); hold on; plot(squeeze(result5(randChan,:,randTrial)));
legend("Using dothewave function", "Manually subtract the baseline", 'Fontsize', 20);
