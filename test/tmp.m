
contra_ipsi = pow_left_group(2,:,dstimes2,:) - pow_left_group(1,:,dstimes2,:);

% nfreqs x dstimes x nchan x subjects
test = permute(contra_ipsi, [2 3 1 4]);
test2 = squeeze(test(:,:,1,:));
test3 = squeeze(mean(test2, 3));
baseline = zeros(size(test));

% permutation statistics with cluster correction (note: use the latest SVN revision of EEGLAB 12)
pvals = std_stat({test baseline}', 'mode', 'fieldtrip', 'fieldtripmethod', 'montecarlo', 'condstats', 'on', 'fieldtripmcorrect', 'cluster'); 
  %
% permutation statistics with FDR correction
pvals = std_stat({ test baseline }', 'method', 'permutation', 'condstats', 'on', 'correctm', 'fdr');  
tmpersp = mean(test,4); % average ERSP for all subjects
tmpersp(pvals{1} > 0.05) = 0; % zero out non-significant values


p = imagesc(dstimes, freqs, tmpersp); 
%caxis([-20 10]);
set(gca, 'ydir', 'normal'); xlabel('Time (ms)'); ylabel('Frequencies (Hz)'); cbar; % plot ERSP

% [clusters, p_values, t_sums, permutation_distribution ] = permutest(test3', ...
%                                                                     baseline', ...
%                                                                     true, .05, 10000, false);
% 

% p_va = tmpersp;
% t = zeros(size(p_va));
% for i = 1:size(p_va,1)-1
%     for j = 1:size(p_va,2)-1
%         if p_va(i,j) - p_va(i,j+1) > 0 ||p_va(i,j) - p_va(i+1,j) > 0; t(i,j) = 1; end
%     end
% end

imagesc(dstimes, freqs, t); 
set(gca, 'ydir', 'normal'); xlabel('Time (ms)'); ylabel('Frequencies (Hz)'); cbar; % plot ERSP

[k, av] = convhull(tmpersp);









