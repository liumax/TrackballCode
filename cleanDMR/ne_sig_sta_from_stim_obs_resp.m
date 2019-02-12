function [sta_sig, ptd, siglevel] = ne_sig_sta_from_stim_obs_resp(sta, locator, stim, nreps, nlags, pval)

% 
%    [sta_sig, siglevel, rand_dist] = ...
%       ne_sig_sta_from_stim_obs_resp(sta, locator, stimulus, nreps, pval)
%
%    computes a randomization test on STA to determine the significant pixels
%    in STA. The test is performed by shifting the spike train, in LOCATOR, by
%    an amount equal to 
% 
%       shiftsize = round( length(locator)/(nreps+1) );
% 
%    and then calculating a random STA. NREPS determines the shift size of the
%    spike train, and also the size of random distribution. PVAL is the
%    significance level of the test.
% 
%    STA_SIG is the significant part of STA at the level PVAL. SIGLEVEL is 
%    the level that was used to determine signifcant parts of STA, and
%    RAND_DIST is the distribution of randomized STA values. It is a vector 
%    of length = NREPS*length(STA(:))
% 
%    Inputs
%    ----------------------------------------------------------
%    sta : spike-triggered average previously estimated (each row
%          represents one vectorized STA)
%    locator : spike train
%    stim : matrix of the stimulus. Each row is one trial.
%    nreps : number of iterations
%    nlags : number of bins before each spike to use for STRF calculation
%    pval : significance level
% 
%    Outputs
%    ----------------------------------------------------------
%    sta_sig : significant parts of sta
%    ptd : peak-to-trough difference of shuffled STAs
%    siglevel : value used to threshold sta
% 
%    Craig Atencio
%    2/1/13
% 
%    Updated 5/29/18 by JS to use quick_calc_sta.m
%    Updated 6/6/18 by JS to get peak-trough difference (ptd) for shuffled
%    data.

if ~exist('pval', 'var')
    pval = 95;
end


shiftsize = round( length(locator)/(nreps+1) );

stalength = size(sta,2);
sta_rand_mat = zeros(size(sta,1),nreps*stalength);
ptd = zeros(size(locator,1), nreps);
     
fprintf('\nCalculating significant STAs...\n')

for i = 1:nreps
    
    fprintf('\nIteration %d of %d', i, nreps)
    shift = i * shiftsize;
    loc_rand = circshift(locator, [0 shift]);
    temp = quick_calc_sta(stim, loc_rand, nlags);
    ptd(:,i) = max(temp, [], 2) - min(temp, [], 2);
    sta_rand_mat(:,(i-1)*stalength+1:i*stalength) = temp;
    
end

lower = (100 - pval)/2;
upper = 100 - lower;

siglevel = prctile(sta_rand_mat,[lower upper],2);

sta_sig = zeros(size(sta));

for j = 1:size(sta,1)
    sta_sig(j,  sta(j,:) < siglevel (j,1) | sta(j,:) > siglevel(j,2)) = ...
        sta(j, sta(j,:) < siglevel (j,1) | sta(j,:) > siglevel(j,2));
end

fprintf('\n')
return;






