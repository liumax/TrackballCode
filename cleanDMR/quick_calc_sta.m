function [sta, stabigmat, spkcountvec] = quick_calc_sta(stimulus, locator, nlags, chunks)

% Quick calculate STA (using matrix multiplication) based on the system's
% memory.

% Inputs:
%   stimulus: (nf * nt) stimulus matrix.
%   locator: (m * nt) spike train matrix, where m = number of neurons.
%   nlags: number of time bins to look at before each spike.
%   chunks: number of chunks to calculate STA with

% Outputs:
%   sta: sta for each neuron (1 per row)
%   stabigmat: sta for each chunk of stimulus/spktrain segment (neurons x
%   sta size x number of chunks)
%   spkcountvec: spike count for each chunk of spike train

% Written by JS, 12/7/16.
% Updated by JS, 6/18/18, included 'chunks' option to manually input 
% number of chunks. Included stabigmat as an output.

narginchk(3,4)

assert(size(stimulus,2) == size(locator,2));

nf = size(stimulus,1);

if ~exist('chunks','var')
    % calculate max array based on memory
    user = memory;
    MaxElements = user.MaxPossibleArrayBytes / 8;
    limit = MaxElements / 5;

    % dimensions of stim matrix
    size1 = size(stimulus,2) - nlags + 1;
    size2 = nlags * nf;
    numelements = size1*size2;

    % determine number of chunks
    chunks = ceil(numelements / limit);
    
end

chunksize = floor(size(stimulus,2) / chunks);

% get chunk indices
startidx = [1 chunksize.*(1:chunks-1)-nlags + 2]; % include previous nlag bins for calculation of STA
endidx = [chunksize.*(1:chunks-1) size(stimulus,2)];

stabigmat = zeros(size(locator,1), nf * nlags, chunks);
spkcountvec = zeros(size(locator, 1), chunks);

fprintf('\nNumber of chunks for STA calculation: %d', chunks)

for i = 1:chunks
    fprintf('\nCalculating chunk %d of %d', i, chunks)  
    [stim, resp] = ne_create_stim_trial_from_stim_matrix(stimulus(:,startidx(i):endidx(i)),...
        locator(:,startidx(i):endidx(i)), nlags);
    spkcountvec(:,i) = sum(locator(:,startidx(i):endidx(i)), 2);
    stabigmat(:,:,i) = resp*stim;
end

sta = sum(stabigmat, 3);

fprintf('\n')

