function [stim, varargout] = ne_create_stim_trial_from_stim_matrix(stimulus, locator, nlags)
% ca_create_stim_trial_from_stim_matrix Stimulus in observation form.
% 
%     [stim, resp] = ca_create_stim_trial_from_stim_matrix(stimulus, locator, nlags)
% 
%     Takes the stimulus matrix, in 'stimulus', and creates a new stimulus, 
%     where each row represents one observation, and the number of columns
%     equals the dimensionality of the stimulus, which is NF X nlags, or the 
%     number of frequencies that will be used to estimate the receptive 
%     field, and nlags is the number of time bins to include.
% 
%     stimulus : NF x Ntrials matrix, NF = number of frequencies, and Ntrials = 
%     total number of time bins in the stimulus.
% 
%     locator : spike train vector or matrix. Each row of locator contains
%     the spike count of a neuron in the time bins.
% 
%     nlags : number of time bins to use for the memory of the receptive
%     field. For AI cells, and 5 ms bins, nlags = 20 is sufficient.
% 
%     stim : stimulus observation matrix. Each row is an observation. Columns
%     represent the stimulus dimensions.
% 
%     resp : spike train matrix that has been adjusted to correspond to stim.

% Updated 7/29/16 by JS.

if ( nargin == 2 )
    nlags = 20; % # of time bins in STA
end

Nsamples = size(stimulus,2);

% Spike train length must equal stimulus length

if ~isempty(locator) && ( size(locator,2) ~= size(stimulus,2) )
   error('Number of trials in locator and stimulus do not match.');
end

Ndim = size(stimulus,1); % # dimensions per time lag
nf = Ndim; % New # of frequencies

if ( nlags > 1 )

   Nsamples = Nsamples - (nlags-1);
   Ndimtotal = nf*nlags; % total stimulus dimensions = #freq's by #time lags

   stim = zeros(Nsamples,Ndimtotal);

   for i = 1:Nsamples
      chunk = stimulus(:,i:i+nlags-1);
      stim(i,:) = chunk(:)'; % make a row vector and assign to stim
   end

else % only happens for visual stimuli
    stim = stimulus;
end


% Open responses and cut according to # of time lags

if ~isempty(locator)
    varargout{1} = locator(:, nlags:length(locator(1,:)));
end


return;


