function params = strf_parameters(strf, trigger)
%strf_parameters - get strf parameters from significant portions
%   of the strf and the rtf (ripple transfer function.
%
% Calculates temporal, spectral modulation transfer
% functions as well as the separability index for
% all strfs.
%
% params = strf_parameters(strf, trigger)
%
% You can use the function plot_strf_parameters.m to
% view the results of the analysis. 
%
%   caa 10/20/02


% parameters to obtain:
% separability index, i.e. the singular and eigen- values
% temporal correlation index
% spectral correlation index
% asymmetry index
% phase-locking index
% best tm
% best sm
% mtfs: tm and sm
% strf energy - use function strf_energy.m
% firing rate
% number of spikes
% duration
% 
% and later, feature selectivity index


params = struct(...
'exp',            [], ...
'site',           [], ...
'chan',           [], ...
'model',          [], ...
'depth',          [], ...
'position',       [], ...
'stim',           [], ...
'atten',          [], ...
'sm',             [], ...
'tm',             [], ...
'mdb',            [], ...
'spl',            [], ...
'n0',             [], ...
'w0',             [], ...
'percent_energy', [], ...
'rfenergy',       [], ...
'tmf',            [], ...
'xmf',            [], ...
'rtf',            [], ...
'singvals',       [], ...
'eigvals',        [], ...
'tci',            [], ...
'sci',            [], ...
'pli',            []);


percent_energy = [80 85 90 95 97.5 100];
energy = zeros(size(percent_energy));

for i = 1:length(strf)
i
   % Define some initial parameters
   p = 0.001;
   n0 = strf(i).n0contra;
   w0 = strf(i).w0contra;
   sm = strf(i).sm;
   tm = strf(i).tm;
   mdb = strf(i).mdb;
   stim = strf(i).stim;
   fs = strf(i).fs;
   t = strf(i).taxis;
   x = log2(strf(i).faxis ./ min(strf(i).faxis));
   dur = (trigger(end)-trigger(1)) / fs; % duration in seconds

   % Get the significant strf
   [rfsig] = significant_strf(strf(i).rfcontra, p, n0, mdb, dur);

   % take the singular value decomposition to decompose
   % the strf into separable subunits
   [u,s,v] = svd(rfsig);
   singvals = sum(s);
   eigvals = singvals .^ 2;
size(u)
size(s)
size(v)
length(sum(s))
pause

   % Get phase-locking index
   pli = phase_locking_index(rfsig, mdb, w0, stim);

   % Get temporal correlation index
   deltat = t(2)-t(1);
   [rt, tshift] = temporal_correlation_index(rfsig, deltat);

   % Get spectral correlation index
   deltax = x(2)-x(1);
   [rx, xshift] = spectral_correlation_index(rfsig, deltax);

   % get strf energy at different power levels
   cumpower = 100 * cumsum(eigvals) / sum(eigvals+eps);

   rfsvd = zeros(size(rfsig,1), size(rfsig,2), length(percent_energy));

   for j = 1:length(percent_energy)

      % Get rid of noise for plotting appearance
      index_percent_power = find(cumpower <= percent_energy(j));
      num_sing_vals = min([max(index_percent_power) min([size(u) size(v) length(singvals)]) ]);

      rftemp = zeros(size(rfsig));
      for jj = 1:num_sing_vals
         rftemp = rftemp + singvals(jj) * u(:,jj) * v(:,jj)';
      end % (for k)

      % Get STRF energy
      [energy(j)] = strf_energy(rftemp, mdb, stim);

      % Save RFs at various SVD power levels
      rfsvd(:,:,j) = rftemp;

   end % (for j)


   dt = diff(t);
   dt = dt(1); % strf temporal resolution
   dx = diff(x);
   dx = dx(1); % strf spectral resolution

   maxtmf = ceil( 1 / dt );
   maxsmf = ceil( 1 / dx );

   ntbins = 100; %ceil(maxtmf / 20);
   ntbins = ntbins + ~mod(ntbins,2); % make odd number of bins

   nfbins = ceil(maxsmf / 0.15); % make resolution 0.15 cycles / octave
   nfbins = nfbins + ~mod(nfbins,2); % make odd number of bins

   dtfreq = maxtmf / ntbins; %size(rtftemp,2); % fft temporal frequency resolution
   dxfreq = maxsmf / nfbins; %size(rtftemp,1); % fft spectral frequency resolution


   for k = 1:size(rfsvd,3)

      % get the ripple transfer function/s
      rfk = fft2(rfsvd(:,:,k), nfbins, ntbins);
      rtftemp = fftshift(abs(rfk));

      % Get the tm frequency vector - will be in cycles per second
      if ( mod(size(rtftemp,2),2) )
         tmf = [-(size(rtftemp,2)-1)/2:(size(rtftemp,2)-1)/2]*dtfreq;
      else
         tmf = [-size(rtftemp,2)/2:(size(rtftemp,2)/2-1)]*dtfreq;
      end

      itmf0 = find(tmf==0);
      itmfp40 = find( tmf <= tm );
      ditmf = max(itmfp40)-itmf0;
      itmf = [itmf0-ditmf:itmf0+ditmf];
      tmf = tmf(itmf);

      % Get the sm frequency vector - will be in cycles per octave
      if ( mod(size(rtftemp,1),2) )
         xmf = [-(size(rtftemp,1)-1)/2:(size(rtftemp,1)-1)/2]*dxfreq;
      else
         xmf = [-size(rtftemp,1)/2:(size(rtftemp,1)/2-1)]*dxfreq;
      end

      ixmf0 = find(xmf == 0);
      ixmf4 = find( xmf <= sm );
      ixmf = ixmf0:max(ixmf4);
      xmf = xmf(ixmf);

      rtf(:,:,k) = rtftemp(ixmf, itmf);

   end % (for k)


   % assign parameters to output struct array
   params(i).exp = strf(i).exp;
   params(i).site = strf(i).site;
   params(i).chan = strf(i).chan;
   params(i).model = strf(i).model;
   params(i).depth = strf(i).depth;
   params(i).position = strf(i).position;
   params(i).stim = strf(i).stim;
   params(i).atten = strf(i).atten;
   params(i).sm = strf(i).sm;
   params(i).tm = strf(i).tm;
   params(i).mdb = strf(i).mdb;
   params(i).spl = strf(i).spl;
   params(i).n0 = n0;
   params(i).w0 = w0;
   params(i).percent_energy = percent_energy(:);
   params(i).rfenergy = energy(:);
   params(i).tmf = tmf;
   params(i).xmf = xmf;
   params(i).rtf = rtf;
   params(i).singvals = singvals;
   params(i).eigvals = eigvals;
   params(i).tci = [tshift(:) rt(:)];
   params(i).sci = [xshift(:) rx(:)];
   params(i).pli = pli;

end % (for i)


