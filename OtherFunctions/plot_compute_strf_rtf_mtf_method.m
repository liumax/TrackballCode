function plot_compute_strf_rtf_mtf_method(strf, trigger, n, flim, tlim)
%STRF_RTF_Method Plot STRF, RTF, and MTFs
%
% Function calls:
%
% plot_compute_strf_rtf_mtf_method(strf, trigger)
%     plot all STRFs. strf is a struct array holding the STRF data and 
%     trigger is a vector of trigger times. trigger is just used to
%     determine the duration of the stimulus.
%
% plot_compute_strf_rtf_mtf_method(strf, trigger, n)
%     plot the nth STRF in strf.
%
% plot_compute_strf_rtf_mtf_method(strf, trigger, flim)
%     plot all STRFs and change the frequency axis limits to flim = [flower fupper]
%     flim must be a 1X2 vector.
%     flim example:  flim = [2000 20000]
%
% plot_compute_strf_rtf_mtf_method(strf, trigger, n, flim)
%     plot the nth STRf with the frequency axis limits in flim
%
% plot_compute_strf_rtf_mtf_method(strf, trigger, flim, tlim)
%     plot all STRFs with the frequency limits specifice in flim 
%     and the time axis limits specified in tlim, in milliseconds.
%     tlim and flim must be 1X2 vectors.
%     tlim example:  tlim = [0 100] or tlim = [-20 100]
%
% plot_compute_strf_rtf_mtf_method(strf, trigger, n, flim, tlim)
%     plot the nth STRF in strf with the frequency limits specifice in flim 
%     and the time axis limits specified in tlim     
%
%   caa 10/28/02

error(nargchk(2,5,5));



% Parameters for plotting STRFs
if ( nargin == 2 )
   flim = [];
   tlim = [];
end

if ( nargin == 3 )
   if ( length(n) == 2 ) % flim was the 3rd input argument
      flim = sort(n);
   elseif ( length(n) == 1 ) % n was the third input argument
      strf = strf(n);
      flim = [];
   else
      error('Bad input arguments.');
   end
   tlim = [];
end

if ( nargin == 4 )
   if ( length(n) == 1 && length(flim) == 2 )
      strf = strf(n);
      flim = sort(flim);
      tlim = [];
   elseif ( length(n) == 2 && length(flim) == 2 )
      tlim = flim;
      flim = n;
      flim = sort(flim);
      tlim = sort(tlim);
   else
      error('Bad input arguments.');
   end
end

if ( nargin == 5 )
   if ( length(n) == 1 && length(flim) == 2 && length(tlim) == 2 )
      strf = strf(n);
      flim = sort(flim);
      tlim = sort(tlim);
   else
      error('Bad input arguments.');
   end
end



% Go through the STRFs and plot them
close all;
for i = 1:length(strf)

   % Define some initial parameters
   exp = strf(i).exp; % experiment date
   site = strf(i).site; % penetration number
   chan = strf(i).chan; % channel of michigan probe
   model = strf(i).model; % neuron number
   position = strf(i).position; % recording position from cortical surface
   p = 0.001; % significance level of STRF
   n0 = strf(i).n0contra; % number of spikes
   w0 = strf(i).w0contra; % firing rate
   sm = strf(i).sm; % maximum spectral modulation frequency of stimulus
   tm = strf(i).tm; % maximum temporal modulation frequency of ripple stimulus
   mdb = strf(i).mdb; % modulation depth of ripple
   stim = strf(i).stim; % ripple type
   fs = strf(i).fs; % A/D sampling rate
   taxis = strf(i).taxis; % time axis of STRF
   faxis = strf(i).faxis; % frequency axis of STRF
   t = strf(i).taxis; % time axis of STRF
   x = log2(strf(i).faxis ./ min(strf(i).faxis)); % get frequency in octaves
   dur = (trigger(end)-trigger(1)) / fs; % duration in seconds

   % Get the significant strf
   [rfsig] = significant_strf(strf(i).rfcontra, p, n0, mdb, dur);

   % Resolution of temporal and spectral axes
   dt = diff(t);
   dt = dt(1); % strf temporal resolution
   dx = diff(x);
   dx = dx(1); % strf spectral resolution in octaves

   maxtmf = ceil( 1 / dt ); % maximum possible modulation frequency based on STRF sampling
   maxsmf = ceil( 1 / dx );

%    dtmf = tm / 20; % desired resolution of Fourier transform
   dtmf = 2;
   ntbins = ceil(maxtmf/dtmf); % number of FFT bins needed for desired resolution
   ntbins = ntbins + ~mod(ntbins,2); % make odd number of bins

% Following is for general case. For cortical DMR, we always use 4 cyc/oct,
% so we want about 0.15 cyc/oct resolution in the RTF
%    dsmf = sm / 20; % resolution for spectral modulation axis
%    nfbins = ceil(maxsmf/dsmf); % number of needed bins
%    nfbins = nfbins + ~mod(nfbins,2); % make odd number of bins

   dsmf = 0.15;
   nfbins = ceil(maxsmf / dsmf); % make resolution 0.15 cycles / octave
   nfbins = nfbins + ~mod(nfbins,2); % make odd number of bins

   % Frequency resolution of the RTF/FFT
   dtfreq = maxtmf / ntbins; %size(rtftemp,2); % fft temporal frequency resolution
   dxfreq = maxsmf / nfbins; %size(rtftemp,1); % fft spectral frequency resolution


   % get the ripple transfer function/s
   rfk = fft2(rfsig, nfbins, ntbins); % 2D FFT
   rtftemp = fftshift(abs(rfk)); % shift FFT so it's centered at 0 frequency

   % Get the tm frequency vector - will be in cycles per second
   if ( mod(size(rtftemp,2),2) )
      tmf = [-(size(rtftemp,2)-1)/2:(size(rtftemp,2)-1)/2]*dtfreq;
   else
      tmf = [-size(rtftemp,2)/2:(size(rtftemp,2)/2-1)]*dtfreq;
   end

   % Get temporal modulation frequencies for RTF/MTF
   itmf0 = find(tmf==0);
   itmfp40 = find( tmf <= tm );
   ditmf = max(itmfp40)-itmf0;
   itmf = [itmf0-ditmf:itmf0+ditmf];
   rtf_tmf = tmf(itmf); % temporal axis for folded RTF


   % Get the sm frequency vector - will be in cycles per octave
   if ( mod(size(rtftemp,1),2) )
      xmf = [-(size(rtftemp,1)-1)/2:(size(rtftemp,1)-1)/2]*dxfreq;
   else
      xmf = [-size(rtftemp,1)/2:(size(rtftemp,1)/2-1)]*dxfreq;
   end

   % Get spectral modulation frequencies for RTF/MTF
   ixmf0 = find(xmf == 0);
   ixmf4 = find( xmf <= sm );
   ixmf = ixmf0:max(ixmf4);
   rtf_xmf = xmf(ixmf); % octave frequency axis for folded RTF

   % RTF for range of modulation frequencies used in stimulus
   rtf = rtftemp(ixmf, itmf);

   % RTF folded about 0 cyc/s.
   [rtf_fold, tmf_mtf, tmtf, smf_mtf, smtf] = rtf2mtf(rtf, rtf_tmf, rtf_xmf);

   % Parameters from MTFs
   [tbmf_temp, tbw6db_temp, tbw3db_temp, twc3db_temp] = mtf_bmf_bw_wc(tmf_mtf, tmtf);
   [sbmf_temp, sbw6db_temp, sbw3db_temp, swc3db_temp] = mtf_bmf_bw_wc(smf_mtf, smtf);



   %********************** Plot the results *******************


   % Plot the STRF
   figure;
   subplot(2,2,1);
   plot_strf_single(strf(i), trigger, flim, tlim);
   box off;


   % Plot the folded RTF. This will include the full range of spectral 
   % modulation frequencies but only positive temporal modulation 
   % frequencies
   subplot(2,2,2);
   imagesc(tmf_mtf, smf_mtf, rtf_fold);
   axis xy;
   box off;
   set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
   xlabel('Temp Mod Freq (cyc/sec)');
   ylabel('Spec Mod Freq (cyc/oct)');
   title(sprintf('%s %.0f %.0f-%.0f', exp, site, chan, model));
%    title(sprintf('%s %.0f %.0f-%.0f [%.0f %.0f] kHz', ...
%       exp, site, chan, model, min(flim)/1000, max(flim)/1000));
%    reds = brewmaps('reds', 21);
%    colormap(flipud(reds));



% The following code plots the frequency axis through the maximum value in
% the STRF. I have commented this out for concision.
%    % Take slice through the maximum value of the STRF
% %    figure;
%    subplot(3,2,2);
%    [irow, icol] = find(rfsig == max(max(rfsig)));
%    slice = rfsig(:, icol);
%    if ( nargin == 4 || nargin == 5 )
%       index = find(faxis > flim(1) & faxis < flim(2));
%       fslice = faxis(index);
%       slice = slice(index);
%    else
%       fslice = faxis;
%    end
%    slice = slice ./ max(slice); % Normalize so maximum value is 1
%    hold on;
%    plot(fslice, slice, 'k-');
%    set(gca,'xscale', 'log');
%    set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
%    ytick = 0:0.25:1;
%    set(gca,'ytick', ytick, 'yticklabel', ytick);
%    xlim([min(fslice) max(fslice)]);
%    ylim([1.1*min(slice) 1.1*max(slice)]);
%    xlabel('Frequency (cyc/sec)');
%    ylabel('Norm. Amp.');
%    title('Cochlear Place');

   % Plot the tMTF
   subplot(2,2,3);
   hold on;
   plot(tmf_mtf, tmtf, 'ko-', 'markerfacecolor', 'k', ...
      'markersize', 3);
   set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
   ytick = 0:0.25:1;
   set(gca,'ytick', ytick, 'yticklabel', ytick);
   xlabel('Temporal Modulation Freq (cyc/sec)');
   ylabel('Norm. Amp.');
   title('tMTF');

   % Plot the sMTF
%    figure;
   subplot(2,2,4);
   hold on;
   plot(smf_mtf, smtf, 'ko-', 'markerfacecolor', 'k', ...
      'markersize', 3);
   set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
   ylim([0 1.1]);
   xtick = 0:4;
   set(gca,'xtick', xtick, 'xticklabel', xtick);
   ytick = 0:0.25:1;
   set(gca,'ytick', ytick, 'yticklabel', ytick);
   xlabel('Spectral Modulation Freq (cyc/oct)');
   ylabel('Norm. Amp.');
   title('sMTF');

   % pause;
   set(gcf,'position', [170 230 500 400]);

end


return







