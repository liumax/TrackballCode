function [rtmtf, rsmtf] = strf_threshold_rtf_mtf(strf, trigger, pval)
% strf_threshold_rtf_mtf(strf, trigger)
% -----------------------------------------------
%
% Calculate RTFs and MTFs for STRFs that have been
% thresholded at different significance levels.
%
% This will help us understand the effect of the 
% significance testing on the RTF and MTF estimates.
%
%
% strf : struct array holding the strf data
%
% trigger : vector of ripple stimulus trigger times, in 
% units of sample number
%
% pval : signficance level at which to threshold the STRF. Can be
% a vector, in which case all elements in pval will be used
%
%   caa 12/2/09


% parameters to obtain:
% best tm
% best sm
% mtfs: tm and sm
% strf energy - use function strf_energy.m

close all;

percent_energy = [80 85 90 95 97.5 100];
energy = zeros(size(percent_energy));

rtmtf = [];
rsmtf = [];


for i = 1:length(strf)

   tmtf_mat = [];
   smtf_mat = [];

   % Define some initial parameters
   n0 = strf(i).n0contra;
   w0 = strf(i).w0contra;
   sm = strf(i).sm;
   tm = strf(i).tm;
   mdb = strf(i).mdb;
   stim = strf(i).stim;
   fs = strf(i).fs;
   t = strf(i).taxis;
   f = strf(i).faxis;
   x = log2(strf(i).faxis ./ min(strf(i).faxis));
   dur = (trigger(end)-trigger(1)) / fs; % duration in seconds
   rfcontra = strf(i).rfcontra;


   for j = 1:length(pval)

      p = pval(j);

      % Get the significant strf
      [rfsig] = significant_strf(rfcontra, p, n0, mdb, dur);

      dt = diff(t);
      dt = dt(1); % strf temporal resolution
      dx = diff(x);
      dx = dx(1); % strf spectral resolution

      maxtmf = ceil( 1 / dt );
      maxsmf = ceil( 1 / dx );

      ntbins = ceil(maxtmf / 1); % this will give us 1 Hz tmtf resolution
      nfbins = ceil(maxsmf / 0.15); % make resolution 0.15 cycles / octave

      dtfreq = maxtmf / ntbins; %size(rtftemp,2); % fft temporal frequency resolution
      dxfreq = maxsmf / nfbins; %size(rtftemp,1); % fft spectral frequency resolution

      % get the ripple transfer function/s
      rfk = fft2(rfsig, nfbins, ntbins);
      rtftemp = fftshift(abs(rfk));


      % Get the tm frequency vector - will be in cycles per second
      if ( mod(size(rtftemp,2),2) )
         tmf = [-(size(rtftemp,2)-1)/2:(size(rtftemp,2)-1)/2]*dtfreq;
      else
         tmf = [-size(rtftemp,2)/2:(size(rtftemp,2)/2-1)]*dtfreq;
      end

      itmf0 = find(tmf==0);
      itmfp40 = find(tmf<=40);
      ditmf = max(itmfp40)-itmf0;
      itmf = [itmf0-ditmf:itmf0+ditmf];
      tmf_rtf = tmf(itmf);

      % Get the sm frequency vector - will be in cycles per octave
      if ( mod(size(rtftemp,1),2) )
         xmf = [-(size(rtftemp,1)-1)/2:(size(rtftemp,1)-1)/2]*dxfreq;
      else
         xmf = [-size(rtftemp,1)/2:(size(rtftemp,1)/2-1)]*dxfreq;
      end

      ixmf0 = find(xmf == 0);
      ixmf4 = find(xmf <= 4);
      ixmf = ixmf0:max(ixmf4);
      xmf_rtf = xmf(ixmf);

      rtf = rtftemp(ixmf, itmf);

      [rtf_fold, tmf_mtf, tmtf, smf_mtf, smtf] = rtf2mtf(rtf, tmf_rtf, xmf_rtf);


      figure;

      % STRF
      % ------------------
      subplot(2,3,1);
      minmin = min(min(rfsig));
      maxmax = max(max(rfsig));
      boundary = max([abs(minmin) abs(maxmax)]);

      imagesc(t*1000,f / 1000,single(rfsig));
      set(gca,'ydir', 'normal');
      set(gca, 'tickdir', 'out', 'ticklength', [0.025 0.025]);
      set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
      title(sprintf('P value = %.5f', p));


      % RTF 
      % ------------------
      subplot(2,3,2);
      minmin = 0; %min(min(rfsig));
      maxmax = max(max(rtf));
      boundary = max([minmin abs(maxmax)]);

      imagesc( tmf_rtf, xmf_rtf, rtf );
      set(gca,'ydir', 'normal');
      set(gca, 'tickdir', 'out', 'ticklength', [0.025 0.025]);
      set(gca, 'clim', [minmin maxmax]);

%    xtick = get(gca,'xtick');
%    xticklabel = round(tmf_rtf(xtick)*10) / 10;
%    set(gca,'xticklabel', xticklabel);
% 
%    ytick = get(gca,'ytick');
%    yticklabel = round(xmf_rtf(ytick)*10) / 10;
%    set(gca,'yticklabel', yticklabel);


      % Folded RTF 
      % ------------------
      subplot(2,3,3);
      imagesc(tmf_mtf, smf_mtf, rtf_fold);
      set(gca, 'tickdir', 'out', 'ticklength', [0.025 0.025]);
      axis xy;
      title('Folded and Summed RTF');


      % temporal MTF
      % ------------------
      subplot(2,3,4);
      hold on;
      plot(tmf_mtf, tmtf, 'ko-');
      set(gca, 'tickdir', 'out', 'ticklength', [0.025 0.025]);
      xlabel('Temp Mod Freq (cyc/sec)');
      title('tMTF');


      % spectral MTF
      % ------------------
      subplot(2,3,5);
      hold on;
      plot(smf_mtf, smtf, 'ko-');
      set(gca, 'tickdir', 'out', 'ticklength', [0.025 0.025]);
      xlabel('Spec Mod Freq (cyc/oct)');
      title('sMTF');


      tmtf_mat = [tmtf_mat; tmtf];
      smtf_mat = [smtf_mat; smtf];

   end % (for j)


   if ( length(pval) > 1 )

      cmb = nchoosek(1:length(pval),2);
      [nrows, ncols] = size(cmb);

      temp = [];
      rtmtf_temp = [];
      rsmtf_temp = [];

      for k = 1:nrows

         temp = corrcoef( tmtf_mat( cmb(k,1),: ), tmtf_mat( cmb(k,2), : ) );
         rtmtf_temp(k) = temp(1,2);

         temp = corrcoef( smtf_mat( cmb(k,1),: ), smtf_mat( cmb(k,2), : ) );
         rsmtf_temp(k) = temp(1,2);

      end % (for k)

      rtmtf = [rtmtf; rtmtf_temp];
      rsmtf = [rsmtf; rsmtf_temp];

   end % (if)


end % (for i)

return;



