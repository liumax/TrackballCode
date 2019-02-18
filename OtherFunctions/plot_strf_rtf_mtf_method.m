function plot_strf_rtf_mtf_method(params)
%
%   along with asymmetry indices and separability indices.
%
% params is obtained from the function strf_parameters.m
%
% rtf_params = strf_ripple_transfer_function(params)
%
%   caa 10/28/02



exp = params(1).exp;
site = params(1).site;
stim = params(1).stim;
atten = params(1).atten;
sm = params(1).sm;
tm = params(1).tm;
mdb = params(1).mdb;
spl = params(1).spl;


% temporal and spectral modulation frequency axes
tmf = params(1).tmf;
tmf = tmf(:)';
rtf_tmf = tmf;
ind0 = find(tmf==0);

xmf = params(1).xmf;
xmf = xmf(:)';
rtf_xmf = xmf;

% functions used for curve fitting. Define as inline functions
mtf_lowpass = inline('k(1) * (1 + (f ./ k(2)).^k(3) ).^-1 + k(4)', 'k', 'f');
% mtf_bandpass = inline('k(1) * (1 + (k(2) ./ f).^k(3) ).^-1 .* (1 + (f ./ k(4)).^k(5) ).^-1 + k(6)', 'k', 'f');
mtf_bandpass = inline('k(1) * ( exp(k(2).*f) + exp(k(3).*f) ) + k(4)', 'k', 'f');

% curve fitting options
options = optimset('TolFun', 1e-8, 'TolX', 1e-8, 'MaxFunEvals', 1000, ...
                   'MaxIter', 1000, 'LevenbergMarquardt', 'on');


% First get all the separability indices

for i = 1:length(params)

   % recording basic information
   chan = params(i).chan;
   model = params(i).model;
   depth = params(i).depth;
   position = params(i).position;
   n0 = params(i).n0;
   w0 = params(i).w0;
   percent_energy = params(i).percent_energy;

   for j = length(percent_energy):length(percent_energy)

      % the ripple transfer function and its folded version
      rtf = params(i).rtf(:,:,j);

      % asymmetry index - power in left versus right quadrants
      asi = ( sum(sum(rtf(:,1:ind0-1))) - sum(sum(rtf(:,ind0+1:end))) ) ./ sum(sum(rtf));

      % calculate quadrant separability - left quadrant
      [uleft, sleft, vleft] = svd(rtf(:,1:ind0-1));
      svals_left = sum(sleft);
      evals_left = svals_left .^ 2;

      eigvals_left = evals_left;
      sepindex_left = evals_left(1) / sum(evals_left);

      % calculate quadrant separability - right quadrant
      [uright, sright, vright] = svd(rtf(:,ind0+1:end));
      svals_right = sum(sright);
      evals_right = svals_right .^ 2;

      eigvals_right = evals_right;
      sepindex_right = evals_right(1) / sum(evals_right);

      [rtf_fold, tmf_mtf, tmtf, smf_mtf, smtf] = rtf2mtf(rtf, rtf_tmf, rtf_xmf);

      [tbmf_temp, tbw6db_temp, tbw3db_temp, twc3db_temp] = mtf_bmf_bw_wc(tmf_mtf, tmtf);
      [sbmf_temp, sbw6db_temp, sbw3db_temp, swc3db_temp] = mtf_bmf_bw_wc(smf_mtf, smtf);

      figure;

%       subplot(4,1,1);
%       imagesc(rtf_tmf, rtf_xmf, rtf);
%       axis xy;
%       set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
%       xlabel('Temp Mod Freq (cyc/sec)');
%       ylabel('Spec Mod Freq (cyc/oct)');
%       title('RTF');

      subplot(2,2,2);
      imagesc(tmf_mtf, smf_mtf, rtf_fold);
      axis xy;
      set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
      xlabel('Temp Mod Freq (cyc/sec)');
      ylabel('Spec Mod Freq (cyc/oct)');
      title(sprintf('%s %.0f %.0f-%.0f %.0f um\nRTF', ...
         exp, site, chan, model, position));

      subplot(2,2,3);
      hold on;
      plot(tmf_mtf, tmtf, 'ko-', 'markerfacecolor', 'k', ...
      'markersize', 2);
      plot([tbmf_temp tbmf_temp], [0 1], 'k-');
      set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
      ytick = 0:0.25:1;
      set(gca,'ytick', ytick, 'yticklabel', ytick);
      xlabel('Temporal Modulation Freq (cyc/sec)');
      ylabel('Norm. Amp.');
      title('tMTF');

      subplot(2,2,4);
      hold on;
      plot(smf_mtf, smtf, 'ko-', 'markerfacecolor', 'k', ...
      'markersize', 2);
      plot([sbmf_temp sbmf_temp], [0 1], 'k-');
      set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
      ytick = 0:0.25:1;
      set(gca,'ytick', ytick, 'yticklabel', ytick);
      xlabel('Spectral Modulation Freq (cyc/oct)');
      ylabel('Norm. Amp.');
      title('sMTF');

%       pause;
      set(gcf,'position', [680   169   308   809]);

   end % (for j)


end

return



%    [tmtf_lowpass_fit_params, r2_tmtf_lowpass_fit] = lsqcurvefit(mtf_lowpass, [50 10 2 100], tmf, tmtf, [0 0 0 0], [inf inf inf inf], options);
%    [tmtf_bandpass_fit_params, r2_tmtf_bandpass_fit] = lsqcurvefit(mtf_bandpass, [50 10 2 10 2 100], tmf, tmtf, [0 0 0 0 0 0], [inf inf inf inf inf inf], options);
% 
%    tmtf_lowpass_fit = mtf_lowpass(tmtf_lowpass_fit_params, tmf);
%    tmtf_bandpass_fit = mtf_bandpass(tmtf_bandpass_fit_params, tmf);


%    [xmtf_lowpass_fit_params, r2_xmtf_lowpass_fit] = lsqcurvefit(mtf_lowpass, [max(xmtf) 1 2 0], xmf, xmtf, [0 0 0 0], [inf 4 inf inf], options);
%    [xmtf_bandpass_fit_params, r2_xmtf_bandpass_fit] = lsqcurvefit(mtf_bandpass, [max(xmtf) 1 2 1 2 0], xmf, xmtf, [0 0 0 0 0 0], [inf 4 inf 4 inf inf], options);
% 
%    xmtf_lowpass_fit = mtf_lowpass(xmtf_lowpass_fit_params, xmf);
%    xmtf_bandpass_fit = mtf_bandpass(xmtf_bandpass_fit_params, xmf);
% 
%    tmtf_lowpass_fit_params
%    r2_tmtf_lowpass_fit
% 
%    tmtf_bandpass_fit_params
%    r2_tmtf_bandpass_fit
% 
%    xmtf_lowpass_fit_params
%    r2_xmtf_lowpass_fit
% 
%    xmtf_bandpass_fit_params
%    r2_xmtf_bandpass_fit









