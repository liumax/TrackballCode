function rtf_params = strf_ripple_transfer_function(params)
%
%   along with asymmetry indices and separability indices.
%
% params is obtained from the function strf_parameters.m
%
% rtf_params = strf_ripple_transfer_function(params)
%
%   caa 10/28/02


rtf_params = struct(...
'exp',                 [], ...
'site',                [], ...
'chan',                [], ...
'model',               [], ...
'depth',               [], ...
'position',            [], ...
'stim',                [], ...
'atten',               [], ...
'sm',                  [], ...
'tm',                  [], ...
'mdb',                 [], ...
'spl',                 [], ...
'n0',                  [], ...
'w0',                  [], ...
'percent_energy',      [], ...
'rtf_tmf',             [], ...
'rtf_xmf',             [], ...
'rtf',                 [], ...
'rtf_asymmetry_index', [], ...
'rtf_left_eigvals',    [], ...
'rtf_right_eigvals',   [], ...
'rtf_left_sep_index',  [], ...
'rtf_right_sep_index', [], ...
'rtf_folded',          [], ...
'tmf',                 [], ...
'tmtf',                [], ...
'best_tmf',            [], ...
'tbw6db',              [], ...
'tbw3db',              [], ...
'twc3db',              [], ...
'xmf',                 [], ...
'xmtf',                [], ...
'best_xmf',            [], ...
'sbw6db',              [], ...
'sbw3db',              [], ...
'swc3db',              []);



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

   % initialize data variables
   asi = cell(size(percent_energy));
   eigvals_left = cell(size(percent_energy));
   sepindex_left = cell(size(percent_energy));
   eigvals_right = cell(size(percent_energy));
   sepindex_right = cell(size(percent_energy));

   rtf_cell = cell(size(percent_energy));
   rtf_folded_cell = cell(size(percent_energy));

   tmf_cell = cell(size(percent_energy));
   tmtf_cell = cell(size(percent_energy));

   tmbf_cell = cell(size(percent_energy));
   tbw6db_cell = cell(size(percent_energy));
   tbw3db_cell = cell(size(percent_energy));
   twc3db_cell = cell(size(percent_energy));


   smf_cell = cell(size(percent_energy));
   smtf_cell = cell(size(percent_energy));

   smbf_cell = cell(size(percent_energy));
   sbw6db_cell = cell(size(percent_energy));
   sbw3db_cell = cell(size(percent_energy));
   swc3db_cell = cell(size(percent_energy));

   for j = 1:length(percent_energy)

      % the ripple transfer function and its folded version
      rtf = params(i).rtf(:,:,j);

      % asymmetry index - power in left versus right quadrants
      asi{j} = ( sum(sum(rtf(:,1:ind0-1))) - sum(sum(rtf(:,ind0+1:end))) ) ./ sum(sum(rtf));

      % calculate quadrant separability - left quadrant
      [uleft, sleft, vleft] = svd(rtf(:,1:ind0-1));
      svals_left = sum(sleft);
      evals_left = svals_left .^ 2;

      eigvals_left{j} = evals_left;
      sepindex_left{j} = evals_left(1) / sum(evals_left);

      % calculate quadrant separability - right quadrant
      [uright, sright, vright] = svd(rtf(:,ind0+1:end));
      svals_right = sum(sright);
      evals_right = svals_right .^ 2;

      eigvals_right{j} = evals_right;
      sepindex_right{j} = evals_right(1) / sum(evals_right);

      [rtf_fold, tmf_mtf, tmtf, smf_mtf, smtf] = rtf2mtf(rtf, rtf_tmf, rtf_xmf);

      rtf_cell{j} = rtf;
      rtf_folded_cell{j} = rtf_fold;

      tmf_cell{j} = tmf_mtf;
      tmtf_cell{j} = tmtf;

      smf_cell{j} = smf_mtf;
      smtf_cell{j} = smtf;

      [tbmf_temp, tbw6db_temp, tbw3db_temp, twc3db_temp] = mtf_bmf_bw_wc(tmf_mtf, tmtf);
      [sbmf_temp, sbw6db_temp, sbw3db_temp, swc3db_temp] = mtf_bmf_bw_wc(smf_mtf, smtf);

      pause

      tbmf_cell{j} = tbmf_temp;
      tbw6db_cell{j} = tbw6db_temp;
      tbw3db_cell{j} = tbw3db_temp;
      twc3db_cell{j} = twc3db_temp;

      sbmf_cell{j} = sbmf_temp;
      sbw6db_cell{j} = sbw6db_temp;
      sbw3db_cell{j} = sbw3db_temp;
      swc3db_cell{j} = swc3db_temp;


      clf;

      subplot(2,2,1);
      imagesc(rtf);
      axis xy;
      title('Original RTF');

      subplot(2,2,2);
      imagesc(rtf_fold);
      axis xy;
      title('Flipped and Summed RTF');

      subplot(2,2,3);
      hold on;
      plot(tmf_mtf, tmtf, 'ko-');
%       plot(tmf, tmtf_lowpass_fit, 'r-');
%       plot(tmf, tmtf_bandpass_fit, 'c-');
      xlabel('Temporal Modulation (cyc/sec)');
      title('RTF collapsed across freq');

      subplot(2,2,4);
      hold on;
      plot(smf_mtf, smtf, 'ko-');
%    plot(xmf, xmtf_lowpass_fit, 'r-');
%    plot(xmf, xmtf_bandpass_fit, 'c-');
      xlabel('Ripple Density (cyc/oct)');
      title('RTF collapsed across time');

      pause(1);

   end % (for j)

   rtf_params(i).exp = exp;
   rtf_params(i).site = site;
   rtf_params(i).chan = chan;
   rtf_params(i).model = model;
   rtf_params(i).depth = depth;
   rtf_params(i).position = position;
   rtf_params(i).stim = stim;
   rtf_params(i).atten = atten;
   rtf_params(i).sm = sm;
   rtf_params(i).tm = tm;
   rtf_params(i).mdb = mdb;
   rtf_params(i).spl = spl;
   rtf_params(i).n0 = n0;
   rtf_params(i).w0 = w0;
   rtf_params(i).percent_energy = percent_energy;
   rtf_params(i).rtf_tmf = rtf_tmf;
   rtf_params(i).rtf_xmf = rtf_xmf;
   rtf_params(i).rtf = rtf_cell;
   rtf_params(i).rtf_asymmetry_index = asi;
   rtf_params(i).rtf_left_eigvals = eigvals_left;
   rtf_params(i).rtf_right_eigvals = eigvals_right;
   rtf_params(i).rtf_left_sep_index = sepindex_left;
   rtf_params(i).rtf_right_sep_index = sepindex_right;
   rtf_params(i).rtf_folded = rtf_folded_cell;

   rtf_params(i).tmf = tmf_cell;
   rtf_params(i).tmtf = tmtf_cell;
   rtf_params(i).best_tmf = tbmf_cell;
   rtf_params(i).tbw6db = tbw6db_cell;
   rtf_params(i).tbw3db = tbw3db_cell;
   rtf_params(i).twc3db = twc3db_cell;

   rtf_params(i).xmf = smf_cell;
   rtf_params(i).xmtf = smtf_cell;
   rtf_params(i).best_xmf = sbmf_cell;
   rtf_params(i).sbw6db = sbw6db_cell;
   rtf_params(i).sbw3db = sbw3db_cell;
   rtf_params(i).swc3db = swc3db_cell;

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









