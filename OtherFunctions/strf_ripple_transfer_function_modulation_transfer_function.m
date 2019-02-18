function rtf_params = strf_ripple_transfer_function(params)
%plot_strf_parameters - plots best temporal, spectral modulation values
%   along with asymmetry indices and separability indices.
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
'xmf',                 [], ...
'tmtf',                [], ...
'xmtf',                [], ...
'best_tmf',            [], ...
'best_xmf',            []);





exp = params(1).exp;
site = params(1).site;
stim = params(1).stim;
atten = params(1).atten;
sm = params(1).sm;
tm = params(1).tm;
mdb = params(1).mdb;
spl = params(1).spl;

layer2 = 200;
layer3 = 400;
layer4 = 800;
layer5 = 1100;
layer6 = 1500;
whitematter = 2000;

layer1label = 100;
layer2label = 300;
layer3label = 600;
layer4label = 950;
layer5label = 1300;
layer6label = 1750;
whitematterlabel = 2200;



% temporal and spectral modulation frequency axes
tmf = params(1).tmf;
tmf = tmf(:)';
rtf_tmf = tmf;
ind0 = find(tmf==0);
tmf = tmf(ind0:end);

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


   chan = params(i).chan;
   model = params(i).model;
   depth = params(i).depth;
   position = params(i).position;
   n0 = params(i).n0;
   w0 = params(i).w0;

   % the ripple transfer function and its folded version
   rtf = params(i).rtf;
   rtf_right = rtf(:,ind0:end);
   rtf_left = rtf(:,1:ind0-1);
   rtfnew = zeros(size(rtf_right));
   rtfnew = rtfnew + rtf_right;
   rtfnew(:,2:end) = rtfnew(:,2:end) + fliplr(rtf_left);
   rtfnew(:,1) = 2*rtfnew(:,1);

   % asymmetry index
   asi = ( sum(sum(rtf(:,1:ind0-1))) - sum(sum(rtf(:,ind0+1:end))) ) ./ sum(sum(rtf));

   % calculate quadrant separability
   [uleft,sleft,vleft] = svd(rtf(:,1:ind0-1));
   singvals_left = sum(sleft);
   eigvals_left = singvals_left .^ 2;
   sepindex_left = 1 - eigvals_left(1) / sum(eigvals_left);

   [uright,sright,vright] = svd(rtf(:,ind0+1:end));
   singvals_right = sum(sright);
   eigvals_right = singvals_right .^ 2;
   sepindex_right = 1 - eigvals_right(1) / sum(eigvals_right);



   % temporal modulation parameters
   tmtf = sum(rtf,1);
   tmtf = tmtf(:)';
   tmtf_left = tmtf(1:ind0-1);
   tmtf_right = tmtf(ind0:end);
   tmtf = tmtf_right + [0 fliplr(tmtf_left)];
   tmtf(1) = 2 * tmtf(1);

   index = find(tmtf == max(tmtf));
   best_tmf = tmf(index);

   % spectral modulation parameters
   xmtf = sum(rtf,2)';
   xmtf = xmtf(:)';

   index = find(xmtf == max(xmtf));
   best_xmf = xmf(index);



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


   clf;

   subplot(2,2,1);
   imagesc(rtf);
   axis xy;
   title('Original RTF');

   subplot(2,2,2);
   imagesc(rtfnew);
   axis xy;
   title('Flipped and Summed RTF');

   subplot(2,2,3);
   hold on;
   plot(tmf, tmtf, 'ko-');
%    plot(tmf, tmtf_lowpass_fit, 'r-');
%    plot(tmf, tmtf_bandpass_fit, 'c-');
   xlabel('Temporal Modulation (cyc/sec)');
   title('RTF collapsed across freq');

   subplot(2,2,4);
   hold on;
   plot(xmf, xmtf, 'ko-');
%    plot(xmf, xmtf_lowpass_fit, 'r-');
%    plot(xmf, xmtf_bandpass_fit, 'c-');
   xlabel('Ripple Density (cyc/oct)');
   title('RTF collapsed across time');

   pause(1);

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
rtf_params(i).rtf_tmf = rtf_tmf;
rtf_params(i).rtf_xmf = rtf_xmf;
rtf_params(i).rtf = rtf;
rtf_params(i).rtf_asymmetry_index = asi;
rtf_params(i).rtf_left_eigvals = eigvals_left;
rtf_params(i).rtf_right_eigvals = eigvals_right;
rtf_params(i).rtf_left_sep_index = sepindex_left;
rtf_params(i).rtf_right_sep_index = sepindex_right;
rtf_params(i).rtf_folded = rtfnew;
rtf_params(i).tmf = tmf;
rtf_params(i).xmf = xmf;
rtf_params(i).tmtf = tmtf;
rtf_params(i).xmtf = xmtf;
rtf_params(i).best_tmf = best_tmf;
rtf_params(i).best_xmf = best_xmf;

end

% figure;
% suptitle(sprintf('%s site%.0f %s', exp, site, stim));
% 
% % First plot the separability index vs. depth
% subplot(2,2,1);
% hold on;
% plot(sepindex, position, 'ko', 'markerfacecolor','k');
% % plot(sepindex95, position, 'kv', 'markerfacecolor','k');
% box on;
% axis([0 1 0 2400]);
% plot([0 1],[layer2 layer2],'k:');
% plot([0 1],[layer3 layer3],'k:');
% plot([0 1],[layer4 layer4],'k:');
% plot([0 1],[layer5 layer5],'k:');
% plot([0 1],[layer6 layer6],'k:');
% plot([0 1],[whitematter whitematter],'k-');
% text([0.9],[layer1label],'L I');
% text([0.9],[layer2label],'L II');
% text([0.9],[layer3label],'L III');
% text([0.9],[layer4label],'L IV');
% text([0.9],[layer5label],'L V');
% text([0.9],[layer6label],'L VI');
% text([0.9],[whitematterlabel],'WM');
% 
% set(gca,'ytick',[0 layer2 layer3 layer4 layer5 layer6 whitematter]);
% xlabel('SI');
% ylabel('Depth (um)');
% title('Separability Index');
% set(gca,'ydir','rev');
% 
% 
% % plot the strf energy vs. depth
% subplot(2,2,2);
% hold on;
% plot(energy, position, 'ko', 'markerfacecolor','k');
% box on;
% axis([0 1 0 2400]);
% plot([0 1],[layer2 layer2],'k:');
% plot([0 1],[layer3 layer3],'k:');
% plot([0 1],[layer4 layer4],'k:');
% plot([0 1],[layer5 layer5],'k:');
% plot([0 1],[layer6 layer6],'k:');
% plot([0 1],[whitematter whitematter],'k-');
% text([0.9],[layer1label],'L I');
% text([0.9],[layer2label],'L II');
% text([0.9],[layer3label],'L III');
% text([0.9],[layer4label],'L IV');
% text([0.9],[layer5label],'L V');
% text([0.9],[layer6label],'L VI');
% text([0.9],[whitematterlabel],'WM');
% 
% set(gca,'ytick',[0 layer2 layer3 layer4 layer5 layer6 whitematter]);
% xlabel('Energy');
% ylabel('Depth (um)');
% title('STRF Energy');
% set(gca,'ydir','rev');
% 
% 
% 
% % plot phase-locking index vs. depth
% subplot(2,2,3);
% hold on;
% plot(pli, position, 'ko', 'markerfacecolor','k');
% box on;
% axis([0 0.5 0 2400]);
% plot([0 1],[layer2 layer2],'k:');
% plot([0 1],[layer3 layer3],'k:');
% plot([0 1],[layer4 layer4],'k:');
% plot([0 1],[layer5 layer5],'k:');
% plot([0 1],[layer6 layer6],'k:');
% plot([0 1],[whitematter whitematter],'k-');
% text([0.45],[layer1label],'L I');
% text([0.45],[layer2label],'L II');
% text([0.45],[layer3label],'L III');
% text([0.45],[layer4label],'L IV');
% text([0.45],[layer5label],'L V');
% text([0.45],[layer6label],'L VI');
% text([0.45],[whitematterlabel],'WM');
% 
% set(gca,'ytick',[0 layer2 layer3 layer4 layer5 layer6 whitematter]);
% xlabel('PLI');
% ylabel('Depth (um)');
% title('Phase-locking index');
% set(gca,'ydir','rev');
% 
% 
% 
% % plot temporal correlation index vs. depth
% subplot(2,2,4);
% hold on;
% plot(tci5, position, 'ko', 'markerfacecolor','k');
% plot(tci10, position, 'rs', 'markerfacecolor','r');
% plot(tci15, position, 'gv', 'markerfacecolor','g');
% legend('5ms', '10ms', '15ms', 0);
% box on;
% axis([0 1 0 2400]);
% plot([0 1],[layer2 layer2],'k:');
% plot([0 1],[layer3 layer3],'k:');
% plot([0 1],[layer4 layer4],'k:');
% plot([0 1],[layer5 layer5],'k:');
% plot([0 1],[layer6 layer6],'k:');
% plot([0 1],[whitematter whitematter],'k-');
% text([0.9],[layer1label],'L I');
% text([0.9],[layer2label],'L II');
% text([0.9],[layer3label],'L III');
% text([0.9],[layer4label],'L IV');
% text([0.9],[layer5label],'L V');
% text([0.9],[layer6label],'L VI');
% text([0.9],[whitematterlabel],'WM');
% % plot([0.75], [2175], 'rs', 'markerfacecolor', 'r');
% % text([0.78],[2175],'10ms');
% 
% set(gca,'ytick',[0 layer2 layer3 layer4 layer5 layer6 whitematter]);
% xlabel('TCI');
% ylabel('Depth (um)');
% title('Temporal correlation index');
% set(gca,'ydir','rev');
% 
% orient landscape
% 
% % print_mfilename(mfilename);
% 
% 
% 
% if ( 0 )
% 
% 
% 
% suptitle(sprintf('%s site%.0f %s', exp, site, stim));
% print_mfilename(mfilename);
% legend('All Evals', '95% Evals', 0);
% 
% % We also want to plot the 
% % (1) phase locking index
% % (2) temporal correlation index
% % (3) strf energy
% % (4) mtfs
% % (5) best temporal modulation frequency
% % (6) best spectral modulation frequency
% 
% 
% pos = [rfstats.position];
% bsm = [rfstats.bestsm];
% btm = [rfstats.besttm];
% bvel = [rfstats.bestvel];
% ai = [rfstats.ai];
% pval = [rfstats.aipval];
% si = [rfstats.si];
% 
% 
% subplot(3,2,5);
% hold on;
% plot(si, pos, 'ko', 'markerfacecolor','k');
% box on;
% axis([0 1 0 2500]);
% plot([0 1],[layer2 layer2],'k:');
% plot([0 1],[layer3 layer3],'k:');
% plot([0 1],[layer4 layer4],'k:');
% plot([0 1],[layer5 layer5],'k:');
% plot([0 1],[layer6 layer6],'k:');
% xlabel('SI');
% ylabel('Depth (um)');
% title('Separability Index');
% set(gca,'ydir','rev');
% 
% 
% subplot(3,2,1);
% hold on;
% plot(bsm,pos,'ko','markerfacecolor','k');
% xlim = get(gca,'xlim');
% axis([xlim(1) xlim(2) 0 2500]);
% plot([xlim(1) xlim(2)],[layer2 layer2],'k:');
% plot([xlim(1) xlim(2)],[layer3 layer3],'k:');
% plot([xlim(1) xlim(2)],[layer4 layer4],'k:');
% plot([xlim(1) xlim(2)],[layer5 layer5],'k:');
% plot([xlim(1) xlim(2)],[layer6 layer6],'k:');
% box on;
% set(gca,'ydir','rev');
% xlabel('Spec Mod (cyc/oct)');
% ylabel('Depth (um)');
% title('Best Spec Mod');
% 
% subplot(3,2,2);
% hold on;
% plot(btm,pos,'ko','markerfacecolor','k');
% xlim = get(gca,'xlim');
% temp = max(abs(xlim));
% axis([-temp temp 0 2500]);
% plot([-temp temp],[layer2 layer2],'k:');
% plot([-temp temp],[layer3 layer3],'k:');
% plot([-temp temp],[layer4 layer4],'k:');
% plot([-temp temp],[layer5 layer5],'k:');
% plot([-temp temp],[layer6 layer6],'k:');
% box on;
% set(gca,'ydir','rev');
% xlabel('Temp Mod (cyc/s)');
% title('Best Temp Mod');
% 
% 
% subplot(3,2,3);
% hold on;
% plot(bvel,pos,'ko','markerfacecolor','k');
% box on;
% xlim = get(gca,'xlim');
% axis([xlim(1) xlim(2) 0 2500]);
% plot([xlim(1) xlim(2)],[layer2 layer2],'k:');
% plot([xlim(1) xlim(2)],[layer3 layer3],'k:');
% plot([xlim(1) xlim(2)],[layer4 layer4],'k:');
% plot([xlim(1) xlim(2)],[layer5 layer5],'k:');
% plot([xlim(1) xlim(2)],[layer6 layer6],'k:');
% xlabel('Velocity (oct/s)');
% ylabel('Depth (um)');
% title('Best Velocity');
% set(gca,'ydir','rev');
% 
% subplot(3,2,4);
% hold on;
% plot(aisig, possig, 'ko', 'markerfacecolor','k');
% plot(ainotsig, posnotsig, 'ko', 'markerfacecolor','none');
% box on;
% maxai = max([abs(min(ai)) max(ai)]);
% xlo = -maxai-0.05;
% xhi = maxai+0.05;
% axis([xlo xhi 0 2500]);
% plot([xlo xhi],[layer2 layer2],'k:');
% plot([xlo xhi],[layer3 layer3],'k:');
% plot([xlo xhi],[layer4 layer4],'k:');
% plot([xlo xhi],[layer5 layer5],'k:');
% plot([xlo xhi],[layer6 layer6],'k:');
% legend('p<=0.05','p>0.05');
% xlabel('AI');
% title('Asymmetry Index');
% set(gca,'ydir','rev');
% 
% subplot(3,2,5);
% hold on;
% plot(si, pos, 'ko', 'markerfacecolor','k');
% box on;
% axis([0 1 0 2500]);
% plot([0 1],[layer2 layer2],'k:');
% plot([0 1],[layer3 layer3],'k:');
% plot([0 1],[layer4 layer4],'k:');
% plot([0 1],[layer5 layer5],'k:');
% plot([0 1],[layer6 layer6],'k:');
% xlabel('SI');
% ylabel('Depth (um)');
% title('Separability Index');
% set(gca,'ydir','rev');
% 
% suptitle(sprintf('%s site%.0f %s', ...
%          rfstats(1).exp, rfstats(1).site, rfstats(1).stim));
% 
% print_mfilename(mfilename);
% orient tall;
% 