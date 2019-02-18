function plot_strf_mtf_ptparams(rtf_params, fra_params)
%plot_freq_resp_area_parameters - plots best temporal, spectral modulation values
%   along with asymmetry indices and separability indices.
%
% rtf_params : holds ripple transfer function parameters. 'rtf_params' is usually
% found in a file such as 2004-1-14-site2-2352um-30db-dmr1-fs18115-rtf-params.mat
%
% rtf_params is obtained from the function
%
%    rtf_params = strf_ripple_transfer_function(params)
%
%   caa 10/28/02

exp = rtf_params(1).exp;
site = rtf_params(1).site;
stim = rtf_params(1).stim;

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



% initialize the variables
fra_position = zeros(1,length(fra_params));
cf = zeros(1,length(fra_params));
threshold = zeros(1,length(fra_params));
latency = zeros(1,length(fra_params));
q10 = zeros(1,length(fra_params));
q20 = zeros(1,length(fra_params));
q30 = zeros(1,length(fra_params));
q40 = zeros(1,length(fra_params));




ptposition = [ptparams.position];
cf = [ptparams.cf];
q = [ptparams.q];
latency = [ptparams.latency];
fr = [ptparams.w0];







% initialize the variables
position = zeros(1,length(rtf_params));
w0 = zeros(1,length(rtf_params));
asi = zeros(1,length(rtf_params));
left_sep_index = zeros(1,length(rtf_params));
right_sep_index = zeros(1,length(rtf_params));
best_tmf = zeros(1,length(rtf_params));
best_smf = zeros(1,length(rtf_params));

tmf = rtf_params(1).tmf;
xmf = rtf_params(1).xmf;

tmtf = zeros(length(rtf_params),length(tmf));
xmtf = zeros(length(rtf_params),length(xmf));

% First get all the separability indices

for i = 1:length(rtf_params)

   position(i) = rtf_params(i).position;
   w0(i) = rtf_params(i).w0;
   asi(i) = rtf_params(i).rtf_asymmetry_index(end);
   left_sep_index(i) = rtf_params(i).rtf_left_sep_index(end);
   right_sep_index(i) = rtf_params(i).rtf_right_sep_index(end);
%    best_tmf(i) = rtf_params(i).best_tmf(end);
%    best_xmf(i) = rtf_params(i).best_xmf(end);


   if ( iscell(rtf_params(i).tmtf) )
      tmtf_temp = rtf_params(i).tmtf{end};
      xmtf_temp = rtf_params(i).xmtf{end};
   else
      tmtf_temp = rtf_params(i).tmtf;
      xmtf_temp = rtf_params(i).xmtf;
   end

   tmtf(i,:) = tmtf_temp ./ max(tmtf_temp);
   xmtf(i,:) = xmtf_temp ./ max(xmtf_temp);

end



% Process the temporal modulation transfer functions
twidth6db = zeros(1,length(position));
twidth3db = zeros(1,length(position));
twc3db = zeros(1,length(position));
tbp3db = [];
tbp6db = [];
tlp3db = [];
tlp6db = [];

for i = 1:length(position)

   tfit = min(tmf):0.001:max(tmf);
   mtf = interp1(tmf, tmtf(i,:), tfit, 'spline');

%    plot(tmf, tmtf(i,:), 'ro-', tfit, mtf, 'g-');
%    axis([-0.5 40.5 0 1.1]);
%    title(sprintf('i = %.0f',i));

   indmax = find(mtf == max(mtf));

   ind3dblow = find(mtf(1:indmax) >= 0.7);
   ind3dbhi = find(mtf(indmax:end) >= 0.7);

   ind6dblow = find(mtf(1:indmax) >= 0.5);
   ind6dbhi = find(mtf(indmax:end) >= 0.5);
   
   twidth3db(i) = tfit(max(ind3dbhi+indmax-1)) - tfit(min(ind3dblow));
   twidth6db(i) = tfit(max(ind6dbhi+indmax-1)) - tfit(min(ind6dblow));
   twc3db(i) = tfit(max(ind3dbhi+indmax-1));

   if ( mtf(min(ind3dblow)) < 0.72 & mtf(max(ind3dbhi+indmax-1)) < 0.72 ) % it's bandpass
      itemp = find(mtf == max(mtf));
      best_tmf(i) = tfit(itemp);
      tbp3db = [tbp3db position(i)];
   else  % it's lowpass
      best_tmf(i) = tfit(max(ind3dbhi+indmax-1)) ./ 2;
      tlp3db = [tlp3db position(i)];
   end

   if ( mtf(min(ind6dblow)) < 0.52 & mtf(max(ind6dbhi+indmax-1)) < 0.52 )
      tbp6db = [tbp6db position(i)];
   else
      tlp6db = [tlp6db position(i)];
   end

end % (for i)



% Process the spectral modulation functions
xwidth6db = zeros(1,length(position));
xwidth3db = zeros(1,length(position));
xwc3db = zeros(1,length(position));
xbp3db = [];
xbp6db = [];
xlp3db = [];
xlp6db = [];

for i = 1:length(position)

   xfit = min(xmf):0.001:max(xmf);
   mtf = interp1(xmf, xmtf(i,:), xfit, 'spline');

%    plot(xmf, xmtf(i,:), 'ro-', xfit, mtf, 'g-');
%    axis([-0.05 4.05 0 1.1]);

   indmax = find(mtf == max(mtf));

   ind3dblow = find(mtf(1:indmax) >= 0.7);
   ind3dbhi = find(mtf(indmax:end) >= 0.7);

   ind6dblow = find(mtf(1:indmax) >= 0.5);
   ind6dbhi = find(mtf(indmax:end) >= 0.5);
   
   xwidth3db(i) = xfit(max(ind3dbhi+indmax-1)) - xfit(min(ind3dblow));
   xwidth6db(i) = xfit(max(ind6dbhi+indmax-1)) - xfit(min(ind6dblow));
   xwc3db(i) = xfit(max(ind3dbhi+indmax-1));

   if ( mtf(min(ind3dblow)) < 0.72 & mtf(max(ind3dbhi+indmax-1)) < 0.72 )  % it's bandpass
      istemp = find(mtf == max(mtf));
      best_smf(i) = xfit(istemp);
      xbp3db = [xbp3db position(i)];
   else  % it's lowpass
      best_smf(i) = xfit(max(ind3dbhi+indmax-1)) ./ 2;
      xlp3db = [xlp3db position(i)];
   end

   if ( mtf(min(ind6dblow)) < 0.52 & mtf(max(ind6dbhi+indmax-1)) < 0.52 )
      xbp6db = [xbp6db position(i)];
   else
      xlp6db = [xlp6db position(i)];
   end

end % (for i)




figure;

suptitle(sprintf('%s site%.0f %s', exp, site, stim));

subplot(2,3,1);
imagesc(tmtf);
xlabel('TMTF (cyc/s)');
upos = unique(position);
for i = 1:length(upos)
   index(i) = min(find(upos(i) == position));
end
set(gca,'xtick',1:2:length(tmf),'xticklabel',tmf(1:2:end));
set(gca,'ytick',index,'yticklabel',upos);
set(gca,'fontsize', 8)


subplot(2,3,2);
f = round(xmf * 10)/10;
imagesc(xmtf);
xlabel('XMTF (cyc/oct)');
set(gca,'xtick',1:4:length(xmf),'xticklabel',(f(1:4:end)));
set(gca,'ytick',index,'yticklabel',upos);
set(gca,'fontsize', 8)


% First plot the cf vs. depth
subplot(2,3,3);
hold on;
plot(cf, ptposition, 'ko', 'markerfacecolor','k');
set(gca,'xscale', 'log');
box on;
axis([1 50 0 2400]);
plot([1 50],[layer2 layer2],'k:');
plot([1 50],[layer3 layer3],'k:');
plot([1 50],[layer4 layer4],'k:');
plot([1 50],[layer5 layer5],'k:');
plot([1 50],[layer6 layer6],'k:');
plot([1 50],[whitematter whitematter],'k-');
text([30],[layer1label],'L I');
text([30],[layer2label],'L II');
text([30],[layer3label],'L III');
text([30],[layer4label],'L IV');
text([30],[layer5label],'L V');
text([30],[layer6label],'L VI');
text([30],[whitematterlabel],'WM');

set(gca,'ytick',[0 layer2 layer3 layer4 layer5 layer6 whitematter]);
xlabel('CF (kHz)');
title('CF');
set(gca,'ydir','rev');



% Plot temporal MTF width vs. depth
% subplot(2,3,4);
% hold on;
% posdata = min(min(upos)):max(max(upos));
% fitmodel = fit(position(:), twidth6db(:), 'smoothingspline', 'SmoothingParam', 0.5e-5);
% [twdata] = feval(fitmodel, posdata);
% plot(twdata, posdata, 'r-', 'linewidth', 2);
% plot(twidth6db, position, 'ko', 'markerfacecolor','k');
% % plot(twidth3db, position, 'ro', 'markerfacecolor','r');
% % legend('6db','3db');
% box on;
% lax = -0.1 * max(twidth6db);
% rax = 1.1 * max(twidth6db);
% axis([lax rax 0 2400]);
% plot([lax rax],[layer2 layer2],'k:');
% plot([lax rax],[layer3 layer3],'k:');
% plot([lax rax],[layer4 layer4],'k:');
% plot([lax rax],[layer5 layer5],'k:');
% plot([lax rax],[layer6 layer6],'k:');
% plot([lax rax],[whitematter whitematter],'k-');
% set(gca,'ytick',[0 layer2 layer3 layer4 layer5 layer6 whitematter]);
% xlabel('TMTF Width (cyc/s)');
% set(gca,'ydir','rev');


subplot(2,3,4);
tposbins = 0:300:2250;
tnbp = hist(tbp3db, tposbins);
tnlp = hist(tlp3db, tposbins);
num = sum([tnbp; tnlp]);
tnbp ./ num;
barh(tposbins, [tnlp(:) tnbp(:)], 'stacked');
% barh(tposbins, tnbp ./ num);
% axis([0 1.05 -250 2450]);
axis([0 max(max([tnbp+tnlp]))+1 -250 2450]);
% set(gca,'xticklabel', tposbins/100);
% set(gca,'fontsize', 8)
children = get(gca,'children');
set(children(1),'facecolor', [0.75 0.75 0.75]);
set(children(2),'facecolor', [0 0 0]);
legend('LowPass', 'BandPass');
set(gca, 'ydir', 'rev');


% Plot spectral MTF width vs. depth
subplot(2,3,5);
hold on;
posdata = min(min(upos)):max(max(upos));
fitmodel = fit(position(:), xwidth3db(:), 'smoothingspline', 'SmoothingParam', 0.5e-5);
[xwdata] = feval(fitmodel, posdata);
plot(xwdata, posdata, 'r-', 'linewidth', 2);
% plot(xwidth6db, position, 'ko', 'markerfacecolor','k');
plot(xwidth3db, position, 'ko', 'markerfacecolor','k');
box on;
lax = -0.1 * max(xwidth6db);
rax = 1.1 * max(xwidth6db);
axis([lax rax 0 2400]);
plot([lax rax],[layer2 layer2],'k:');
plot([lax rax],[layer3 layer3],'k:');
plot([lax rax],[layer4 layer4],'k:');
plot([lax rax],[layer5 layer5],'k:');
plot([lax rax],[layer6 layer6],'k:');
plot([lax rax],[whitematter whitematter],'k-');
set(gca,'ytick',[0 layer2 layer3 layer4 layer5 layer6 whitematter]);
xlabel('XMTF Width (cyc/oct)');
set(gca,'ydir','rev');


% plot sharpness of tuning vs. depth
subplot(2,3,6);
hold on;
posdata = min(min(ptposition)):max(max(ptposition));
fitmodel = fit(ptposition(:), q(:), 'smoothingspline', 'SmoothingParam', 0.5e-5);
[qdata] = feval(fitmodel, posdata);
plot(qdata, posdata, 'r-', 'linewidth', 2);
plot(q, fra_position, 'ko', 'markerfacecolor','k');
% plot(q20, fra_position, 'rs', 'markerfacecolor','r');
% plot(q30, fra_position, 'gv', 'markerfacecolor','g');
% plot(q40, fra_position, 'cd', 'markerfacecolor','c');
set(gca,'xscale', 'log');
% legend('Q10', 'Q20', 'Q30', 'Q40', 0);
box on;
rax = 15;
lax = 0.25;
axis([0.25 rax 0 2400]);
plot([0.25 rax],[layer2 layer2],'k:');
plot([0.25 rax],[layer3 layer3],'k:');
plot([0.25 rax],[layer4 layer4],'k:');
plot([0.25 rax],[layer5 layer5],'k:');
plot([0.25 rax],[layer6 layer6],'k:');
plot([0.25 rax],[whitematter whitematter],'k-');
set(gca,'ytick',[0 layer2 layer3 layer4 layer5 layer6 whitematter]);
xlabel('Q');
title('Tuning');
set(gca,'ydir','rev');




orient landscape;

print_mfilename(mfilename);





function mtf_paper_figure3_mtf_position_example(rtf_params)
%mtf_paper_figure3_mtf_position_example - plots best temporal, spectral modulation values
% along with asymmetry indices and separability indices.
%
% rtf_params : holds ripple transfer function parameters. 'rtf_params' is usually
% found in a file such as 2004-1-14-site2-2352um-30db-dmr1-fs18115-rtf-params.mat
%
% rtf_params is obtained from the function
%
%    rtf_params = strf_ripple_transfer_function(params)
%
% exportfig(gcf,'fig3_laminar_mtf_position_example1.eps','color','rgb','width', 3.6, 'height', 8.25, 'fontmode', 'fixed', 'fontsize', 8);
%
% Followed by:
% subplotspace('horizontal', -10);
% subplotspace('vertical', 5);
%
%   caa 10/28/02

pos = [rtf_params.position];
[temp, indpos] = sort(pos);
rtf_params = rtf_params(indpos);

pos = [rtf_params.position];
indpos = find(pos <= 2000);
rtf_params = rtf_params(indpos);

exp = rtf_params(1).exp;
site = rtf_params(1).site;
stim = rtf_params(1).stim;

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

% initialize the variables
position = zeros(1,length(rtf_params));
w0 = zeros(1,length(rtf_params));
asi = zeros(1,length(rtf_params));
left_sep_index = zeros(1,length(rtf_params));
right_sep_index = zeros(1,length(rtf_params));
best_tmf = zeros(1,length(rtf_params));
best_smf = zeros(1,length(rtf_params));

tmf = rtf_params(1).tmf;
xmf = rtf_params(1).xmf;

tmtf = zeros(length(rtf_params),length(tmf));
xmtf = zeros(length(rtf_params),length(xmf));

% First get all the separability indices

for i = 1:length(rtf_params)

   position(i) = rtf_params(i).position;
   w0(i) = rtf_params(i).w0;
   asi(i) = rtf_params(i).rtf_asymmetry_index(end);
   left_sep_index(i) = rtf_params(i).rtf_left_sep_index(end);
   right_sep_index(i) = rtf_params(i).rtf_right_sep_index(end);
%    best_tmf(i) = rtf_params(i).best_tmf(end);
%    best_xmf(i) = rtf_params(i).best_xmf(end);


   if ( iscell(rtf_params(i).tmtf) )
      tmtf_temp = rtf_params(i).tmtf{end};
      xmtf_temp = rtf_params(i).xmtf{end};
   else
      tmtf_temp = rtf_params(i).tmtf;
      xmtf_temp = rtf_params(i).xmtf;
   end

   tmtf(i,:) = tmtf_temp ./ max(tmtf_temp);
   xmtf(i,:) = xmtf_temp ./ max(xmtf_temp);

end



% Process the temporal modulation transfer functions
twidth6db = zeros(1,length(position));
twidth3db = zeros(1,length(position));
twc3db = zeros(1,length(position));
tbp3db = [];
tbp6db = [];
tlp3db = [];
tlp6db = [];

for i = 1:length(position)

   tfit = min(tmf):0.001:max(tmf);
   mtf = interp1(tmf, tmtf(i,:), tfit, 'spline');

%    plot(tmf, tmtf(i,:), 'ro-', tfit, mtf, 'g-');
%    axis([-0.5 40.5 0 1.1]);
%    title(sprintf('i = %.0f',i));

   indmax = find(mtf == max(mtf));

   ind3dblow = find(mtf(1:indmax) >= 0.7);
   ind3dbhi = find(mtf(indmax:end) >= 0.7);

   ind6dblow = find(mtf(1:indmax) >= 0.5);
   ind6dbhi = find(mtf(indmax:end) >= 0.5);
   
   twidth3db(i) = tfit(max(ind3dbhi+indmax-1)) - tfit(min(ind3dblow));
   twidth6db(i) = tfit(max(ind6dbhi+indmax-1)) - tfit(min(ind6dblow));
   twc3db(i) = tfit(max(ind3dbhi+indmax-1));

   if ( mtf(min(ind3dblow)) < 0.72 & mtf(max(ind3dbhi+indmax-1)) < 0.72 ) % it's bandpass
      itemp = find(mtf == max(mtf));
      best_tmf(i) = tfit(itemp);
      tbp3db = [tbp3db position(i)];
   else  % it's lowpass
      best_tmf(i) = tfit(max(ind3dbhi+indmax-1)) ./ 2;
      tlp3db = [tlp3db position(i)];
   end

   if ( mtf(min(ind6dblow)) < 0.52 & mtf(max(ind6dbhi+indmax-1)) < 0.52 )
      tbp6db = [tbp6db position(i)];
   else
      tlp6db = [tlp6db position(i)];
   end

end % (for i)



% Process the spectral modulation functions
xwidth6db = zeros(1,length(position));
xwidth3db = zeros(1,length(position));
xwc3db = zeros(1,length(position));
xbp3db = [];
xbp6db = [];
xlp3db = [];
xlp6db = [];

for i = 1:length(position)

   xfit = min(xmf):0.001:max(xmf);
   mtf = interp1(xmf, xmtf(i,:), xfit, 'spline');

%    plot(xmf, xmtf(i,:), 'ro-', xfit, mtf, 'g-');
%    axis([-0.05 4.05 0 1.1]);

   indmax = find(mtf == max(mtf));

   ind3dblow = find(mtf(1:indmax) >= 0.7);
   ind3dbhi = find(mtf(indmax:end) >= 0.7);

   ind6dblow = find(mtf(1:indmax) >= 0.5);
   ind6dbhi = find(mtf(indmax:end) >= 0.5);
   
   xwidth3db(i) = xfit(max(ind3dbhi+indmax-1)) - xfit(min(ind3dblow));
   xwidth6db(i) = xfit(max(ind6dbhi+indmax-1)) - xfit(min(ind6dblow));
   xwc3db(i) = xfit(max(ind3dbhi+indmax-1));

   if ( mtf(min(ind3dblow)) < 0.72 & mtf(max(ind3dbhi+indmax-1)) < 0.72 )  % it's bandpass
      istemp = find(mtf == max(mtf));
      best_smf(i) = xfit(istemp);
      xbp3db = [xbp3db position(i)];
   else  % it's lowpass
      best_smf(i) = xfit(max(ind3dbhi+indmax-1)) ./ 2;
      xlp3db = [xlp3db position(i)];
   end

   if ( mtf(min(ind6dblow)) < 0.52 & mtf(max(ind6dbhi+indmax-1)) < 0.52 )
      xbp6db = [xbp6db position(i)];
   else
      xlp6db = [xlp6db position(i)];
   end

end % (for i)




figure;

% Plot temporal modulation transfer function data

dt = diff(tmf);
dt = dt(1);
dtindex = ceil(1/dt);

subplot(4,2,1);
imagesc(tmtf);
upos = unique(position);
for i = 1:length(upos)
   index(i) = min(find(upos(i) == position));
end
xtick = [1 8*dtindex+1 8*2*dtindex+1 8*3*dtindex+1 8*4*dtindex+1 length(tmf)];
set(gca,'xtick', xtick, 'xticklabel', tmf(xtick));
set(gca,'ytick',index,'yticklabel',upos);
set(gca,'tickdir', 'out');
xlabel('tMTF (cyc/s)', 'fontname', 'arial', 'fontsize', 8);
ylabel('Position (um)', 'fontname', 'arial', 'fontsize', 8);
title('Temporal', 'fontname', 'arial', 'fontsize', 8);
set(gca,'fontname', 'arial', 'fontsize', 8);
colorbar;



% First plot the best temporal modulation frequency vs. depth
subplot(4,2,3);
hold on;
plot(best_tmf, position, 'ko', 'markerfacecolor','k');

posdata = min(min(position)):max(max(position));
fitmodel = fit(position(:), best_tmf(:), 'smoothingspline', 'SmoothingParam', 0.25e-5);
[fitdata] = feval(fitmodel, posdata);
plot(fitdata, posdata, 'k-', 'linewidth', 1);

box on;
% lax = -0.1 * max(twidth6db);
% rax = 1.1 * max(twidth6db);
lax = -2;
rax = 42;
axis([lax rax 0 2200]);
plot([lax rax],[layer2 layer2],'k:');
plot([lax rax],[layer3 layer3],'k:');
plot([lax rax],[layer4 layer4],'k:');
plot([lax rax],[layer5 layer5],'k:');
plot([lax rax],[layer6 layer6],'k:');
plot([lax rax],[whitematter whitematter],'k-');
tax = 35;
text([tax],[layer1label],'I', 'fontname', 'arial', 'fontsize', 8);
text([tax],[layer2label],'II', 'fontname', 'arial', 'fontsize', 8);
text([tax],[layer3label],'III', 'fontname', 'arial', 'fontsize', 8);
text([tax],[layer4label],'IV', 'fontname', 'arial', 'fontsize', 8);
text([tax],[layer5label],'V', 'fontname', 'arial', 'fontsize', 8);
text([tax],[layer6label],'VI', 'fontname', 'arial', 'fontsize', 8);
% text([tax],[whitematterlabel],'WM');

set(gca, 'xtick', [0:5:40], 'xticklabel', [0:5:40]);
set(gca,'ytick',[0 layer2 layer3 layer4 layer5 layer6 whitematter]);
set(gca,'ydir','rev');
set(gca,'tickdir', 'out');
xlabel('bTMF (cyc/s)', 'fontname', 'arial', 'fontsize', 8);
ylabel('Position (um)', 'fontname', 'arial', 'fontsize', 8);
set(gca,'fontname', 'arial', 'fontsize', 8);



% % width of temporal mtf vs. depth
% subplot(4,2,5);
% hold on;
% plot(twidth3db, position, 'ko', 'markerfacecolor','k');
% posdata = min(min(position)):max(max(position));
% fitmodel = fit(position(:), twidth3db(:), 'smoothingspline', 'SmoothingParam', 0.25e-5);
% [twdata] = feval(fitmodel, posdata);
% plot(twdata, posdata, 'k-', 'linewidth', 1);
% % plot(twidth6db, position, 'ro', 'markerfacecolor','r');
% % legend('3db');
% box on;
% % lax = -0.1 * max(twidth6db);
% % rax = 1.1 * max(twidth6db);
% lax = -2;
% rax = 42;
% axis([lax rax 0 2200]);
% plot([lax rax],[layer2 layer2],'k:');
% plot([lax rax],[layer3 layer3],'k:');
% plot([lax rax],[layer4 layer4],'k:');
% plot([lax rax],[layer5 layer5],'k:');
% plot([lax rax],[layer6 layer6],'k:');
% plot([lax rax],[whitematter whitematter],'k-');
% 
% set(gca, 'xtick', [0:5:40], 'xticklabel', [0:5:40]);
% set(gca,'ytick',[0 layer2 layer3 layer4 layer5 layer6 whitematter]);
% set(gca,'ydir','rev');
% set(gca,'tickdir', 'out');
% xlabel('tMTF Width (cyc/s)', 'fontname', 'arial', 'fontsize', 8);
% ylabel('Position (um)', 'fontname', 'arial', 'fontsize', 8);
% set(gca,'fontname', 'arial', 'fontsize', 8);



% subplot(4,2,7);
% tposbins = [0 150 450 750 1050 1350 1650 1950 2250];
% tnbp = histc([tbp3db 3000], tposbins);
% tnlp = histc([tlp3db 3000], tposbins);
% bincenters = [0 300 600 900 1200 1500 1800 2100 2400];
% bar(bincenters, [tnlp(:) tnbp(:)], 'stacked');
% axis([-100 2250 0 max(max([tnbp+tnlp]))+1]);
% set(gca,'xtick', bincenters, 'xticklabel', bincenters);
% set(gca,'tickdir', 'out');
% xlabel('Position (um)', 'fontname', 'arial', 'fontsize', 8);
% ylabel('# Units', 'fontname', 'arial', 'fontsize', 8);
% set(gca,'fontsize', 8);
% children = get(gca,'children');
% set(children(1),'facecolor', [0.75 0.75 0.75]);
% set(children(2),'facecolor', [0 0 0]);
% legend('LowPass', 'BandPass');
% set(gca,'fontname', 'arial', 'fontsize', 8);




% Plot spectral modulation transfer function data

dx = diff(xmf);
dx = dx(1);
dxindex = ceil(1/dx);

if ( strcmp(exp,'2004-1-14') )
   temp = smooth(xmtf(2,:));
   xmtf(2,:) = temp(:)';
end

subplot(4,2,2);
f = round(xmf * 10)/10;
imagesc(xmtf);
xtick = [1 dxindex+1 2*dxindex+1 3*dxindex+1 length(xmf)];
set(gca,'xtick', xtick,'xticklabel', [0 1 2 3 4]);
set(gca, 'ytick', index, 'yticklabel', '');
set(gca, 'tickdir', 'out');
xlabel('sMTF (cyc/oct)', 'fontname', 'arial', 'fontsize', 8);
title('Spectral', 'fontname', 'arial', 'fontsize', 8);
set(gca,'fontname', 'arial', 'fontsize', 8);


% plot best spectral modulation frequency vs. depth
subplot(4,2,4);
hold on;
plot(best_smf, position, 'ko', 'markerfacecolor','k');
posdata = min(min(position)):max(max(position));
fitmodel = fit(position(:), best_smf(:), 'smoothingspline', 'SmoothingParam', 0.25e-5);
[fitdata] = feval(fitmodel, posdata);
plot(fitdata, posdata, 'k-', 'linewidth', 1);

box on;
lax = -0.2;
rax = 4.2;
% lax = -0.1 * max(xwidth6db);
% rax = 1.1 * max(xwidth6db);
axis([lax rax 0 2200]);
plot([lax rax],[layer2 layer2],'k:');
plot([lax rax],[layer3 layer3],'k:');
plot([lax rax],[layer4 layer4],'k:');
plot([lax rax],[layer5 layer5],'k:');
plot([lax rax],[layer6 layer6],'k:');
plot([lax rax],[whitematter whitematter],'k-');

set(gca, 'xtick', [0:0.5:4], 'xticklabel', [0:0.5:4]);
set(gca,'ytick',[0 layer2 layer3 layer4 layer5 layer6 whitematter], 'yticklabel', '');
set(gca,'ydir','rev');
set(gca,'tickdir', 'out');
xlabel('bSMF (cyc/oct)', 'fontname', 'arial', 'fontsize', 8);
set(gca,'fontname', 'arial', 'fontsize', 8);




% subplot(4,2,6);
% hold on;
% plot(xwidth3db, position, 'ko', 'markerfacecolor','k');
% posdata = min(min(position)):max(max(position));
% fitmodel = fit(position(:), xwidth3db(:), 'smoothingspline', 'SmoothingParam', 0.1e-5);
% [xwdata] = feval(fitmodel, posdata);
% plot(xwdata, posdata, 'k-', 'linewidth', 1);
% % plot(xwidth6db, position, 'ro', 'markerfacecolor','r');
% box on;
% % lax = -0.1 * max(xwidth6db);
% % rax = 1.1 * max(xwidth6db);
% lax = -0.2;
% rax = 4.2;
% axis([lax rax 0 2200]);
% plot([lax rax],[layer2 layer2],'k:');
% plot([lax rax],[layer3 layer3],'k:');
% plot([lax rax],[layer4 layer4],'k:');
% plot([lax rax],[layer5 layer5],'k:');
% plot([lax rax],[layer6 layer6],'k:');
% plot([lax rax],[whitematter whitematter],'k-');
% 
% set(gca, 'xtick', [0:0.5:4], 'xticklabel', [0:0.5:4]);
% set(gca,'ytick',[0 layer2 layer3 layer4 layer5 layer6 whitematter], 'yticklabel', '');
% set(gca,'ydir','rev');
% set(gca,'tickdir', 'out');
% xlabel('sMTF Width (cyc/oct)', 'fontname', 'arial', 'fontsize', 8);
% set(gca,'fontname', 'arial', 'fontsize', 8);



% subplot(4,2,8);
% xposbins = [0 150 450 750 1050 1350 1650 1950 2250];
% xnbp = histc(xbp3db, xposbins);
% xnlp = histc(xlp3db, xposbins);
% bincenters = [0 300 600 900 1200 1500 1800 2100 2400];
% bar(bincenters, [xnlp(:) xnbp(:)], 'stacked');
% axis([-100 2250 0 max(max([xnbp+xnlp]))+1]);
% set(gca,'xtick', bincenters, 'xticklabel', bincenters);
% set(gca, 'tickdir', 'out');
% xlabel('Position (um)', 'fontname', 'arial', 'fontsize', 8);
% set(gca,'fontsize', 8)
% children = get(gca,'children');
% set(children(1),'facecolor', [0.75 0.75 0.75]);
% set(children(2),'facecolor', [0 0 0]);
% set(gca,'fontname', 'arial', 'fontsize', 8);

orient tall;

print_mfilename(mfilename);



