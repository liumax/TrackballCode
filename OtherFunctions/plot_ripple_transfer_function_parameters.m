function plot_ripple_transfer_function_parameters(rtf_params)
%plot_ripple_transfer_function_parameters - plots best temporal, spectral modulation values
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


temp = [rtf_params.position];
[temp, sort_index] = sort(temp);
rtf_params = rtf_params(sort_index);

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
   best_tmf(i) = rtf_params(i).best_tmf(end);
   best_xmf(i) = rtf_params(i).best_xmf(end);


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

   if ( mtf(min(ind3dblow)) < 0.72 & mtf(max(ind3dbhi+indmax-1)) < 0.72 )
      tbp3db = [tbp3db position(i)];
   else
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

   if ( mtf(min(ind3dblow)) < 0.72 & mtf(max(ind3dbhi+indmax-1)) < 0.72 )
      xbp3db = [xbp3db position(i)];
   else
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

% First plot the best temporal modulation frequency vs. depth
subplot(2,5,1);
hold on;
plot(best_tmf, position, 'ko', 'markerfacecolor','k');
% set(gca,'xscale', 'log');
box on;
lax = -5;
rax = 40;
axis([lax rax 0 2400]);
plot([lax rax],[layer2 layer2],'k:');
plot([lax rax],[layer3 layer3],'k:');
plot([lax rax],[layer4 layer4],'k:');
plot([lax rax],[layer5 layer5],'k:');
plot([lax rax],[layer6 layer6],'k:');
plot([lax rax],[whitematter whitematter],'k-');
tax = 30;
text([tax],[layer1label],'L I');
text([tax],[layer2label],'L II');
text([tax],[layer3label],'L III');
text([tax],[layer4label],'L IV');
text([tax],[layer5label],'L V');
text([tax],[layer6label],'L VI');
text([tax],[whitematterlabel],'WM');

set(gca,'ytick',[0 layer2 layer3 layer4 layer5 layer6 whitematter]);
xlabel('bTMF (Hz)');
ylabel('Depth (um)');
set(gca,'ydir','rev');


subplot(2,5,2);
hold on;
for i = 1:length(position)
   plot(tmf, 2400 - position(i) + 150*tmtf(i,:) - mean(150*tmtf(i,:)),'k-');
   p(i) = 2400 - position(i);
end
xlabel('TMTF (cyc/s)');
set(gca,'ytick',(unique(p)),'yticklabel',fliplr(unique(position)));
axis([0 40 0 2400]);
% set(gca,'ydir','rev');




dt = diff(tmf);
dt = dt(1);
dtindex = ceil(1/dt);

subplot(2,5,3);
imagesc(tmtf);
xlabel('TMTF (cyc/s)');
upos = unique(position);
for i = 1:length(upos)
   index(i) = min(find(upos(i) == position));
end
xtick = [1 8*dtindex+1 8*2*dtindex+1 8*3*dtindex+1 8*4*dtindex+1 length(tmf)];
set(gca,'xtick', xtick, 'xticklabel', tmf(xtick));
set(gca,'ytick',index,'yticklabel',upos);
set(gca,'fontsize', 8)



% width of temporal mtf vs. depth
subplot(2,5,4);
hold on;
plot(twidth6db, position, 'ko', 'markerfacecolor','k');
plot(twidth3db, position, 'ro', 'markerfacecolor','r');
legend('6db','3db');
box on;
lax = -0.1 * max(twidth6db);
rax = 1.1 * max(twidth6db);
axis([lax rax 0 2400]);
plot([lax rax],[layer2 layer2],'k:');
plot([lax rax],[layer3 layer3],'k:');
plot([lax rax],[layer4 layer4],'k:');
plot([lax rax],[layer5 layer5],'k:');
plot([lax rax],[layer6 layer6],'k:');
plot([lax rax],[whitematter whitematter],'k-');
set(gca,'ytick',[0 layer2 layer3 layer4 layer5 layer6 whitematter]);
xlabel('TMTF Width (cyc/s)');
set(gca,'ydir','rev');


subplot(2,5,5);
tposbins = 0:200:2400;
tnbp = histc([tbp3db 3000], tposbins);
tnlp = histc([tlp3db 3000], tposbins);
bar(tposbins, [tnlp(:) tnbp(:)], 'stacked');
axis([-100 2500 0 max(max([tnbp+tnlp]))+1]);
set(gca,'xtick', [0 600 1200 1800 2400], 'xticklabel', [0 600 1200 1800 2400]);
set(gca,'fontsize', 8);
children = get(gca,'children');
set(children(1),'facecolor', [0.75 0.75 0.75]);
set(children(2),'facecolor', [0 0 0]);
legend('LP3dB', 'BP3dB');


% plot best spectral modulation frequency vs. depth
subplot(2,5,6);
hold on;
plot(best_xmf, position, 'ko', 'markerfacecolor','k');
box on;
lax = -0.5;
rax = 4;
axis([lax rax 0 2400]);
plot([lax rax],[layer2 layer2],'k:');
plot([lax rax],[layer3 layer3],'k:');
plot([lax rax],[layer4 layer4],'k:');
plot([lax rax],[layer5 layer5],'k:');
plot([lax rax],[layer6 layer6],'k:');
plot([lax rax],[whitematter whitematter],'k-');
set(gca,'ytick',[0 layer2 layer3 layer4 layer5 layer6 whitematter]);
xlabel('bSMF (cyc/oct)');
ylabel('Depth (um)');
set(gca,'ydir','rev');


subplot(2,5,7);
hold on;
for i = 1:length(position)
%    plot(xmf, 150*xmtf(i,:)+position(i)-75,'k-');
   plot(xmf, 2400 - position(i) + 150*xmtf(i,:) - mean(150*xmtf(i,:)),'k-');
   p(i) = 2400 - position(i);
end
xlabel('XMTF (cyc/oct)');
set(gca,'ytick',(unique(p)),'yticklabel',fliplr(unique(position)));
axis([0 4 0 2400]);
% set(gca,'ydir','rev');



dx = diff(xmf);
dx = dx(1);
dxindex = ceil(1/dx);

subplot(2,5,8);
f = round(xmf * 10)/10;
imagesc(xmtf);
xlabel('XMTF (cyc/oct)');
xtick = [1 dxindex+1 2*dxindex+1 3*dxindex+1 length(xmf)];
set(gca,'xtick', xtick,'xticklabel', [0 1 2 3 4]); %f(xtick) );
set(gca, 'ytick', index, 'yticklabel', upos);
set(gca, 'fontsize', 8)


subplot(2,5,9);
hold on;
plot(xwidth6db, position, 'ko', 'markerfacecolor','k');
plot(xwidth3db, position, 'ro', 'markerfacecolor','r');
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


subplot(2,5,10);
xposbins = 0:200:2400;
xnbp = histc(xbp3db, xposbins);
xnlp = histc(xlp3db, xposbins);
bar(xposbins, [xnlp(:) xnbp(:)], 'stacked');
axis([-100 2500 0 max(max([xnbp+xnlp]))+1]);
xlabel('Depth (um)');
set(gca,'xtick', [0 600 1200 1800 2400], 'xticklabel', [0 600 1200 1800 2400]);
set(gca,'fontsize', 8)
children = get(gca,'children');
set(children(1),'facecolor', [0.75 0.75 0.75]);
set(children(2),'facecolor', [0 0 0]);

orient landscape;

print_mfilename(mfilename);




% % firing rate vs. depth
% subplot(2,6,6);
% hold on;
% plot(w0, position, 'ko', 'markerfacecolor','k');
% box on;
% lax = -0.1*max(w0);
% rax = 1.1 * max(w0);
% axis([lax rax 0 2400]);
% plot([lax rax],[layer2 layer2],'k:');
% plot([lax rax],[layer3 layer3],'k:');
% plot([lax rax],[layer4 layer4],'k:');
% plot([lax rax],[layer5 layer5],'k:');
% plot([lax rax],[layer6 layer6],'k:');
% plot([lax rax],[whitematter whitematter],'k-');
% 
% set(gca,'ytick',[0 layer2 layer3 layer4 layer5 layer6 whitematter]);
% xlabel('Firing Rate (sp/s)');
% set(gca,'ydir','rev');
% 
% 
% % rtf asymmetry index vs. depth
% subplot(2,6,12);
% hold on;
% plot(asi, position, 'ko', 'markerfacecolor','k');
% % set(gca,'xscale', 'log');
% box on;
% lax = -1.1;
% rax = 1.1;
% axis([lax rax 0 2400]);
% plot([lax rax],[layer2 layer2],'k:');
% plot([lax rax],[layer3 layer3],'k:');
% plot([lax rax],[layer4 layer4],'k:');
% plot([lax rax],[layer5 layer5],'k:');
% plot([lax rax],[layer6 layer6],'k:');
% plot([lax rax],[whitematter whitematter],'k-');
% 
% set(gca,'ytick',[0 layer2 layer3 layer4 layer5 layer6 whitematter]);
% xlabel('RTF Asymmetry Index');
% set(gca,'ydir','rev');
% 
% orient landscape





% Here's some code that I previously used but discontinued. I'm putting it 
% here so I don't have to reinvent the wheel again if I need it.
% 
% 
% 
% mtf_bandpass = inline('k(1) ./ (1 + (f ./ k(2)).^(2.*k(3)) ) ./ (1 + (k(4) ./ f).^(2.*k(5)) ) + k(6)', 'k', 'f');
% 
% % amplitude constants
% k1bp = 0.5;
% k1lbbp = 0;
% k1ubbp = 1;
% 
% % lowpass frequencies
% k2bp = 16;
% k2lbbp = 0;
% k2ubbp = 40;
% 
% % lowpass order
% k3bp = 2;
% k3lbbp = 0;
% k3ubbp = 6;
% 
% % high pass frequencies
% k4bp = 4;
% k4lbbp = 0;
% k4ubbp = 40;
% 
% % high pass order
% k5bp = 2;
% k5lbbp = 0;
% k5ubbp = 6;
% 
% % dc offset
% k6bp = 0.75;
% k6lbbp = 0.25;
% k6ubbp = 1;
% 
% 
% mtf_lowpass = inline('k(1) ./ (1 + (f ./ k(2)).^(2.*k(3)) ) + k(4)', 'k', 'f');
% 
% % amplitude constants
% k1lp = 1;
% k1lblp = 0;
% k1ublp = inf;
% 
% % lowpass frequencies
% k2lp = 1;
% k2lblp = 0;
% k2ublp = 40;
% 
% % lowpass order
% k3lp = 2;
% k3lblp = 0;
% k3ublp = 6;
% 
% % dc offset
% k4lp = 0.5;
% k4lblp = 0;
% k4ublp = 1;
% 
% bp = zeros(1,length(position));
% bppos = [];
% lp = zeros(1,length(position));
% lppos = [];
% 
% %    warning off MATLAB:divideByZero;
% %    [klp, r2_tmtf_lowpass] = ...
% %       lsqcurvefit(mtf_lowpass, [k1lp k2lp k3lp k4lp], tmf, tmtf(i,:), ...
% %       [k1lblp k2lblp k3lblp k4lblp], ...
% %       [k1ublp k2ublp k3ublp k4ublp]);
% % 
% %    tmtf_lowpass_fit = mtf_lowpass(klp, tmf);
% %    lowpass_fit = mtf_lowpass(klp, tfit);
% % 
% % 
% %    warning off MATLAB:divideByZero;
% %    [kbp, r2_tmtf_bandpass] = ...
% %       lsqcurvefit(mtf_bandpass, [k1bp k2bp k3bp k4bp k5bp k6bp], tmf, tmtf(i,:), ...
% %       [k1lbbp k2lbbp k3lbbp k4lbbp k5lbbp k6lbbp], ...
% %       [k1ubbp k2ubbp k3ubbp k4ubbp k5ubbp k6ubbp]);
% % 
% %    tmtf_bandpass_fit = mtf_bandpass(kbp, tmf);
% %    bandpass_fit = mtf_bandpass(kbp, tfit);
% % 
% % 
% %    sselp = (tmtf(i,:)-tmtf_lowpass_fit) * (tmtf(i,:)-tmtf_lowpass_fit)';
% %    ssebp = (tmtf(i,:)-tmtf_bandpass_fit) * (tmtf(i,:)-tmtf_bandpass_fit)';
% 
% %    if ( sselp < ssebp )
% %       lp(i) = 1;
% %       lppos = [lppos position(i)];
% %       width(i) = klp(2);
% %    else
% %       bp(i) = 1;
% %       bppos = [bppos position(i)];
% %       width(i) = kbp(2) - kbp(4);
% %    end
% 
% % %    plot(tmf, tmtf(i,:), 'ko-', tfit, mtf, 'r-', tfit, smooth(mtf,3), 'b-');
% %    plot(tmf, tmtf(i,:), 'ro-', tfit, mtf, 'g-');
% % %    plot(tmf, tmtf(i,:), 'ko-', tfit, lowpass_fit, 'r-', tfit, bandpass_fit, 'b-', tfit, mtf, 'g-');
% % %    subplot(2,1,2);
% % %    plot(tmf, tmtf(i,:), 'ko-', tmf, smooth(tmtf(i,:),2), 'bo-');
% %    title(sprintf('sse lp = %.5f, sse bp = %.5f', sselp, ssebp));
% %    axis([-0.5 40.5 0 1.1]);
% %    pause
% 
% 
% mtf_bandpass = inline('k(1) ./ (1 + (f ./ k(2)).^(2.*k(3)) ) ./ (1 + (k(4) ./ f).^(2.*k(5)) ) + k(6)', 'k', 'f');
% % amplitude constants
% k1bp = 0.5;k1lbbp = 0;k1ubbp = 1;
% % lowpass frequencies
% k2bp = 3;k2lbbp = 0;k2ubbp = 4;
% % lowpass order
% k3bp = 2;k3lbbp = 0;k3ubbp = 6;
% % high pass frequencies
% k4bp = 1;k4lbbp = 0;k4ubbp = 4;
% % high pass order
% k5bp = 2;k5lbbp = 0;k5ubbp = 6;
% % dc offset
% k6bp = 0.75;k6lbbp = 0.25;k6ubbp = 1;
% 
% 
% mtf_lowpass = inline('k(1) ./ (1 + (f ./ k(2)).^(2.*k(3)) ) + k(4)', 'k', 'f');
% % amplitude constants
% k1lp = 1; k1lblp = 0; k1ublp = inf;
% % lowpass frequencies
% k2lp = 1; k2lblp = 0; k2ublp = 4;
% % lowpass order
% k3lp = 2; k3lblp = 0; k3ublp = 6;
% % dc offset
% k4lp = 0.5; k4lblp = 0; k4ublp = 1;
% 
% 
% %    warning off MATLAB:divideByZero;
% %    [klp, r2_xmtf_lowpass] = ...
% %       lsqcurvefit(mtf_lowpass, [k1lp k2lp k3lp k4lp], xmf, xmtf(i,:), ...
% %       [k1lblp k2lblp k3lblp k4lblp], ...
% %       [k1ublp k2ublp k3ublp k4ublp]);
% %    xmtf_lowpass_fit = mtf_lowpass(klp, xmf);
% % 
% % 
% %    warning off MATLAB:divideByZero;
% %    [kbp, r2_xmtf_bandpass] = ...
% %       lsqcurvefit(mtf_bandpass, [k1bp k2bp k3bp k4bp k5bp k6bp], xmf, xmtf(i,:), ...
% %       [k1lbbp k2lbbp k3lbbp k4lbbp k5lbbp k6lbbp], ...
% %       [k1ubbp k2ubbp k3ubbp k4ubbp k5ubbp k6ubbp]);
% %    xmtf_bandpass_fit = mtf_bandpass(kbp, xmf);
% 
% %    sselp = (xmtf(i,:)-xmtf_lowpass_fit) * (xmtf(i,:)-xmtf_lowpass_fit)';
% %    ssebp = (xmtf(i,:)-xmtf_bandpass_fit) * (xmtf(i,:)-xmtf_bandpass_fit)';
% % 
% %    if ( sselp < ssebp )
% %       lppos = [lppos position(i)];
% %       width(i) = klp(2);
% %    else
% %       bppos = [bppos position(i)];
% %       width(i) = kbp(2) - kbp(4);
% %    end
% % 
% %    plot(xmf, xmtf(i,:), 'ko-', xmf, xmtf_lowpass_fit, 'ro-', xmf, xmtf_bandpass_fit, 'bo-');
% %    title(sprintf('sse lp = %.5f, sse bp = %.5f', sselp, ssebp));
% %    axis([-0.05 4.05 0 1.1]);
% %    pause
% 

