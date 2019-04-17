load('nonlinearity_demo.mat');
figure('Position',[363 314 1098 436]);
cmap = cschemes('rdbu', 100);
colormap(cmap);

x0 = filtstr.x0;
tempstim = stim_mat(x0:x0+25-1,:);

sta = filtstr.v_sta;
mid1 = filtstr.v1;
mid2 = filtstr.v2;

[xprior, xposterior] = ne_sta_stimulus_projection(sta, spktrain, tempstim);
[px,pspk,pxspk,xbincenter] = calc_px_pspk_pxspk(xprior,xposterior, 15);
filtstr.pspk = pspk;

xbins = xbincenter;
nl = pspk .* pxspk ./ px;

subplot(121)
quick_plot_sta(filtstr.v_sta, 'faxis', faxis);

subplot(122); hold on
maxmax = max(nl);
plot(xbins, nl, 'ko-', 'markerfacecolor', 'k');
minx = min(xbins);
maxx = max(xbins);
plot([minx-1 maxx+1], [pspk pspk], 'k--');
xlim([minx-2 maxx+2]);
ylim([0 maxmax]);
xlabel('Projection (SD)');
ylabel('P(spk|x)');




