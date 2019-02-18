function batch_plot_sta_mid_rtf(filtstr, filtopt)

if filtopt == 1
    gaussian = fspecial('gaussian', [2 2], 1);
end

mapscheme = 'rdbu';


for i = 1:length(filtstr)
    
    [tmf, xmf, sta_rtf, mid1_rtf, mid2_rtf] = calc_sta_mid_rtf(filtstr(i), filtopt);

    figure;

    time = 1000 .* filtstr(i).time; % change to ms
    freq = filtstr(i).freq / 1000;

    ytick = [1 floor(length(freq)/2) length(freq)];
    yticklabel = round(freq(ytick)*10)/10;

    xtick = [1 floor(length(time)/2)+1 length(time)];
    xticklabel = fliplr(round( time(xtick) ));

    site = filtstr(i).site;
    exp = filtstr(i).exp;
    unit = filtstr(i).unit;

    sta = filtstr(i).v_sta;
    v1 = filtstr(i).v1;
    v2 = filtstr(i).v2;

    if filtopt == 1
        sta = imfilter(sta,gaussian);
        v1 = imfilter(v1,gaussian);
        v2 = imfilter(v2,gaussian);
    end


    subplot(2,3,1);
    imagesc( sta );
    axis xy;
    minmin = min(min(sta));
    maxmax = max(max(sta));
    boundary = max([abs(minmin) abs(maxmax)]);
    set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
    set(gca,'xtick', xtick, 'xticklabel', xticklabel);
    set(gca,'ytick', ytick, 'yticklabel', yticklabel);
    set(gca, 'tickdir', 'out', 'ticklength', [0.025 0.025]);
    cmap = cschemes(mapscheme,15);
    colormap(gca, cmap);
    title('STA');
    ylabel('Frequency (kHz)');



    subplot(2,3,2);
    imagesc( v1 );
    axis xy;
    minmin = min(min(v1));
    maxmax = max(max(v1));
    boundary = max([abs(minmin) abs(maxmax)]);
    set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
    set(gca,'xtick', xtick, 'xticklabel', xticklabel);
    set(gca,'ytick', ytick, 'yticklabel', yticklabel);
    set(gca, 'tickdir', 'out', 'ticklength', [0.025 0.025]);
    set(gca, 'yticklabel', []);
    colormap(gca, cmap);
    title('MID1');
    xlabel('Time before spike (ms)');

    subplot(2,3,3);
    imagesc( v2 );
    axis xy;
    minmin = min(min(v2));
    maxmax = max(max(v2));
    boundary = max([abs(minmin) abs(maxmax)]);
    set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);
    set(gca,'xtick', xtick, 'xticklabel', xticklabel);
    set(gca,'ytick', ytick, 'yticklabel', yticklabel);
    set(gca, 'tickdir', 'out', 'ticklength', [0.025 0.025]);
    set(gca, 'yticklabel', []);
    cmap = cschemes(mapscheme,15);
    colormap(gca, cmap);
    title('MID2');

    subplot(2,3,4);
    imagesc(tmf,xmf,sta_rtf);
    axis xy
    set(gca, 'tickdir', 'out', 'ticklength', [0.025 0.025]);
    colormap(gca, 'jet');
    ylabel('spectral modulation (cyc/oct)')
    title('STA RTF')

    subplot(2,3,5);
    imagesc(tmf,xmf,mid1_rtf);
    axis xy
    set(gca, 'tickdir', 'out', 'ticklength', [0.025 0.025]);
    colormap(gca, 'jet');
    xlabel('temporal modulation (Hz)')
    title('MID1 RTF')

    subplot(2,3,6);
    imagesc(tmf,xmf,mid2_rtf);
    axis xy
    set(gca, 'tickdir', 'out', 'ticklength', [0.025 0.025]);
    colormap(gca, 'jet');
    title('MID2 RTF')
   

    suptitle(sprintf('%s-site%s-unit%s', exp, site, unit));
    set(gcf, 'Position', [100 100 900 600]);

    export_fig(sprintf('temp%d.pdf',i), '-transparent', '-nocrop')
    close all
end
    