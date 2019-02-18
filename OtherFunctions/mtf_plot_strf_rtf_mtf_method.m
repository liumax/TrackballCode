function varargout = mtf_plot_strf_rtf_mtf_method(strf,trigger)
% mtf_plot_strf_rtf_mtf_method Estimate/Plot RTF,MTF from STA
% 
%     mtf_plot_strf_rtf_mtf_method(strf,trigger)
% 
%     strf : struct element or struct array holding strf data. Only the
%     first element of strf is used, strf(1), for estimating mtf data
% 
%     trigger : vector of trigger times, in sample number.
% 
%     The sta in strf is thresholded at the p = 0.01 level, and then the
%     rtf is estimated from the significant sta.
% 
%     mtfs are derived from the rtf, and the best modulation frequencies
%     are calculated.
%
%     res = mtf_plot_strf_rtf_mtf_method(strf,trigger) returns the sta,
%     rtf, and mtf estimates.




%Define some initial parameters
exp = strf(1).exp;
site = strf(1).site;
chan = strf(1).chan;
model = strf(1).model;
position = strf(1).position;
pval = 0.01;
n0 = strf(1).n0contra;
w0 = strf(1).w0contra;
sm = strf(1).sm;
tm = strf(1).tm;
mdb = strf(1).mdb;
stim = strf(1).stim;
fs = strf(1).fs;
taxis = strf(1).taxis;
faxis = strf(1).faxis;
dur = (trigger(end)-trigger(1)) / fs; % duration in seconds


% Get the significant strf
[sta] = significant_strf(strf(1).rfcontra, pval, n0, mdb, dur);


% Compute RTF
[rtf_tmf, rtf_smf, rtf] = sta2rtf(sta, taxis, faxis, tm, sm);


% Fold RTF about TM = 0 axis
[rtf_fold, tmf_mtf, tmtf, smf_mtf, smtf] = rtf2mtf(rtf, rtf_tmf, rtf_smf);


% Compute Temporal/Spectral MTFs
[tbmf, tbw6db, tbw3db, twc3db] = mtf_bmf_bw_wc(tmf_mtf, tmtf);
[sbmf, sbw6db, sbw3db, swc3db] = mtf_bmf_bw_wc(smf_mtf, smtf);



% Plot the results: sta, rtf, and mtfs
figure;
colormap jet
subplot(2,2,1);
imagesc(taxis*1000, faxis/1000, sta);
axis xy;
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
xlabel('Time (ms)');
ylabel('Freq (kHz)');
ht = title(sprintf('%s site%.0f chan%.0f-model%.0f %.0f um', ...
 exp, site, chan, model, position),'Interpreter','none');


subplot(2,2,2);
imagesc(rtf_tmf, rtf_smf, rtf);
axis xy;
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
xlabel('Temp Mod Freq (cyc/sec)');
ylabel('Spec Mod Freq (cyc/oct)');
title(sprintf('mfile = %s', mfilename), 'Interpreter', 'none');


subplot(2,2,3);
hold on;
plot(tmf_mtf, tmtf, 'ko-', 'markerfacecolor', 'k', ...
'markersize', 2);
plot([tbmf tbmf], [0 1], 'k-');
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
ytick = 0:0.25:1;
set(gca,'ytick', ytick, 'yticklabel', ytick);
xlabel('Temporal Modulation Freq (cyc/sec)');
ylabel('Norm. Amp.');
title(sprintf('tMTF: %.2f Hz',tbmf));


subplot(2,2,4);
hold on;
plot(smf_mtf, smtf, 'ko-', 'markerfacecolor', 'k', ...
'markersize', 2);
plot([sbmf sbmf], [0 1], 'k-');
set(gca,'tickdir', 'out', 'ticklength', [0.025 0.025]);
ytick = 0:0.25:1;
set(gca,'ytick', ytick, 'yticklabel', ytick);
xlabel('Spectral Modulation Freq (cyc/oct)');
ylabel('Norm. Amp.');
title(sprintf('sMTF: %.2f cyc/oct',sbmf));

if ( nargout == 1 )
    res.exp = exp;
    res.site = site;
    res.chan = chan;
    res.model = model;
    res.position = position;
    res.n0 = n0;
    res.w0 = w0;
    res.sm = sm;
    res.tm = tm;
    res.mdb = mdb;
    res.stim = stim;
    res.fs = fs;
    res.taxis = taxis;
    res.faxis = faxis;
    res.sta = sta;
    res.tmfRtf = rtf_tmf;
    res.smfRtf = rtf_smf;
    res.rtf = rtf;
    res.tmfMtf = tmf_mtf;
    res.smfMtf = smf_mtf;
    res.rtfFold = rtf_fold;
    res.tmtf = tmtf;
    res.smtf = smtf;
    res.tbmf = tbmf;
    res.tbw6db = tbw6db;
    res.tbw3db = tbw3db;
    res.twc3db = twc3db;
    res.sbmf = sbmf;
    res.sbw6db = sbw6db;
    res.sbw3db = sbw3db;
    res.swc3db = swc3db;

    varargout{1} = res;
end

return;





