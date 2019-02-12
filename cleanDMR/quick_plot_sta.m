function quick_plot_sta(sta, faxis, nlags)

if ~exist('nf','var')
    load('I:\Ripple_Noise\downsampled_for_MID\rn1-500flo-40000fhi-4SM-40TM-40db-96khz-48DF-10min_DFt1_DFf5_param.mat', 'faxis')
end
if ~exist('nlags','var')
    nlags = 20;
end

nf = length(faxis);

assert(numel(sta) == nf*nlags)

if size(sta,1) ~= nf
    sta = reshape(sta, nf, nlags);
end

ytick = 20:20:length(faxis);
ylab = round(faxis(20:20:length(faxis))/1000);

imagesc(sta)

xlabel('Time before spike (ms)')
ylabel('Frequency (kHz)')

set(gca, 'xtick',0:5:20, 'xticklabel',fliplr(0:25:100))
set(gca,'ydir', 'normal');
set(gca, 'ytick',ytick,'yticklabel',ylab)
set(gca, 'tickdir', 'out', 'ticklength', [0.025 0.025]);

boundary = max(abs(sta(:)));
set(gca, 'clim', [-1.05*boundary-eps 1.05*boundary+eps]);

tickpref;