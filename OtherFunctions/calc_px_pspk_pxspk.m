function [px,pspk,pxspk,varargout] = calc_px_pspk_pxspk(xprior,xposterior, numbins)

% ntrials = length(xposterior);

if ~exist('numbins','var')
    numbins = 15;
end

xmn = mean(xprior);
xstd = std(xprior);

xprior_scaled = (xprior - xmn) ./ xstd;
xbins_edges = linspace(min(xprior_scaled), max(xprior_scaled), numbins + 1);

xposterior_scaled = (xposterior - xmn) ./ xstd;

nx = histcounts(xprior_scaled, xbins_edges);
px = nx ./ sum(nx); % p(x)
px = px(:);

nxspk = histcounts(xposterior_scaled, xbins_edges);
pxspk = nxspk ./ sum( nxspk ); % p(x|spk)
pxspk = pxspk(:);

pspk = length(xposterior)/length(xprior);

varargout{1} = xbins_edges;
varargout{2} = xprior_scaled;


end

