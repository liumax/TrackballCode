function [rtf_fold, tmf_mtf, tmtf, smf_mtf, smtf] = rtf2mtf(rtf, tmf, xmf)
%RTF2MTF - Modulation transfer functions from ripple transfer functions
%
% [rtf_fold, tmf, tmtf, smf, smtf] = rtf2mtf(rtf, tmf, xmf)
% ---------------------------------------------------------
%
% rtf : absolute value of 2D FFT of STRF
% tmf : temporal modulation frequency of rtf
% xmf : spectral modulation frequency of rtf
%
% rtf_fold : rtf folded about the tmf == 0 axis
% tmtf : temporal MTF
% tmf_mtf : temporal modulation frequency axis for rtf_fold
% smf_mtf : spectral modulation frequency axis for rtf_fold
% smtf : spectral MTF
%
% caa 12/2/09



% temporal and spectral modulation frequency axes
tmf = tmf(:)';
rtf_tmf = tmf;
ind0 = find(tmf==0);
tmf_mtf = tmf(ind0:end);

smf_mtf = xmf;
smf_mtf = smf_mtf(:)';

% the folded version of the ripple transfer function
rtf_right = rtf(:,ind0:end);
rtf_left = rtf(:,1:ind0-1);
rtf_fold = zeros(size(rtf_right));
rtf_fold = rtf_fold + rtf_right;
rtf_fold(:,2:end) = rtf_fold(:,2:end) + fliplr(rtf_left);
rtf_fold(:,1) = 2*rtf_fold(:,1);


% temporal modulation transfer function
tmtf = sum(rtf,1); % sum across spectral modulation frequency
tmtf = tmtf(:)';
tmtf_left = tmtf(1:ind0-1);
tmtf_right = tmtf(ind0:end);
tmtf = tmtf_right + [0 fliplr(tmtf_left)];
tmtf(1) = 2 * tmtf(1);
tmtf = tmtf ./ max(tmtf);


% spectral modulation transfer function
smtf = sum(rtf,2)'; % sum across temporal modulation frequency
smtf = smtf(:)';
smtf = smtf ./ max(smtf);

return;








