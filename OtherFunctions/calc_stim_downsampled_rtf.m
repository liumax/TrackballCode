function [RDds, FMds] = calc_stim_downsampled_rtf(paramfile, trigger, edges, DF)

% paramfile = paramfile of stimulus
% binsize = trigger of one such recording
% edges = bin size edges

% Written by JS, 12/11/17

load(paramfile, 'RD', 'FM', 'Fs')

N = 32;
idx = 1:N:length(trigger)*N + 1;

RDtot = zeros(N * 1000, length(trigger));
FMtot = zeros(N * 1000, length(trigger));


for i = 1:length(trigger)
    
    tempidx = idx(i:i+1);
    
    RDi = RD(tempidx(1):tempidx(2));
    RDint = interp10(RDi,3);
    RDtot(:,i) = RDint(1:end-1);
    FMi = FM(tempidx(1):tempidx(2));
    FMint = interp10(FMi,3);
    FMtot(:,i) = FMint(1:end-1);
        
end

RDtot = RDtot(:)';
FMtot = FMtot(:)';

selectidx = round(linspace(1,length(RDtot), length(edges)-1));

RD = RDtot(selectidx);
FM = FMtot(selectidx);


RDds = RD(1:DF:end);
FMds = FM(1:DF:end);