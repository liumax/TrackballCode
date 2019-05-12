function [tmf, xmf, rtf] = sta2rtf(v, taxis, faxis, MaxFM, MaxRD, Display)
% STA2RTF - Ripple transfer function from an STA, or other filter
%
% [tmf, xmf, rtf] = sta2rtf(v, taxis, faxis, MaxFM, MaxRD, Display)
% -----------------------------------------------------------------
%
% v : filter. May be either STA, MID1, or MID2. Must be in correct shape (nf x nlags). 
% taxis   : Time Axis
% faxis   : Frequency Axis
% MaxFm   : Maximum Modulation Rate for Experiment
% MaxRD   : Maximum Ripple Density for Experiment
% Display : Display : 'y' or 'n'
%           Default : 'n'
%
% tmf : temporal modulation frequency axis (cycles / s).
% xmf : spectral modulaiton frequency axis (cycles / octave).
% rtf : ripple transfer function. A matrix where spectral modulation
%    varies along the rows and temporal modulation varies along the
%    columns of the matrix.
%
% caa

if ( nargin < 6 )
   Display = 'n';
end

FsT = 1/(taxis(3)-taxis(2));
FsX = 1/log2(faxis(2)/faxis(1));


dtf = MaxFM / 50; % temporal modlation! 190512 adjusted Max for data
ntbins = ceil(FsT / dtf); % sampling rate over temp res
ntbins = ntbins + ~rem(ntbins,2); % add a bin if we have an even number

dff = MaxRD / 100; % spectral modlation! 190512 adjusted Max for data
nfbins = ceil(FsX / dff); % sampling rate over freq resolution
nfbins = nfbins + ~rem(nfbins,2); % add a bin to get an odd number

%Computing Ripple Transfer Function
[nf, nt] = size(v);
% N2 = 2^nextpow2(nt)*2+1;
% N1 = 2^nextpow2(nf)*2+1;
% 
% N2 = 2^nextpow2(nt)*2+1;
% N1 = 2^nextpow2(nf)+1;

rtf = fft2(v,nfbins,ntbins);
rtf = fftshift( real( rtf.*conj(rtf) ) );


% Get the tm frequency vector - will be in cycles per second
if ( mod(size(rtf,2),2) )
   tmf = ( -(size(rtf,2)-1)/2:(size(rtf,2)-1)/2 ) / size(rtf,2) * FsT;
else
   tmf = ( -size(rtf,2)/2:(size(rtf,2)/2-1) ) / size(rtf,2) * FsT;
end


% Get the sm frequency vector - will be in cycles per octave
if ( mod(size(rtf,1),2) )
   xmf = [-(size(rtf,1)-1)/2:(size(rtf,1)-1)/2]/size(rtf,1) * FsX;
else
   xmf = [-size(rtf,1)/2:(size(rtf,1)/2-1)]/size(rtf,1) * FsX;
end


%Discarding Unecessary Data Samples
index1 = find( xmf >= 0 & xmf <= MaxRD );
index2 = find(tmf <= MaxFM & tmf >= -MaxFM);
rtf = rtf(index1,index2);
xmf = xmf(index1);
tmf = tmf(index2);


if ( strcmp(Display,'y') )
   figure;
   imagesc(tmf, xmf, rtf )
   axis xy;
   set(gca,'tickdir', 'out');
end


return;

