

function medWFs = extractMedianWFs(clu, st, Fs, datPath, dataType, dataSize, chanMap, gain)
% medWFs = extractMedianWFs(clu, st, cids, cgs, Fs, datPath, dataType, dataSize, chanMap, gain)
%
% - medWFs [nClusters, nChannels, nSamples] will be median waveforms of 
%    every cluster represented in clu
% -- note: if you desire only "good", try this:
% >> cluSub = clu(ismember(clu, cids(cgs==2))); % and same for spike times
%
% - clu [nSpikes,1] cluster identities %MAX spClu
% - st [nSpikes,1] spike times in trodes sampletime 190327 edit %MAX spTime
% - Fs [1,1] sampling frequency %MAX 30000
% - datPath [string] filename %MAX (for phy.dat) E:\maxKilo\maxFsiXcorr5
% - dataType [string] e.g. 'int16' %Max 'int16'
% - dataSize [2,1] e.g. [nChannelsInFile nSamples] %Max 64 channels, 
% - chanMap [vector] note: assume this is zero-indexed
% - gain [1,1] to convert to uV if desired, otherwise use 1

nWFsToLoad = 1000; % I find this to be enough, but 100 not enough. Could try other values.

% window is -0.5 to 1.25ms
wfWin = -round(0.5/1000*Fs):round(1.25/1000*Fs); 
nWFsamps = numel(wfWin);

nChInFile = dataSize(1);
nSamp = dataSize(2);
mmf = memmapfile(datPath, 'Format', {dataType, [nChInFile nSamp], 'x'});

cids = unique(clu);

nClu = length(cids);
nCh = length(chanMap);

medWFs = zeros(nClu, nCh, nWFsamps);

for q = 1:nClu
    fprintf(1, 'cluster %d (%d/%d)\n', cids(q), q, nClu);
    theseST = st(clu==cids(q));
    %now lets add a limiter to avoid issues of exceeding bounds
    theseST((theseST + wfWin(end) > length(mmf.Data.x))) = [];
    theseST((theseST + wfWin(1) < 1)) = [];
    
    nWFsToLoad = min(nWFsToLoad, length(theseST));
    extractST = round(theseST(randperm(length(theseST), nWFsToLoad)));
    
    theseWF = zeros(nWFsToLoad, nCh, nWFsamps);
    for n=1:nWFsToLoad
%         disp(strcat('q:',num2str(q),'n:',num2str(n)))
%         disp(strcat('Value:',num2str(extractST(n)),',Size:',num2str(length(mmf.Data.x))))
        tempWF = mmf.Data.x(1:nChInFile,extractST(n)+wfWin(1):extractST(n)+wfWin(end));
        theseWF(n,:,:) = tempWF(chanMap,:);
    end
    
    
    medWF = squeeze(median(double(theseWF),1));
    medWFuV = medWF.*gain;
    
    medWFs(q,:,:) = medWFuV;
    
end