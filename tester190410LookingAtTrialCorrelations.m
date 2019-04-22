%This is a case study to look at different units and see what their
%responses look like
cd Z:\Max\maxKilo\DMRTuning\190206ML181105F_RAudStr3633pen2rec1tuningDMR
load('190206ML181105F_RAudStr3633pen2rec1tuningDMRAnalysis.mat')

%pairs include:
%FSI MSN
%265 270
%273 272, 277, 84, 40
%283 40, 15, 95


%lets start with 265 and 270. 

dataFSI = s.clu265;
dataMSN = s.clu270;

respFSI = dataFSI.LatPeakBin;
respMSN = dataMSN.LatPeakBin;

sigCutoff = 0.05;
figure
counter = 1;
subplot = @(m,n,p) subtightplot (m, n, p, [0.02 0.02], [0.05 0.05], [0.05 0.05]);
for i = 1:size(respFSI,1)
    for j = 1:size(respFSI,2)
        if i == 1
            subplot(size(respFSI,1),size(respFSI,2),j)
        else
            subplot(size(respFSI,1),size(respFSI,2),size(respFSI,2)*(i-1) + j)
        end
        hold on
        
        tarPull = respFSI{i,j};
        binDiffFSI = tarPull.BinnedSpikesTone - tarPull.BinnedSpikesToneBase;
        
        tarPull = respMSN{i,j};
        binDiffMSN = tarPull.BinnedSpikesTone - tarPull.BinnedSpikesToneBase;
        bigStore(counter:counter + length(binDiffFSI) - 1,:) = [binDiffFSI-mean(binDiffFSI),binDiffMSN-mean(binDiffMSN)];
        plot(binDiffFSI,binDiffMSN,'k.')
        [b,bintr,bintjm] = gmregress(binDiffFSI,binDiffMSN,sigCutoff);
        plot([min(binDiffFSI) max(binDiffFSI)],[min(binDiffFSI)*b(2) + b(1) max(binDiffFSI)*b(2) + b(1)],'r')
        title(num2str(b(2)))
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        slopeStore(size(respFSI,2)*(i-1)+j) = b(2);
        counter = counter + length(binDiffFSI);
    end
end


%try shuffling trial.
shuffVect = [2:1:20];
shuffVect(end+1) = 1;
sigCutoff = 0.05;
figure
counter = 1;
subplot = @(m,n,p) subtightplot (m, n, p, [0.02 0.02], [0.05 0.05], [0.05 0.05]);
for i = 1:size(respFSI,1)
    for j = 1:size(respFSI,2)
        if i == 1
            subplot(size(respFSI,1),size(respFSI,2),j)
        else
            subplot(size(respFSI,1),size(respFSI,2),size(respFSI,2)*(i-1) + j)
        end
        hold on
        
        tarPull = respFSI{i,j};
        binDiffFSI = tarPull.BinnedSpikesTone - tarPull.BinnedSpikesToneBase;
        
        tarPull = respMSN{i,j};
        binDiffMSN = tarPull.BinnedSpikesTone - tarPull.BinnedSpikesToneBase;
        bigStoreShuff(counter:counter + length(binDiffFSI) - 1,:) = [binDiffFSI(shuffVect)-mean(binDiffFSI),binDiffMSN-mean(binDiffMSN)];
        plot(binDiffFSI(shuffVect),binDiffMSN,'k.')
        [b,bintr,bintjm] = gmregress(binDiffFSI(shuffVect),binDiffMSN,sigCutoff);
        plot([min(binDiffFSI) max(binDiffFSI)],[min(binDiffFSI)*b(2) + b(1) max(binDiffFSI)*b(2) + b(1)],'r')
        title(num2str(b(2)))
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        slopeStore(size(respFSI,2)*(i-1)+j) = b(2);
        counter = counter + length(binDiffFSI);
    end
end
