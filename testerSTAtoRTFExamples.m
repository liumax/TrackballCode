


masterFolder = pwd;
masterDir = dir;
masterDir = {masterDir.name};
%this code here is to remove folders reported in dir like "." or ".."
masterIndex = strfind(masterDir,'DMRData');
masterIndex = find(not(cellfun('isempty', masterIndex)));


%masterFolders is a listing of all folders that I want to target for batch
%analysis
masterFolders = masterDir(masterIndex);

[findString] = functionCellStringFind(masterFolders,'.pdf');
masterFolders(findString) = [];
[findString] = functionCellStringFind(masterFolders,'.fig');
masterFolders(findString) = [];


numFolders = length(masterFolders);

for i = 1:numFolders
    load(masterFolders{i})
    fileName = masterFolders{i};
    fileName = fileName(1:end-11);
    stepSize = 599.8333/length(stimulus);
    for i = 1:size(sta,1)
        [tmf, xmf, rtf] = sta2rtf(reshape(sta(i,:),100,[]), [stepSize:stepSize:stepSize*100], faxis, 40, 4, 'n');
        
        hFig = figure;
        set(hFig, 'Position', [10 10 1000 1000])
        subplot = @(m,n,p) subtightplot (m, n, p, [0.07 0.07], [0.05 0.05], [0.05 0.05]);
        subplot(2,2,1)
        imagesc([-0.1:0.001:-0.001],faxis,reshape(sta(i,:),100,[]))
        colormap('parula')
        xlabel('Time (s)')
        ylabel('Freq (Hz)')
        set(gca,'YDir','normal')
        title(fileName,'Interpreter','none')
        
        subplot(2,2,3)
        imagesc([-0.1:0.001:-0.001],faxis,reshape(sta_sig(i,:),100,[]))
        colormap('parula')
        xlabel('Time (s)')
        ylabel('Freq (Hz)')
        set(gca,'YDir','normal')
        title(num2str(sum(spikeArray(i,:))))
        
        subplot(2,2,2)
        imagesc(tmf, xmf, rtf )
        colormap('parula')
        axis xy
        xlabel('Temporal Modulation')
        ylabel('Spectral Modulation')
        title('RTF')
        
        subplot(2,2,4)
        temp1 = flip(rtf(:,1:16),2);
        temp2 = rtf(:,16:end);
        flipRTF = temp1 + temp2;
        
        imagesc(tmf(16:end), xmf, flipRTF )
        colormap('parula')
        axis xy
        xlabel('Temporal Modulation')
        ylabel('Spectral Modulation')
        title('RTF (folded)')
        
        spikeGraphName = strcat(fileName,num2str(i),'RTF');
        savefig(hFig,spikeGraphName);

        %save as PDF with correct name
        set(hFig,'Units','Inches');
        pos = get(hFig,'Position');
        set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        print(hFig,spikeGraphName,'-dpdf','-r0')
    end
    save(fullfile(pwd,strcat(fileName,'STA_RTF')),'sta','faxis','tmf', 'xmf', 'rtf')
%     close all
end


