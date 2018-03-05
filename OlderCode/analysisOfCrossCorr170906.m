
%actually, I'm not sure I see the point in doing an ANOVA. How about we try
%and do spike triggered averages instead? (cross correlograms). This way we
%can see if we can replicate the berke results. This will require a bit
%more of a deep dive into datasets. 

cd 170906TestingCrossCorr

tester = what;
tester = tester.mat;

numFiles = length(tester);
bigInd = 1

for i = 1:numFiles
    %open the targeted file
    load(tester{i})
    
    
    %check for any PV cells.
    pvCheck = find(s.MasterSheet(:,1) == 1)
    %if pv cells, do cross correlograms
    if length(pvCheck) > 0
        %find the number of cells
        numCells = length(s.MasterSheet(:,1));
        corrStore = cell(numCells,numCells);
        for j = 1: numCells
            j
            for k = 1:numCells
                cell1 = s.(s.DesignationName{j}).SpikeTimes;
                cell2 = s.(s.DesignationName{k}).SpikeTimes;
                newInd = 1;
                for fIND = 1:length(cell1)
                    subVal = cell2-cell1(fIND);
                    subVal(subVal>0.050) = [];
                    subVal(subVal<-0.050) = [];
                    if length(subVal) > 0
                        tempStore(newInd:newInd+length(subVal)-1) = subVal;
                        newInd = newInd + length(subVal);
                    end
                end
                corrStore{j,k} = tempStore;
                tempStore = [];
                indStore(j,k) = newInd;
                %also store information about what kinds of cells
                cellStore(j,k,1) = s.MasterSheet(j,1) ;
                cellStore(j,k,2) = + s.MasterSheet(k,1);
                %store whether same or different channel
                if s.DesignationArray(k,1) == s.DesignationArray(j,1)
                    sameStore(j,k) = 1;
                else
                    sameStore(j,k) = 0;
                end
            end
        end
        C = combnk([1:numCells],2);
        for combInd = 1:length(C)
            bigStore{bigInd} = corrStore{C(combInd,1),C(combInd,2)};
            infoStore(bigInd,1) = cellStore(C(combInd,1),C(combInd,2),1);
            infoStore(bigInd,2) = cellStore(C(combInd,1),C(combInd,2),2);
            infoStore(bigInd,3) = sameStore(C(combInd,1),C(combInd,2));
            infoStore(bigInd,4) = C(combInd,1);
            infoStore(bigInd,5) = C(combInd,2);
            bigInd = bigInd + 1;
        end
        
        %now plot out histogram
        subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.01], [0.03 0.05], [0.03 0.01]);

        hFig = figure
        set(hFig, 'Position', [10 80 1240 850])
        histVect = [-0.05:0.001:0.05];
        for j = 1:numCells
            for k = 1:numCells
                subplot(numCells,numCells,(j-1)*numCells + k)
                tempHold = hist(corrStore{j,k},histVect);
                if s.MasterSheet(j,1) == 1 & s.MasterSheet(k,1) == 1 %PV vs PV
                    if s.DesignationArray(j,1) == s.DesignationArray(k,1)%same tetrode
                        bar(histVect,tempHold,'m','EdgeColor','magenta')
    %                     disp('PVPVSAME')
                    else
                        bar(histVect,tempHold,'r','EdgeColor','red')
    %                     disp('PVPVDIFF')
                    end
                elseif s.MasterSheet(j,1) == 1 & s.MasterSheet(k,1) == 0 | s.MasterSheet(j,1) == 0 & s.MasterSheet(k,1) == 1%PV vs MSN
                    if s.DesignationArray(j,1) == s.DesignationArray(k,1) %same tetrode
                        bar(histVect,tempHold,'c','EdgeColor','cyan')
    %                     disp('PVMSNSAME')
                    else
                        bar(histVect,tempHold,'g','EdgeColor','green')
    %                     disp('PVMSNDIFF')
                    end
                elseif s.MasterSheet(j,1) == 0 & s.MasterSheet(k,1) == 0 %MSN vs MSN
                    if s.DesignationArray(j,1) == s.DesignationArray(k,1) %same tetrode
                        bar(histVect,tempHold,'b','EdgeColor','blue')
    %                     disp('MSNMSNSAME')
                    else
                        bar(histVect,tempHold,'k','EdgeColor','black')
    %                     disp('MSNMSNDIFF')
                    end
                end
    %             bar(histVect,tempHold,'r')
                xlim([-0.05 0.05])
                set(gca,'ytick',[])
                set(gca,'yticklabel',[])
                if j == 1
                    if s.MasterSheet(k,1) == 1
                        title(s.DesignationName{k},'Color', 'r')
                    else
                        title(s.DesignationName{k})
                    end

                end
                if j~= numCells
                    set(gca,'xtick',[])
                    set(gca,'xticklabel',[])
                end
            end
        end

        s = [];
        
    end
    
end
% 
% bigInd = 1;
% 
% %below is code just for the analysis side of things. 
% for i = 1:numFiles
%     %open the targeted file
%     load(tester{i})
%     
%     
%     %check for any PV cells.
%     pvCheck = find(s.MasterSheet(:,1) == 1)
%     %if pv cells, do cross correlograms
%     if length(pvCheck) > 0
%         %find the number of cells
%         numCells = length(s.MasterSheet(:,1));
%         corrStore = cell(numCells,numCells);
%         for j = 1: numCells
%             j
%             for k = 1:numCells
%                 cell1 = s.(s.DesignationName{j}).SpikeTimes;
%                 cell2 = s.(s.DesignationName{k}).SpikeTimes;
%                 newInd = 1;
%                 for fIND = 1:length(cell1)
%                     subVal = cell2-cell1(fIND);
%                     subVal(subVal>0.050) = [];
%                     subVal(subVal<-0.050) = [];
%                     if length(subVal) > 0
%                         tempStore(newInd:newInd+length(subVal)-1) = subVal;
%                         newInd = newInd + length(subVal);
%                     end
%                 end
%                 corrStore{j,k} = tempStore;
%                 tempStore = [];
%                 indStore(j,k) = newInd;
%                 %also store information about what kinds of cells
%                 cellStore(j,k,1) = s.MasterSheet(j,1) ;
%                 cellStore(j,k,2) = + s.MasterSheet(k,1);
%                 %store whether same or different channel
%                 if s.DesignationArray(k,1) == s.DesignationArray(j,1)
%                     sameStore(j,k) = 1;
%                 else
%                     sameStore(j,k) = 0;
%                 end
%             end
%         end
%         C = combnk([1:numCells],2);
%         for combInd = 1:length(C)
%             bigStore{bigInd} = corrStore{C(combInd,1),C(combInd,2)};
%             infoStore(bigInd,1) = cellStore(C(combInd,1),C(combInd,2),1);
%             infoStore(bigInd,2) = cellStore(C(combInd,1),C(combInd,2),2);
%             infoStore(bigInd,3) = sameStore(C(combInd,1),C(combInd,2));
%             infoStore(bigInd,4) = C(combInd,1);
%             infoStore(bigInd,5) = C(combInd,2);
%             bigInd = bigInd + 1;
%         end
% 
%         s = [];
%         
%     end
%     
% end

%now we should have a bigStore array full of goodies, and a designation
%guide from infostore. 

%basic classify of cell types by summing

cellSum = infoStore(:,1) + infoStore(:,2);

%pull MSN MSN interactions

msnFind = find(cellSum == 0);

msnInd = 1;
for i = 1:length(msnFind)
    insertLength = length(bigStore{msnFind(i)});
    msnGroup(msnInd:msnInd+insertLength-1) = bigStore{msnFind(i)};
    msnInd = msnInd+insertLength;
end

histVect = [-0.05:0.001:0.05];
msnHist = hist(msnGroup,histVect);

figure
bar(histVect,msnHist)
title('MSN to MSN No selection')

%find pv vs pv.
pvFind = find(cellSum == 2);

pvInd = 1;
for i = 1:length(pvFind)
    insertLength = length(bigStore{pvFind(i)});
    pvGroup(pvInd:pvInd+insertLength-1) = bigStore{pvFind(i)};
    pvInd = pvInd+insertLength;
end

histVect = [-0.05:0.001:0.05];
pvHist = hist(pvGroup,histVect);

figure
bar(histVect,pvHist)
title('PV to PV No selection')

% find pv to MSN and align to the PV

mixFind = find(cellSum == 1 & infoStore(:,3) == 1);
mixFindDiff = find(cellSum == 1 & infoStore(:,3) == 0);

%make sure to align to PV cell firing! this means that PV cell has to go
%first, so if PV cell is second, then need to make negative of the other
%thing.

mixInd = 1;
for i = 1:length(mixFind)
    insertLength = length(bigStore{mixFind(i)});
    mixID = infoStore(mixFind(i),1);
    if mixID == 0
        mixGroup(mixInd:mixInd + insertLength - 1) = -bigStore{mixFind(i)};
    elseif mixID == 1
        mixGroup(mixInd:mixInd + insertLength - 1) = bigStore{mixFind(i)};
    end
    
    mixInd = mixInd + insertLength;
end
histVect = [-0.05:0.001:0.05];
mixHist = hist(mixGroup,histVect);

figure
bar(histVect,mixHist)
title('PV to MSN Same selection')


mixDiffInd = 1;
for i = 1:length(mixFindDiff)
    insertLength = length(bigStore{mixFindDiff(i)});
    mixDiffID = infoStore(mixFindDiff(i),1);
    if mixDiffID == 0
        mixDiffGroup(mixDiffInd:mixDiffInd + insertLength - 1) = -bigStore{mixFindDiff(i)};
    elseif mixDiffID == 1
        mixDiffGroup(mixDiffInd:mixDiffInd + insertLength - 1) = bigStore{mixFindDiff(i)};
    end
    
    mixDiffInd = mixDiffInd + insertLength;
end
histVect = [-0.05:0.001:0.05];
mixDiffHist = hist(mixDiffGroup,histVect);

figure
bar(histVect,mixDiffHist)
title('PV to MSN Different selection')
%soooo it looks like I cant see any suppression at all. 

%I think its time to end this diversion for now. Ultimately, should try
%applying this to the tuning curve data or something. 