
numUnits = length(s.DesignationName); %pulls number of units
numCross = nchoosek(numUnits,2); %pulls # of combinations 
crossDesig = nchoosek([1:1:numUnits],2);
crossCorrWindow = 0.05; %window in ms that I want to look at with cross correlograms
crossCorrSampling = 30000; %number of bins for the crossCorr window. Currently set to sample at sampling freq
crossCorrVector = [-crossCorrWindow:2*crossCorrWindow/crossCorrSampling:crossCorrWindow];

crossCorrStore = cell(numCross,1);
zeroVals = zeros(numCross,1);

for i = 1:numCross
    x = s.(s.DesignationName{crossDesig(i,1)}).SpikeTimes; %pulls spike times from one unit
    y = s.(s.DesignationName{crossDesig(i,2)}).SpikeTimes; %pulls spike times from another unit
    
    if ~isempty(length(intersect(x,y)))
        zeroVals(i) = length(intersect(x,y));
    end
    
    [X Y] = meshgrid(x,y); 
    %generates two grids: one with all values of x tiled in length(y) rows, and
    %all values of y, tiled in length(x) columns. This generates two matrices,
    %both of which are length(y) x length(x). 

    test1 = X-Y; %subtract each matrix both ways to generate symmetry
    test2 = Y- X; %subtract each matrix both ways to generate symmetry
    smallVals = [test1(abs(test1)<crossCorrWindow);test2(abs(test2)<crossCorrWindow)]; %isolate only values below 50ms
    crossCorrStore{i} = hist(smallVals,crossCorrVector)/2;%plot with bin size = sampling rate of Trodes. /2 acts as fudge factor so correct # comes out
end

%find all squares
squareBase = [1:1:20]'; %generates a number of bases for squares
squareSize = arrayfun(@(x) x^2, squareBase)-numCross; %squares those numbers, then subtracts the number of cross correlograms
nearSquare = find(squareSize>=0,1,'first'); %finds nearest larger square
subPlotSize = squareBase(nearSquare);

make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.02], [0.01 0.05], [0.01 0.01]);
if ~make_it_tight,  clear subplot;  end

hFig = figure;
set(hFig, 'Position', [10 80 1240 850])

for i = 1:numCross
    subplot(subPlotSize,subPlotSize,i)
    bar(crossCorrVector,crossCorrStore{i})
    hold on
    bar(0,zeroVals(i),'r')
    xlim([-crossCorrWindow,crossCorrWindow])
    title(strcat(s.DesignationName{crossDesig(i,1)},'vs',s.DesignationName{crossDesig(i,2)}))
end
set(findall(gcf,'-property','FontSize'),'FontSize',8) %8 is minimum legible size


