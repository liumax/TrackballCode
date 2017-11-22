
function [zRaster,baselineMean,baselineSTD] = functionZScore(rasterInput,zeroPoint,numTrials);

%first, find out how the raster input is oriented

rasterSize = size(rasterInput);
orientFind = find(rasterSize == numTrials);
%if orientFind == 2, that means that individual columns are trials. If 1,
%then rows are trials
if orientFind == 1
    newRaster = rasterInput';
else
    newRaster = rasterInput;
end

%now pull baseline period and linearize
baselineVector = reshape(newRaster(1:zeroPoint,:),1,[]);
%calculate mean
baselineMean = mean(baselineVector);
baselineSTD = std(baselineVector);

%make z-scored raster
zRaster = (rasterInput - baselineMean)./baselineSTD;




end