%%This function is meant to plot raster information (in the form of one
%%column of times, and one column of trial numbers) out as vertical lines. 



%inputs: times, trials

function [] = rasterPlot(times,trials)


for rasterInd = 1:length(times)
    plot([times(rasterInd) times(rasterInd)],[trials(rasterInd) trials(rasterInd) + 1])
end

end