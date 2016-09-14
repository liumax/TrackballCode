%% This is meant to be generalized code for making rasters
%spikeTimes should be an array of n x 1, where n is the number of spikes.

%alignTimes should be a n x 1 vector, with n being the number of events to
%align. 

%rasterWindow is a two element vector, [x1 x2] specifying the limits of
%what you want in the raster.

function [rasters] = functionBasicRaster(spikeTimes,alignTimes,rasterWindow);

rasters = zeros(100000,2);
holder = zeros(size(spikeTimes,1),2);
rasterCounter = 1;

for alignCounter = 1:size(alignTimes,1)
    holder = spikeTimes;%allocates spikeTimes in.
    holder(:,1) = spikeTimes(:,1)-alignTimes(alignCounter); %subtracts alignment time
    holder(holder(:,1)<rasterWindow(1) | holder(:,1)>rasterWindow(2),:) = []; %removes data entries that are larger than the targeted window
    holder(:,2) = alignCounter;
    holderSize = size(holder); %sizes remaining entries
    if holderSize(1) > 0
        rasters(rasterCounter:rasterCounter + holderSize(1) - 1,:) = holder; %inserts these entries. 
        rasterCounter = rasterCounter + holderSize(1); %updates counter
    end
end
rasters(rasters(:,1) == 0,:) = [];

end