%This is code meant to compare tuning from before and after.


%pull all mat files.
names = what;
names = names.mat;

%exclude the sounds and TMP
[findString] = functionCellStringFind(names,'Sound');
names(findString) = [];

[findString] = functionCellStringFind(names,'TMP');
names(findString) = [];

%select for analysis files and tuning files.
[findString] = functionCellStringFind(names,'Analysis');
names = names(findString);

[findString] = functionCellStringFind(names,'Tuning');
names = names(findString);

%pull information from analysis files.
for i = 1:length(names)
    targetName = names{i};
    animalNames{i} =  targetName(1:9);
    dates(i) = str2num(targetName(11:16));
end

%now we need to figure out which animals are the same, and figure out which
%recording came before/after.

uniqueNames = unique(animalNames);
desigArray = zeros(length(names),2);

for i = 1:length(uniqueNames)
    findStrings = strfind(animalNames,uniqueNames{i});
    Index = find(not(cellfun('isempty', findStrings)));
    desigArray(Index,1) = i;
    desigArray(Index,2) = [1:length(Index)];
end
desigArray(:,3) = [1:length(names)];


%now that we have a desig array to guide the extraction of files. lets do
%it! Go animal by animal

for i = 1:length(uniqueNames)
    i
    %find the target files
    targetInds = desigArray(desigArray(:,1) == i,3);
    %open the pre file. 
    load(names{targetInds(1)})
    prePhotoAverage = s.Processed.PhotoAverages;
    preRaster = s.Processed.RiseRaster;
    preDBSort = s.SoundData.DBSort;
    preMatrix = s.SoundData.TrialMatrix;
    preVelAverage = s.Processed.VelAverages;
    %load the post file
    load(names{targetInds(2)})
    postPhotoAverage = s.Processed.PhotoAverages;
    postRaster = s.Processed.RiseRaster;
    postDBSort = s.SoundData.DBSort;
    postMatrix = s.SoundData.TrialMatrix;
    postVelAverage = s.Processed.VelAverages;
    
    %pull out general things
    rasterWindow = s.Parameters.RasterWindow;
    
    rasterAxis = zeros(rasterWindow(2)-rasterWindow(1)+1,2);
    rasterAxis(:,1) = [rasterWindow(1):rasterWindow(2)];
    %find a second in photometry sample time
    rasterAxis(:,2) = [0:length(postPhotoAverage)/8:length(postPhotoAverage)];
    rasterAxis(1,2) = 1;
    
    numFreqs = size(s.Processed.PhotoAverages,2);
    uniqueFreqs = unique(postMatrix(:,2));
    
    
    % create a set of ticks for labeling any display of octaves
    if isfield(s.SoundData,'WhiteNoise')
        if s.SoundData.WhiteNoise == 1
            totalOctaves = round(log2(uniqueFreqs(end)/uniqueFreqs(2)));

            %This then makes an array of the full octave steps I've made
            octaveRange = zeros(totalOctaves + 1,2);
            octaveRange(1,1) = uniqueFreqs(2);
            for ind = 1:totalOctaves
                octaveRange (ind+1,1) = octaveRange(ind,1)*2;
            end
            %next, I find the positions from uniqueFreqs that match octaveRange
            for ind = 1:size(octaveRange,1);
                octaveRange(ind,2) = find(round(uniqueFreqs) == octaveRange(ind,1));
            end
            %add in the white noise
            octaveRange(2:end+1,:) = octaveRange(1:end,:);
            octaveRange(1,1) = 0;
            octaveRange(1,2) = 1;
        else soundData.WhiteNoise == 0;
            totalOctaves = round(log2(uniqueFreqs(end)/uniqueFreqs(1)));

            %This then makes an array of the full octave steps I've made
            octaveRange = zeros(totalOctaves + 1,2);
            octaveRange(1,1) = uniqueFreqs(1);
            for ind = 1:totalOctaves
                octaveRange (ind+1,1) = octaveRange(ind,1)*2;
            end
            %next, I find the positions from uniqueFreqs that match octaveRange
            for ind = 1:size(octaveRange,1);
                octaveRange(ind,2) = find(round(uniqueFreqs) == octaveRange(ind,1));
            end
        end
    else
        totalOctaves = round(log2(uniqueFreqs(end)/uniqueFreqs(1)));

        %This then makes an array of the full octave steps I've made
        octaveRange = zeros(totalOctaves + 1,2);
        octaveRange(1,1) = uniqueFreqs(1);
        for ind = 1:totalOctaves
            octaveRange (ind+1,1) = octaveRange(ind,1)*2;
        end
        %next, I find the positions from uniqueFreqs that match octaveRange
        for ind = 1:size(octaveRange,1);
            octaveRange(ind,2) = find(round(uniqueFreqs) == octaveRange(ind,1));
        end
    end
    
    %now lets do some basic subtractions
    subPhotoAverage = postPhotoAverage - prePhotoAverage;
    subVelAverage = postVelAverage - preVelAverage;
    
    %now lets find max and min values across both sets.
    preMaxPhot = max(max(max(prePhotoAverage)));
    preMinPhot = min(min(min(prePhotoAverage)));
    postMaxPhot = max(max(max(postPhotoAverage)));
    postMinPhot = min(min(min(postPhotoAverage)));
    %do the same for vels
    preMaxVel = max(max(max(preVelAverage)));
    preMinVel = min(min(min(preVelAverage)));
    postMaxVel = max(max(max(postVelAverage)));
    postMinVel = min(min(min(postVelAverage)));
    
    %find full maximum and minimum
    fullMaxPhot = max([preMaxPhot,postMaxPhot]);
    fullMinPhot = min([preMinPhot,postMinPhot]);
    
    fullMaxVel = max([preMaxVel,postMaxVel]);
    fullMinVel = min([preMinVel,postMinVel]);
    
    %find max and min for velocity and photometry subtractions
    subMaxPhot = max(max(max(subPhotoAverage)));
    subMinPhot = min(min(min(subPhotoAverage)));
    subMaxVel = max(max(max(subVelAverage)));
    subMinVel = min(min(min(subVelAverage)));
    
    %time to do some plots.
    
    subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.04], [0.03 0.05], [0.03 0.01]);
    
    %first, plot out changes to photometry signal. We will plot the pre,
    %the post, and then the subtraction
    
    hFig = figure;
    
    %plot PRE
    subplot(3,3,1)
    imagesc(squeeze(prePhotoAverage(:,:,1))',[fullMinPhot fullMaxPhot])
    colormap('parula')
    colorbar
    set(gca,'XTick',rasterAxis(:,2));
    set(gca,'XTickLabel',rasterAxis(:,1));
    set(gca,'YTick',octaveRange(:,2));
    set(gca,'YTickLabel',octaveRange(:,1));
    title({'Pre',uniqueNames{i},'60DB'})
    
    subplot(3,3,4)
    imagesc(squeeze(prePhotoAverage(:,:,2))',[fullMinPhot fullMaxPhot])
    colormap('parula')
    colorbar
    set(gca,'XTick',rasterAxis(:,2));
    set(gca,'XTickLabel',rasterAxis(:,1));
    set(gca,'YTick',octaveRange(:,2));
    set(gca,'YTickLabel',octaveRange(:,1));
    title({'Pre',uniqueNames{i},'80DB'})
    
    subplot(3,3,7)
    imagesc(squeeze(prePhotoAverage(:,:,3))',[fullMinPhot fullMaxPhot])
    colormap('parula')
    colorbar
    set(gca,'XTick',rasterAxis(:,2));
    set(gca,'XTickLabel',rasterAxis(:,1));
    set(gca,'YTick',octaveRange(:,2));
    set(gca,'YTickLabel',octaveRange(:,1));
    title({'Pre',uniqueNames{i},'100DB'})
    
    %PLOT POST
    subplot(3,3,2)
    imagesc(squeeze(postPhotoAverage(:,:,1))',[fullMinPhot fullMaxPhot])
    colormap('parula')
    colorbar
    set(gca,'XTick',rasterAxis(:,2));
    set(gca,'XTickLabel',rasterAxis(:,1));
    set(gca,'YTick',octaveRange(:,2));
    set(gca,'YTickLabel',octaveRange(:,1));
    title({'Post',uniqueNames{i},'60DB'})
    
    subplot(3,3,5)
    imagesc(squeeze(postPhotoAverage(:,:,2))',[fullMinPhot fullMaxPhot])
    colormap('parula')
    colorbar
    set(gca,'XTick',rasterAxis(:,2));
    set(gca,'XTickLabel',rasterAxis(:,1));
    set(gca,'YTick',octaveRange(:,2));
    set(gca,'YTickLabel',octaveRange(:,1));
    title({'Post',uniqueNames{i},'80DB'})
    
    subplot(3,3,8)
    imagesc(squeeze(postPhotoAverage(:,:,3))',[fullMinPhot fullMaxPhot])
    colormap('parula')
    colorbar
    set(gca,'XTick',rasterAxis(:,2));
    set(gca,'XTickLabel',rasterAxis(:,1));
    set(gca,'YTick',octaveRange(:,2));
    set(gca,'YTickLabel',octaveRange(:,1));
    title({'Post',uniqueNames{i},'100DB'})
    
    %PLOT SUBTRACTION
    subplot(3,3,3)
    imagesc(squeeze(subPhotoAverage(:,:,1))',[subMinPhot subMaxPhot])
    colormap('parula')
    colorbar
    set(gca,'XTick',rasterAxis(:,2));
    set(gca,'XTickLabel',rasterAxis(:,1));
    set(gca,'YTick',octaveRange(:,2));
    set(gca,'YTickLabel',octaveRange(:,1));
    title({'Sub',uniqueNames{i},'60DB'})
    
    subplot(3,3,6)
    imagesc(squeeze(subPhotoAverage(:,:,2))',[subMinPhot subMaxPhot])
    colormap('parula')
    colorbar
    set(gca,'XTick',rasterAxis(:,2));
    set(gca,'XTickLabel',rasterAxis(:,1));
    set(gca,'YTick',octaveRange(:,2));
    set(gca,'YTickLabel',octaveRange(:,1));
    title({'Sub',uniqueNames{i},'80DB'})
    
    subplot(3,3,9)
    imagesc(squeeze(subPhotoAverage(:,:,3))',[subMinPhot subMaxPhot])
    colormap('parula')
    colorbar
    set(gca,'XTick',rasterAxis(:,2));
    set(gca,'XTickLabel',rasterAxis(:,1));
    set(gca,'YTick',octaveRange(:,2));
    set(gca,'YTickLabel',octaveRange(:,1));
    title({'Sub',uniqueNames{i},'100DB'})
    
    
    hFig = figure;
    subplot(3,3,1)
    imagesc(squeeze(prePhotoAverage(:,:,1))')
    colormap('parula')
    colorbar
    set(gca,'XTick',rasterAxis(:,2));
    set(gca,'XTickLabel',rasterAxis(:,1));
    set(gca,'YTick',octaveRange(:,2));
    set(gca,'YTickLabel',octaveRange(:,1));
    title({'Pre',uniqueNames{i},'60DB'})
    
    subplot(3,3,4)
    imagesc(squeeze(prePhotoAverage(:,:,2))')
    colormap('parula')
    colorbar
    set(gca,'XTick',rasterAxis(:,2));
    set(gca,'XTickLabel',rasterAxis(:,1));
    set(gca,'YTick',octaveRange(:,2));
    set(gca,'YTickLabel',octaveRange(:,1));
    title({'Pre',uniqueNames{i},'80DB'})
    
    subplot(3,3,7)
    imagesc(squeeze(prePhotoAverage(:,:,3))')
    colormap('parula')
    colorbar
    set(gca,'XTick',rasterAxis(:,2));
    set(gca,'XTickLabel',rasterAxis(:,1));
    set(gca,'YTick',octaveRange(:,2));
    set(gca,'YTickLabel',octaveRange(:,1));
    title({'Pre',uniqueNames{i},'100DB'})
    
    subplot(3,3,2)
    imagesc(squeeze(postPhotoAverage(:,:,1))')
    colormap('parula')
    colorbar
    set(gca,'XTick',rasterAxis(:,2));
    set(gca,'XTickLabel',rasterAxis(:,1));
    set(gca,'YTick',octaveRange(:,2));
    set(gca,'YTickLabel',octaveRange(:,1));
    title({'Post',uniqueNames{i},'60DB'})
    
    subplot(3,3,5)
    imagesc(squeeze(postPhotoAverage(:,:,2))')
    colormap('parula')
    colorbar
    set(gca,'XTick',rasterAxis(:,2));
    set(gca,'XTickLabel',rasterAxis(:,1));
    set(gca,'YTick',octaveRange(:,2));
    set(gca,'YTickLabel',octaveRange(:,1));
    title({'Post',uniqueNames{i},'80DB'})
    
    subplot(3,3,8)
    imagesc(squeeze(postPhotoAverage(:,:,3))')
    colormap('parula')
    colorbar
    set(gca,'XTick',rasterAxis(:,2));
    set(gca,'XTickLabel',rasterAxis(:,1));
    set(gca,'YTick',octaveRange(:,2));
    set(gca,'YTickLabel',octaveRange(:,1));
    title({'Post',uniqueNames{i},'100DB'})
    
    subplot(3,3,3)
    imagesc(squeeze(subPhotoAverage(:,:,1))')
    colormap('parula')
    colorbar
    set(gca,'XTick',rasterAxis(:,2));
    set(gca,'XTickLabel',rasterAxis(:,1));
    set(gca,'YTick',octaveRange(:,2));
    set(gca,'YTickLabel',octaveRange(:,1));
    title({'Sub',uniqueNames{i},'60DB'})
    
    subplot(3,3,6)
    imagesc(squeeze(subPhotoAverage(:,:,2))')
    colormap('parula')
    colorbar
    set(gca,'XTick',rasterAxis(:,2));
    set(gca,'XTickLabel',rasterAxis(:,1));
    set(gca,'YTick',octaveRange(:,2));
    set(gca,'YTickLabel',octaveRange(:,1));
    title({'Sub',uniqueNames{i},'80DB'})
    
    subplot(3,3,9)
    imagesc(squeeze(subPhotoAverage(:,:,3))')
    colormap('parula')
    colorbar
    set(gca,'XTick',rasterAxis(:,2));
    set(gca,'XTickLabel',rasterAxis(:,1));
    set(gca,'YTick',octaveRange(:,2));
    set(gca,'YTickLabel',octaveRange(:,1));
    title({'Sub',uniqueNames{i},'100DB'})
end














