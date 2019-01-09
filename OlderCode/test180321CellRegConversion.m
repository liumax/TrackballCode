%this is code to extract filters for the purposes of using cellreg. 




files = what;
files = files.mat;

for i = 1:length(files)
    load(files{i})
    imSize = size(neuron_results.Cn);
    numCells = size(neuron_results.A,2);
    allFiltersMat = zeros(numCells,imSize(1),imSize(2));
    for j = 1:numCells
        allFiltersMat(j,:,:) = reshape(neuron_results.A(:,j),imSize(1),imSize(2));
    end
    targetName = files{i}(1:end-4);
    targetName = strcat(targetName,'ROIs.mat');
    save(targetName,'allFiltersMat')
end







