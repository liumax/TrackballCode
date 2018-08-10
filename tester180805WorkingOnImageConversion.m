
cd /Users/maxliu/Desktop/180805ImagesForDeepCut
%extract names
fullSet = dir;
names = {fullSet.name};
names = names';
%remove crap names
names([1,2]) = [];
%now we need to make sure things are in the right order! 
[cs,index] = sort_nat(names);

%now we need to do image conversion 
for i = 1:length(cs)
    I=imread(cs{i});
    if i < 10
        imName = strcat('img000',num2str(i),'.png');
    elseif i < 100 & i >= 10
        imName = strcat('img00',num2str(i),'.png');
    elseif i < 1000 & i >=100
        imName = strcat('img0',num2str(i),'.png');
    elseif i > 1000
        imName = strcat('img',num2str(i),'.png');
    end
    imwrite(I,imName);
end





