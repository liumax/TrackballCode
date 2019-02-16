

% bigString = [];
strInd = 1;
for i = 28:38
    if i < 10
        newStr = strcat('ffmpeg -i "Trial',{'     '},num2str(i),'.mpg" -filter:v "crop=611:656:55:150" Trial',num2str(i),'Arena1.mp4 & ffmpeg -i "Trial',{'     '},num2str(i),'.mpg" -filter:v "crop=611:656:666:150" Trial',num2str(i),'Arena2.mp4 & ');
%         newStr = newStr{1};
    elseif i >= 10
        newStr = strcat('ffmpeg -i "Trial',{'    '},num2str(i),'.mpg" -filter:v "crop=611:656:55:150" Trial',num2str(i),'Arena1.mp4 & ffmpeg -i "Trial',{'    '},num2str(i),'.mpg" -filter:v "crop=611:656:666:150" Trial',num2str(i),'Arena2.mp4 & ');
%         newStr = newStr{1};
    end
    newStr = newStr{1};
    newLength = length(newStr);
    bigString(strInd:strInd+newLength-1) = newStr;
    strInd = strInd + newLength;
end
