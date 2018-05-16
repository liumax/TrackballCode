function [timestamps] = readCameraModuleTimeStamps(filename)  
%[timestamps] = readCameraModuleTimeStamps(filename)
%filename-- a string containing the name of the .videoTimeStames file
%timestamps-- a vector containing the timestamps of the camera frames (in uint32 format)

clockRate = 30000;  %default clock rate
fid = fopen(filename,'r');
headerText = fread(fid,200,'char');
headerText = headerText';
endHeaderLoc = strfind(headerText,'<End settings>');

if (~isempty(endHeaderLoc))
    headersize = endHeaderLoc+14;
    clockRateLoc  = strfind(headerText,'Clock rate:');
    if (~isempty(clockRateLoc))
        clockRate = str2num(char(strtok(headerText(clockRateLoc+12:end))));
    end
    
else
    headersize = 0;
end
frewind(fid);
junk = fread(fid,headersize,'char');

timestamps = fread(fid,inf,'uint32=>double',0)/clockRate;


fclose(fid);


