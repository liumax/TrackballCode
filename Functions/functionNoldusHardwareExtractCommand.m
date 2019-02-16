function [timeStore,stateStore] = functionNoldusHardwareExtractCommand(fname);


fid = fopen(fname);
%extract header info
tline = fgets(fid); 
findColon = strfind(tline,';');
headerNum = str2num(tline(findColon(1)+2:findColon(2)-2));
%blow past header
for i = 1:headerNum-1
    tline = fgets(fid);
end
bigCount = 1;
while ischar(tline) 
    %pull new line
    tline = fgets(fid);
    if tline == -1
        break
    end
    %check for command. 
    findColon = strfind(tline,';');
    %command is found between third and fourth semicolon
    tarText = tline(findColon(3)+2:findColon(4)-2);
    if strcmp(tarText,'command')
        %store time
        timeStore(bigCount) = str2num(tline(findColon(1)+1:findColon(2)-1));
        tarText = tline(findColon(end-1)+2:findColon(end)-2);
        findSpace = strfind(tarText,' ');
        stateStore{bigCount} = tarText(findSpace(end)+1:end);
        bigCount = bigCount + 1;
    end
end










































