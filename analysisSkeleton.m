%This is meant to be the analysis code for pairing experiments. This code
%should analyze the tuning curves before and after, the playing of long
%tones before and after, and the actual changes that may occur with laser
%light itself. 

%Establishes the folder and subfolders for analysis. Adds all folders to
%path for easy access to files.
homeFolder = pwd;
subFolders = genpath(homeFolder);
addpath(subFolders)
subFoldersCell = strsplit(subFolders,';')';


%First thing is to import the sound file data.

soundName = strcat(fileName,'.mat');
soundFile = open(soundName);



%Next I need to find the pulses!

%extract DIO filenames
[D1FileName] = functionFileFinder(subFoldersCell,'DIO','D1');
D1FileName = D1FileName{1};
[D2FileName] = functionFileFinder(subFoldersCell,'DIO','D2');
D2FileName = D2FileName{1};

%extract actual data
[DIO1Data] = readTrodesExtractedDataFile(D1FileName);
[DIO2Data] = readTrodesExtractedDataFile(D2FileName);

%remove just time and state data. convert to double for easier processing
DIO1Data = double(DIO1Data.fields(1).data,DIO1Data.fields(2).data);
DIO2Data = double(DIO2Data.fields(1).data,DIO2Data.fields(2).data);

%pull all differences in state for DIO data (diff subtracts x(2)-x(1)). +1
%fudge factor adjusts for indexing. MUST BE DOUBLE!!
DIO1Diff = find(diff(DIO1Data(:,2))==1)+1;
DIO1High = find(DIO1Data(:,2) == 1);
%finds all points which are down to up states, extracts these times.
DIO1True = intersect(DIO1Diff,DIO1High);
DIO1True = DIO1Data(DIO1True,1);
%finds differences between time points
DIO1TrueDiff = diff(DIO1True);

%%%VARIABLES HERE%%%
signalITI = 20; %signal TTL ITI in msecs

%%%Time Periods in Recording%%%
baselineTimes = 
tuningFirstTimes = 
presentationFirstTimes = 
pairingTimes = 
presentationSecondTimes = 
tuningSecondTimes = 

baselineTTLs = 
tuningFirstTTLs = 
presentationFirstTTLs = 
pairingTTLs = 
presentationSecondTTLs = 
tuningSecondTTLs = 






