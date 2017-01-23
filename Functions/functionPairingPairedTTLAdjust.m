%This code is supposed to edit the TTL pulses for pairing sessions. Pairing
%sessions should have an increased number of TTL pulses, so this is here to
%adjust for that. This will remove all excess TTL pulses.


function [masterStruct] = functionPairingPairedTTLAdjust(masterStruct,ttlName5,...
    soundFile,soundName5);
TTLs = masterStruct.TTLs.(ttlName5);
numTTLs = size(TTLs,1);
expectedTTLs = soundFile.PairingRepetitions;
if numTTLs == 2* expectedTTLs %in this condition, the laser tones are clearly
    %separated from the tone TTLs, which means that there are double the
    %number of TTLs as expected from the number of tone presentations.
    
    %find all TTL diffs.
    TTLdiff = diff(TTLs);
    %pulls the appropriate ITI from the sound file
    targetITI = soundFile.(soundName5).LaserTriggerPulseITI;
    %since things arent perfect, adds buffer so things are acceptable.
    acceptableRange = [0.9 1.1];
    acceptableRange = acceptableRange*targetITI;
    %finds all TTLs that fit the acceptable range. NOTE THAT THIS CODE WILL
    %BREAK IF THE LASER CUE IS CLOSE ENOUGH TO THE TONE (so that it picks
    %up that ITI as well)
    laserITIFinder = find(TTLdiff >acceptableRange(1) & TTLdiff<acceptableRange(2));
    %removes all TTLs associated with the laser.
    for i = size(laserITIFinder,1):-1:1
        TTLs(laserITIFinder(i):laserITIFinder(i)+1) = [];
    end
    %check to make sure its okay
    if size(TTLs,1) ~= soundFile.PairingRepetitions
        error('TTL SUBTRACTION FAILED')
    end
    %replaces TTL signals in the masterstructure.
    masterStruct.TTLs.(ttlName5) = TTLs;

elseif numTTLs == 3*expectedTTLs/2 %This is the scenario in which laser is 
    %turned on at the same time as the tone.
    %find all TTL diffs.
    TTLdiff = diff(TTLs);
    %pulls the appropriate ITI from the sound file
    targetITI = soundFile.(soundName5).LaserTriggerPulseITI;
    %since things arent perfect, adds buffer so things are acceptable.
    acceptableRange = [0.9 1.1];
    acceptableRange = acceptableRange*targetITI;
    %finds all TTLs that fit the acceptable range. NOTE THAT THIS CODE WILL
    %BREAK IF THE LASER CUE IS CLOSE ENOUGH TO THE TONE (so that it picks
    %up that ITI as well)
    laserITIFinder = find(TTLdiff >acceptableRange(1) & TTLdiff<acceptableRange(2));
    %removes all TTLs associated with the laser.
    for i = size(laserITIFinder,1):-1:1
        TTLs(laserITIFinder(i)) = [];
    end
    
    %check to make sure its okay
    if size(TTLs,1) ~= soundFile.PairingRepetitions
        error('TTL SUBTRACTION FAILED')
    end
    
    %replaces TTL signals in the masterstructure.
    masterStruct.TTLs.(ttlName5) = TTLs;
else
    error('SOMETHING IS WRONG WITH TTL INPUTS: INCORRECT NUMBER OF INPUTS')
end

end