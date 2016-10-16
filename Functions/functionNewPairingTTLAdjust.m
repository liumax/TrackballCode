%This code edits the TTL information from pairing, since the laser signal
%will make alignment more complicated than aligning to every TTL pulse. It
%does so by determining which TTLs are close to each other by the
%LaserTriggerPulseITI, and then removing those. This code isolates two
%cases: where laser is coincident with the tone, and therefore the number
%of TTLs is twice the number of tone presentations, or when the laser is
%outside of the tone, in which case there are thrice as many TTLs.

%Inputs: 
%s: structured array
%soundName: name of the target sound structure. Should be soundNames{3}

%Outputs: 
%s: structured array. Adjusted s.SoundData.soundName.


function [s] = functionNewPairingTTLAdjust(s,soundName);
TTLs = s.TTLs.(soundName);
numTTLs = size(TTLs,1);
expectedTTLs = s.SoundData.(soundName).ToneRepetitions;
%find all TTL diffs.
TTLdiff = diff(TTLs);
%pulls the appropriate ITI from the sound file
targetITI =  s.SoundData.(soundName).LaserTriggerPulseITI;
%since things arent perfect, adds buffer so things are acceptable.
acceptableRange = [0.9 1.1];
acceptableRange = acceptableRange*targetITI;
%finds all TTLs that fit the acceptable range. NOTE THAT THIS CODE WILL
%BREAK IF THE LASER CUE IS CLOSE ENOUGH TO THE TONE (so that it picks
%up that ITI as well)
laserITIFinder = find(TTLdiff >acceptableRange(1) & TTLdiff<acceptableRange(2));


if numTTLs == 2* expectedTTLs %This condition identifies trials in which the laser and the tone are aligned. 
%As a result, we produce twice the expected TTLs, since there is a double
%TTL pulse. 
    
    %removes all TTLs associated with the laser.
    for i = size(laserITIFinder,1):-1:1
        TTLs(laserITIFinder(i)+1) = [];
    end
    %check to make sure its okay
    if size(TTLs,1) ~= expectedTTLs
        error('TTL SUBTRACTION FAILED')
    end
    %replaces TTL signals in the masterstructure.
    s.TTLs.(strcat(soundName,'Original')) = s.TTLs.(soundName);
    s.TTLs.(soundName) = TTLs;

elseif numTTLs == 3*expectedTTLs %This is the scenario in which laser is not coincident with the tone.
    laserTTLs = sort(reshape([laserITIFinder,laserITIFinder+1],length(laserITIFinder)*2,1));
    %removes all TTLs associated with the laser.
    for i = size(laserTTLs,1):-1:1
        TTLs(laserTTLs(i)) = [];
    end
    
    %check to make sure its okay
    if size(TTLs,1) ~= expectedTTLs
        error('TTL SUBTRACTION FAILED')
    end
    
    %replaces TTL signals in the masterstructure.
    s.TTLs.(strcat(soundName,'Original')) = s.TTLs.(soundName);
    s.TTLs.(soundName) = TTLs;

else
    error('SOMETHING IS WRONG WITH TTL INPUTS: INCORRECT NUMBER OF INPUTS')
end

end