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


function [s] = functionABBAPairingTTLAdjust(s,soundName);
TTLs = s.TTLs.(soundName);
numTTLs = size(TTLs,1);
expectedTTLs = s.SoundData.LaserRepetitions;
%find all TTL diffs.
TTLdiff = diff(TTLs);
%pulls the appropriate ITI from the sound file
targetITI =  s.SoundData.LaserITI;
%since things arent perfect, adds buffer so things are acceptable.
acceptableRange = [0.9 1.1];
acceptableRange = acceptableRange*targetITI;
%finds all TTLs that fit the acceptable range. NOTE THAT THIS CODE WILL
%BREAK IF THE LASER CUE IS CLOSE ENOUGH TO THE TONE (so that it picks
%up that ITI as well)
laserITIFinder = find(TTLdiff >acceptableRange(1) & TTLdiff<acceptableRange(2));


if numTTLs == 3* expectedTTLs %This condition identifies trials in which the laser and the tone are aligned. 
%As a result, we produce 3x the correct number of TTLs, since there are two
%tones, and one of them gets double the TTLs expected.
    
    %removes all TTLs associated with the laser.
    for i = size(laserITIFinder,1):-1:1
        TTLs(laserITIFinder(i)) = [];
    end
    %check to make sure its okay
    if size(TTLs,1) ~= 2*expectedTTLs
        error('TTL SUBTRACTION FAILED')
    end
    %replaces TTL signals in the masterstructure.
    s.TTLs.(strcat(soundName,'Original')) = s.TTLs.(soundName);
    s.TTLs.(soundName) = TTLs;
    disp('Corrected TTLs for Coincident Pairing')

elseif numTTLs == 4*expectedTTLs %This is the scenario in which laser is not coincident with the tone.
    laserTTLs = sort(reshape([laserITIFinder,laserITIFinder+1],length(laserITIFinder)*2,1));
    %removes all TTLs associated with the laser.
    for i = size(laserTTLs,1):-1:1
        TTLs(laserTTLs(i)) = [];
    end
    
    %check to make sure its okay
    if size(TTLs,1) ~= 2*expectedTTLs
        error('TTL SUBTRACTION FAILED')
    end
    
    %replaces TTL signals in the masterstructure.
    s.TTLs.(strcat(soundName,'Original')) = s.TTLs.(soundName);
    s.TTLs.(soundName) = TTLs;
    disp('Corrected TTLs for Pre/Post Pairing')
elseif numTTLs == 2*expectedTTLs %This is scenario with no laser
    laserTTLs = sort(TTLs);
    disp('No Correction for TTLs Needed')
else
    disp(strcat('Actual:',num2str(numTTLs)))
    disp(strcat('Expected:',num2str(expectedTTLs)))
    error('SOMETHING IS WRONG WITH TTL INPUTS: INCORRECT NUMBER OF INPUTS')
end

end