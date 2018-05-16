%
%function [] = make_ripple_trigger(outfile, filelength, fs, nbits, initbuffer, outbuffer)
%
%       FILE NAME       : make_trigger
%       DESCRIPTION     : Generates an 'int16' trigger file
%
% outfile : Output file name. The output will be signed word
%           format, i.e. .sw format.
%
% filelength : File Length (Number of Samples). If filelength 
%              is a string then make_trigger.m will open the 
%              file specified as the string and find its length.
%              In this case the file must be of type float. You should
%              use the *.bin file for this input.
%
% fs : sampling rate of the sound file
%
% nbits : signed integer resolution of the sound file. Usually 16.
%
% initbuffer : Number of initial zeros pre-appended to a trigger 
%              file to maintain consistency with stimulus files.
%
% outbuffer : Number of zeros to append to a trigger 
%              file to maintain consistency with stimulus files.
%
% caa 9/16/03

function [] = make_ripple_trigger(outfile, filelength, fs, nbits, initbuffer, outbuffer)

if ( nargin ~= 6 )
   error('You need 6 input args.');
end

% if ( nargin == 2 )
%    nbits = 16;
%    initbuffer = 0;
%    outbuffer = 0;
% end
% 
% if ( nargin == 3 )
%    initbuffer = 0;
%    outbuffer = 0;
% end
% 
% if ( nargin == 4 )
%    outbuffer = 0;
% end


if ( nbits == 16 )
   dtype = 'int16';
elseif ( nbits == 20 )
   dtype = 'bit20';
elseif ( nbits == 24 )
   dtype = 'bit24';
else
   error('nbits must be either 16, 20, or 24.');
end


if ( isstr(filelength) )  % find length of sound file
   totlength = 0;
   fid = fopen(filelength,'r');
   while ~feof(fid)
      temp = fread(fid, 1024*128, 'float');
      totlength = totlength + length(temp);
   end
   fclose(fid);
else
   totlength = filelength;
end % if


% specify and create the trigger exactly

triglength = 32000;  % inter-trigger time MUST BE 32000 samples

pulselength = round(0.05*fs) + mod(round(0.05*fs), 2);  % pulse-length will be 50 ms

trigger = 2^(nbits-1)/2 * [ones(1, pulselength/2) zeros(1,pulselength/2) zeros(1, triglength-pulselength)]; %180320 max edit. eliminate negative. 
% trigger = 2^(nbits-1)/2 * [ones(1, pulselength/2) -1*ones(1,pulselength/2) zeros(1, triglength-pulselength)];

% make the special first trigger, which will have 3 distinct pulses
% firsttrigger = 2^(nbits-1)/2 * [ones(1, pulselength/2) -1*ones(1,pulselength/2) zeros(1, triglength-pulselength)];
% firsttrigger(5000:(5000+pulselength-1)) = 2^(nbits-1)/2 * [ones(1, pulselength/2) -1*ones(1,pulselength/2)];
% firsttrigger(10000:(10000+pulselength-1)) = 2^(nbits-1)/2 * [ones(1, pulselength/2) -1*ones(1,pulselength/2)];


firsttrigger = round( 2^(nbits-1)/2 * [ones(1, pulselength/2) zeros(1,pulselength/2) ...
                           ones(1, pulselength/2) zeros(1,pulselength/2) ...
                           ones(1, pulselength/2) zeros(1,pulselength/2) ...
                           zeros(1, triglength-3*pulselength)]);

pulselength
length(trigger)
length(firsttrigger)
pause


% Do a little error check since the triggers are so essential
if ( (length(trigger) ~= 32000) || (length(firsttrigger) ~= 32000) )
   error('The triggers must be separated by 32000 samples. No exceptions.');
end


%Writing Triggers Blocks of Length M
numblocks = floor(totlength / triglength);


%Opening Output File
fidout = fopen(outfile,'w');


if ( initbuffer ~= 0 )
   fwrite(fidout, zeros(1,initbuffer), dtype);
end

% we're going to use one tripple trigger at the beginning of the ripple
for k = 1:numblocks

   if ( k==1 ) % write out the tripple trigger

      fwrite(fidout, firsttrigger, dtype);
 
   else % write out the normal triggers

      fwrite(fidout, trigger, dtype);

   end

end % (for k)


% Writing the last block
if ( totlength - numblocks*triglength > pulselength )
   trigger = zeros(1, totlength - numblocks*triglength);
   trigger(1:pulselength)= 2^(nbits-1)/2 * ones(1, pulselength);
else
   trigger = 2^(nbits-1)/2 * ones(1, totlength - numblocks*triglength);
end

if ( ~isempty(trigger) )
   fwrite(fidout, trigger, dtype);
end


if ( outbuffer ~= 0 )
   fwrite(fidout, zeros(1,outbuffer), dtype);
end


% Closing Output File
fclose('all');


