function [tmresp] = get_tm_data(thresh, tmparams)


base = 'C:\Users\craig\Data\Gideon\stimuli';
filename = 'gideon-dmr-500flo-30000fhi-4SM-100TM-40db-96khz-48DF-15min.spr';
envfile = fullfile(base, filename);

% Triggers:
% 1-300 : tuning curve triggers (300 stimuli)
% 301-500 : temporal modulation triggers (200 stimuli)
% 501-3200 : dmr triggers (2699 triggers)


% [tcresp] = gideon_freq_resp_area(stimAndSpikes, tcparams);

[tmresp] = gideon_temporal_modulation(stimAndSpikes, tmparams);

% strf = gideon_calculate_strf(stimAndSpikes, envfile);


return;




function [tmresp] = temporal_modulation_response(thresh, tmparams)
%gideon_temporal_modulation   Temporal modulation responses
%
%   [tmresp] = gideon_freq_resp_area(stimAndSpikes, tmparams)
%
%   stimAndSpikes : cell array of Gideon's data
%
%   tmparams : temporal modulation parameters from sound that was played
%
%   caa 7/31/13


fs = 1500;


% Triggers:
% Hope they are already included ...

deltatrig = trigger(2) - trigger(1);


tm = tmparams(:,3);
sm = tmparams(:,4);

tmresp = thresh;
fields = {'input_chan', 'return_chan', 'current'};
tmresp = rmfield(tmresp, fields);

for i = 1:length(thresh)

   cellInd = stimAndSpikes{i}.cellInd;
   spiketimes = stimAndSpikes{i}.spikeTimes;
   spiketimes = spiketimes(spiketimes < max(trigger)+2*deltatrig);

   resp = cell(length(tm), 1);  % make the cell array that holds the responses

   for ii = 1:length(trigger)

      if ( ii ~= length(trigger) )
         index = find( spiketimes >= trigger(ii) & spiketimes < trigger(ii+1) );
         if ( ~isempty(index) )
            resp{ii} = spiketimes(index) - trigger(ii);
         end
      else % for the last tone presentation only
         index = find( spiketimes >= trigger(ii) & spiketimes < trigger(ii)+deltatrig );
         if ( ~isempty(index) )
            resp{ii} = spiketimes(index) - trigger(ii);
         end
      end
   
   end % (for ii)


   tmresp(i).tmparams = tmparams;
   tmresp(i).tm = tm;
   tmresp(i).sm = sm;
   tmresp(i).resp = resp;

   fprintf('Finished %.0f of %.0f\n', i, length(thresh));

end % (for i) 

fprintf('Finished calculating temporal modulation responses.\n');

return;




function [tcresp] = gideon_freq_resp_area(stimAndSpikes, tcparams)
%gideon_freq_resp_area   Tuning curve responses
%
%   [tcresp] = gideon_freq_resp_area(stimAndSpikes, tcparams)
%
%   stimAndSpikes : cell array of Gideon's data
%
%   tcparams : tuning curve parameters from sound that was played
%
%   tcresp : struct array holding the tuning curve responses
%
%   caa 7/31/13



fs = 1500;


% Triggers:
% 1-300 : tuning curve triggers (300 stimuli)
% 301-500 : temporal modulation triggers (200 stimuli)
% 501-3200 : dmr triggers (2699 triggers)

% Get the trigger vector. Trigger times are in seconds
trigger = stimAndSpikes{1}.stimTimes;

% The first pulse in the initial triple trigger was missed
% remove the first two pulses, and insert the correct the 
% first trigger time
trigger = trigger(1:300);
deltatrig = trigger(2) - trigger(1);


freq = tcparams(:,4);
atten = tcparams(:,5);


for i = 2:length(stimAndSpikes)

   cellInd = stimAndSpikes{i}.cellInd;
   spiketimes = stimAndSpikes{i}.spikeTimes;
   spiketimes = spiketimes(spiketimes < max(trigger)+2*deltatrig);

   resp = cell(length(freq), 1);  % make the cell array that holds the responses

   for ii = 1:length(trigger)

      if ( ii < length(trigger) )
         index = find( spiketimes >= trigger(ii) & spiketimes < trigger(ii+1) );

         if ( ~isempty(index) )
            resp{ii} = spiketimes(index) - trigger(ii);
         end
      else % for the last tone presentation only
         index = find( spiketimes >= trigger(ii) & spiketimes < trigger(ii)+deltatrig );
         if ( ~isempty(index) )
            resp{ii} = spiketimes(index) - trigger(ii);
         end
      end
   
   end % (for ii)


   tcresp(i-1).cellInd = cellInd;
   tcresp(i-1).spiketimes = spiketimes;
   tcresp(i-1).trigger = trigger;
   tcresp(i-1).fs = fs;
   tcresp(i-1).tcparams = tcparams;
   tcresp(i-1).freq = freq;
   tcresp(i-1).atten = atten;
   tcresp(i-1).resp = resp;


   fprintf('Finished %.0f of %.0f\n', i, length(2:length(stimAndSpikes) ))

end % (for i) 

fprintf('Finished calculating tuning curves.\n');

return;



function strf = gideon_calculate_strf(stimAndSpikes, envfile, tbefore, tafter)
% 
% Simple script to call the strf function strfdbcalc.m
%
% strf = calculate_strf(spk, trigger, binaural, envfile, tbefore, tafter); 
%
% spk : struct array holding spike time data for the strf
%    calculations.
%
% trigger : array holding trigger times for the strf 
%    calculations.
%
% binaural : optional. Specifies if the stimulus was delivered
%            binaurally. 1 = yes, 0 = no. Default is 0.
%
% envfile : optional argument. Specifies the envelope file (*.spr)
%             to use for reverse correlation. The default is
%             #4, listed below.
%
% tbefore and tafter are optional arguments. They specify how much
% time you wish to process before a spike and how much time after.
% They must be in seconds. Example tafter = 0.050 and tbefore = 0.200,
% which means calculate an strf using stimulus data 50 ms before a
% spike occurs and 200 ms before a spike occurred. The default
% is tbefore = 0.200 and tafter = 0.050.
%
%   caa 7/31/13



if ( nargin < 5 )
   T1 = 0.050; % 50 ms 
   T2 = 0.200; % 200 ms
elseif ( nargin == 6 )
   T1 = tafter;
   T2 = tbefore;
else
   error('You need either 4 or 6 input args.');
end


i1 = findstr(envfile,'-');
i2 = findstr(envfile,'SM');
i3 = findstr(envfile,'TM');
i4 = findstr(envfile,'db');

sm = str2num(envfile(i2(end)-1));
tm = str2num(envfile(i3(end)-2:i3(end)-1));
mdb = str2num(envfile(i4(end)-2:i4(end)-1));

atten = 30;
spl = 70;
fs = 1500;

stype = 'dmr';

if ( ~isempty(findstr(stype,'dmr')) )
   snd = 'MR';
elseif ( ~isempty(findstr(stype,'rn')) )
   snd = 'RN';
else
   error('stim in struct array "spk" must be "dmr" or "rn".');
end

modtype = 'dB';
nblocks = 500;  % update display every 500 blocks


% Triggers:
% 1-300 : tuning curve triggers
% 301-500 : temporal modulation triggers
% 501-3200 : dmr triggers

% Get DMR triggers
trigger = stimAndSpikes{1}.stimTimes;

% The first pulse in the initial triple trigger was missed
% remove the first two pulses, and insert the correct the 
% first trigger time
trigger = trigger(501:end);
trigger = trigger(3:end);
trigger = [trigger(1)-500/1500 trigger];



% Triggers are in seconds. Convert to sample number
trigger = round(trigger * fs);  % sec * sample/sec
trigger = trigger(:)';


for i = 2:length(stimAndSpikes)

   cellInd = stimAndSpikes{i}.cellInd;
   spiketimes = stimAndSpikes{i}.spikeTimes;


   % Spiketimes are in seconds. Convert to sample number
   spet = round( spiketimes * fs ); % sec * sample/sec
   spet = spet(:)';

   [taxis, faxis, rf1, ~, pp, w01, ~, n01, ~, spln] = ...
      rtwstrfdb(envfile, T1, T2, spet, trigger, fs, spl, mdb, modtype, snd, nblocks);


   % Save data to struct array
   strf(i-1).cellInd = cellInd;
   strf(i-1).spiketimes = spiketimes;
   strf(i-1).trigger = trigger;
   strf(i-1).fs = fs;
   strf(i-1).envfile = envfile;
   strf(i-1).atten = atten;
   strf(i-1).sm = sm;
   strf(i-1).tm = tm;
   strf(i-1).mdb = mdb;
   strf(i-1).spl = spl;
   strf(i-1).taxis = taxis;
   strf(i-1).faxis = faxis;
   strf(i-1).n0contra = n01;
   strf(i-1).w0contra = w01;
   strf(i-1).rfcontra = rf1;
   strf(i-1).pp = pp;
   strf(i-1).spln = spln;

   fprintf('%.0f of %.0f completed.\n', i, length(stimAndSpikes) );
   pause(0.5);

end % (for)





