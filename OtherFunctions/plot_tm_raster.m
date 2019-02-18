function plot_tm_raster(tmresp)
%plot_tm_raster - Raster of temporal modulation data
%
%   plot_tm_raster(tcresp)
%
%   tmresp : struct array holding temporal modulation response data.
%
%
%   caa 9/1/13

% Will need to check if spike times are in seconds or milliseconds

% Get some preliminary data.
tm = tmresp(1).tm; 
sm = tmresp(1).sm;

t = unique(tm);
t = t(:)';

s = unique(sm);
s = s(:)';

nreps = length(tm) / length(t); % use to find nreps, since length(s) == 1


% Labels for the y-axis
tmlabel = unique(round(100 * [t(1:5:length(t)) t(end)])/100);
tmtickmark = unique([1:nreps*5:length(t)*nreps length(t)*nreps]);

% For the time axis
ttick = [0:250:1500];
tlabel = [0:250:1500];


for n = 1:length(tmresp)

   exp = tmresp(n).exp;
   site = tmresp(n).site;
   chan = tmresp(n).chan;
   position = tmresp(n).position;
   stim = tmresp(n).stim;


   % plot the dot raster
   rast_figure = figure;

   for m = 1:length(s) % length(s) == 1 for Gideon's data

      spiketimes = [];

      % first get spike data for the raster at one amplitude level
      raster = cell(1, length(t)*nreps);
      vert_psth = zeros(1,length(t));

      for i = 1:length(t)
 
         i1 = find( tm==t(i) & sm==s(m) );

         for j = 1:length(i1)
            trial_spikes = tmresp(n).resp{i1(j)} * 1000;
            raster{(i-1)*nreps + j} = trial_spikes;
            spiketimes = [spiketimes trial_spikes(:)'];

            vert_psth(i) = vert_psth(i) + length(trial_spikes(trial_spikes<1000));
         end % (for j)

      end % (for i)


      subplot(2,3,[1 2]);
      hold on;

      % Gray box showing time when tone is 'ON'
      p = patch([1 1 1000 1000],[0 length(t)*nreps length(t)*nreps 0],[0.9 0.9 0.9]);
      set(p,'linestyle', 'none');

      % Plot spiketimes as dots, 'o'
      for i = 1:length(raster)
      % raster{i}
         hold on;
         plot(raster{i}, i * ones(size(raster{i})), 'ko', 'markersize', 1.25, 'markerfacecolor', 'k');  

         set(gca,'xtick', ttick, 'xticklabel', tlabel, 'fontsize', 7, 'fontname', 'arial');

         set(gca,'ytick', tmtickmark, 'yticklabel', tmlabel, 'fontname', 'arial', 'fontsize', 7);

      end % (for i)

      set(gca,'tickdir','out');
      axis([0 1500 0 length(t)*nreps]); % ISI was 500 ms
      xlabel('Time (ms)');
      ylabel('Temp Mod Freq (Hz)');


      subplot(2,3, [4 5]);
      edges = 0:10:1500;
      n = histc(spiketimes, edges);
      hb = bar(edges, n, 'histc');
      set(hb, 'facecolor', 0 * ones(1,3), 'edgecolor', 0*ones(1,3));
      xlim([0 1500]);
      tickpref;
      xlabel('Time (ms)');
      ylabel('Num of Spikes');



      subplot(2,3, [3]);
      hb = barh(vert_psth);
      set(hb, 'facecolor', 0 * ones(1,3), 'edgecolor', 0*ones(1,3));
      ylim([1 20]);
      tickpref;
      set(gca,'xticklabel', '', 'yticklabel', '');
%       xlabel('Time (ms)');
%       ylabel('Num of Spikes');

   end % (for m)

   suptitle( sprintf('cellInd = [%.0f %.0f %.0f %.0f]', cellInd) );
   orient tall;
   print_mfilename(mfilename);

end % (for n)


return;








