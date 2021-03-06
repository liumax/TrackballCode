function plot_ripple_transfer_function_params_topography(xproj, exp)
% plot_freq_resp_area_params_topography - plot the topography
% of frequency response area parameters, such as CF, threshold,
% latency, and Q10,20,30, and Q40.
%
% plot_freq_resp_area_params_topography(xproj)
%
% xproj : an nx2 matrix, with the first column holding the
%         site or penetration number, and the second column
%         holds the projection of the penetration along an
%         axis orthogonal to the line connecting the tips
%         of the posterior and anterior ectosylvian sulci.
%         Thus the second column is in mm and specifies the
%         distance along this ventral-dorsal axis where the
%         given site falls.
%
% exp : option string specifying the experiment to be plotted.
%       Default is no string. Example: exp = '2003-5-27'
%
% caa 9/16/03

if ( nargin == 0 )
   error('You need to input the dorsal-ventral projection results');
end


for i = 1:length(xproj(:,1))
   dsite(i).name = ['site' num2str(xproj(i,1))];
end

layer1pos = [0 200];
layer2pos = [200 400];
layer3pos = [400 800];
layer4pos = [800 1100];
layer5pos = [1100 1500];
layer6pos = [1500 2000];
whitematterpos = [2000 2400];

% The following layer variables will be cell arrays 
% with the following fields:
%
% (1) dorsal-ventral postion
% (2) best TMF
% (3) tMTF width 3dB
% (4) tMTF width 6dB
% (5) tMTF 3dB shape - 1 = bp, 0 = lp
% (6) tMTF 6dB shape - 1 = bp, 0 = lp
% (7) best SMF
% (8) sMTF width 3dB
% (9) sMTF width 6dB
% (10) sMTF 3dB shape - 1 = bp, 0 = lp
% (11) sMTF 6dB shape - 1 = bp, 0 = lp


% set up some variables for some very detailed histograms
% of the filter shapes for different mtfs

% variables to hold data over all layers
tbp3db = [];
tlp3db = [];
tbp6db = [];
tlp6db = [];
xbp3db = [];
xlp3db = [];
xbp6db = [];
xlp6db = [];

% variables to hold data, but separated into different layers
layer1tbp3db = [];
layer1tlp3db = [];
layer1tbp6db = [];
layer1tlp6db = [];
layer1xbp3db = [];
layer1xlp3db = [];
layer1xbp6db = [];
layer1xlp6db = [];

layer2tbp3db = [];
layer2tlp3db = [];
layer2tbp6db = [];
layer2tlp6db = [];
layer2xbp3db = [];
layer2xlp3db = [];
layer2xbp6db = [];
layer2xlp6db = [];

layer3tbp3db = [];
layer3tlp3db = [];
layer3tbp6db = [];
layer3tlp6db = [];
layer3xbp3db = [];
layer3xlp3db = [];
layer3xbp6db = [];
layer3xlp6db = [];

layer4tbp3db = [];
layer4tlp3db = [];
layer4tbp6db = [];
layer4tlp6db = [];
layer4xbp3db = [];
layer4xlp3db = [];
layer4xbp6db = [];
layer4xlp6db = [];

layer5tbp3db = [];
layer5tlp3db = [];
layer5tbp6db = [];
layer5tlp6db = [];
layer5xbp3db = [];
layer5xlp3db = [];
layer5xbp6db = [];
layer5xlp6db = [];

layer6tbp3db = [];
layer6tlp3db = [];
layer6tbp6db = [];
layer6tlp6db = [];
layer6xbp3db = [];
layer6xlp3db = [];
layer6xbp6db = [];
layer6xlp6db = [];


layer1 = [];
layer2 = [];
layer3 = [];
layer4 = [];
layer5 = [];
layer6 = [];
whitematter = [];

for i = 1:length(dsite)

   site = dsite(i).name;
   filename = [site '\*-' 'dmr1' '-*-rtf-params.mat'];
   dfile = dir(filename);

   if ( ~isempty(dfile) )

      infile = dfile.name;
      load([site '\' infile], 'rtf_params');

      dvpos = xproj(i,2);

      for j = 1:length(rtf_params)

         position = rtf_params(j).position;

         tmf = rtf_params(j).tmf;
         tmtf = rtf_params(j).tmtf ./ max(rtf_params(j).tmtf);

         xmf = rtf_params(j).xmf;
         xmtf = rtf_params(j).xmtf ./ max(rtf_params(j).xmtf);

         temp{1}{1} = dvpos;


         % ********************************************************
         % Process temporal modulation transfer function parameters
         % ********************************************************

         temp{1}{2} = rtf_params(j).best_tmf;

         tfit = min(tmf):0.001:max(tmf);
         mtf = interp1(tmf, tmtf, tfit, 'spline');

         indmax = find(mtf == max(mtf));

         ind3dblow = find(mtf(1:indmax) >= 0.7);
         ind3dbhi = find(mtf(indmax:end) >= 0.7);

         ind6dblow = find(mtf(1:indmax) >= 0.5);
         ind6dbhi = find(mtf(indmax:end) >= 0.5);
   
         temp{1}{3} = tfit(max(ind3dbhi+indmax-1)) - tfit(min(ind3dblow));
         temp{1}{4} = tfit(max(ind6dbhi+indmax-1)) - tfit(min(ind6dblow));

         % 3dB filter analysis
         if ( mtf(min(ind3dblow)) < 0.72 & mtf(max(ind3dbhi+indmax-1)) < 0.72 )

            temp{1}{5} = 1;

            if ( (position >= layer1pos(1)) & (position < layer6pos(2)) )
               tbp3db = [tbp3db dvpos];
            end

            if ( (position >= layer1pos(1)) & (position < layer1pos(2)) ) % layer 1
               layer1tbp3db = [layer1tbp3db dvpos];
            elseif ( (position >= layer2pos(1)) & (position < layer2pos(2)) ) % layer 2
               layer2tbp3db = [layer2tbp3db dvpos];
            elseif ( (position >= layer3pos(1)) & (position < layer3pos(2)) ) % layer 3
               layer3tbp3db = [layer3tbp3db dvpos];
            elseif ( (position >= layer4pos(1)) & (position < layer4pos(2)) ) % layer 4
               layer4tbp3db = [layer4tbp3db dvpos];
            elseif ( (position >= layer5pos(1)) & (position < layer5pos(2)) ) % layer 5
               layer5tbp3db = [layer5tbp3db dvpos];
            elseif ( (position >= layer6pos(1)) & (position < layer6pos(2)) ) % layer 6
               layer6tbp3db = [layer6tbp3db dvpos];
            elseif ( (position >= whitematterpos(1)) & (position < whitematterpos(2)) )
            end

         else

            temp{1}{5} = 0;

            if ( (position >= layer1pos(1)) & (position < layer6pos(2)) )
               tlp3db = [tlp3db dvpos];
            end

            if ( (position >= layer1pos(1)) & (position < layer1pos(2)) ) % layer 1
               layer1tlp3db = [layer1tlp3db dvpos];
            elseif ( (position >= layer2pos(1)) & (position < layer2pos(2)) ) % layer 2
               layer2tlp3db = [layer2tlp3db dvpos];
            elseif ( (position >= layer3pos(1)) & (position < layer3pos(2)) ) % layer 3
               layer3tlp3db = [layer3tlp3db dvpos];
            elseif ( (position >= layer4pos(1)) & (position < layer4pos(2)) ) % layer 4
               layer4tlp3db = [layer4tlp3db dvpos];
            elseif ( (position >= layer5pos(1)) & (position < layer5pos(2)) ) % layer 5
               layer5tlp3db = [layer5tlp3db dvpos];
            elseif ( (position >= layer6pos(1)) & (position < layer6pos(2)) ) % layer 6
               layer6tlp3db = [layer6tlp3db dvpos];
            elseif ( (position >= whitematterpos(1)) & (position < whitematterpos(2)) )
            end

         end


         % 6dB filter analysis
         if ( mtf(min(ind6dblow)) < 0.52 & mtf(max(ind6dbhi+indmax-1)) < 0.52 )

            temp{1}{6} = 1;

            if ( (position >= layer1pos(1)) & (position < layer6pos(2)) )
               tbp6db = [tbp6db dvpos];
            end

            if ( (position >= layer1pos(1)) & (position < layer1pos(2)) ) % layer 1
               layer1tbp6db = [layer1tbp6db dvpos];
            elseif ( (position >= layer2pos(1)) & (position < layer2pos(2)) ) % layer 2
               layer2tbp6db = [layer2tbp6db dvpos];
            elseif ( (position >= layer3pos(1)) & (position < layer3pos(2)) ) % layer 3
               layer3tbp6db = [layer3tbp6db dvpos];
            elseif ( (position >= layer4pos(1)) & (position < layer4pos(2)) ) % layer 4
               layer4tbp6db = [layer4tbp6db dvpos];
            elseif ( (position >= layer5pos(1)) & (position < layer5pos(2)) ) % layer 5
               layer5tbp6db = [layer5tbp6db dvpos];
            elseif ( (position >= layer6pos(1)) & (position < layer6pos(2)) ) % layer 6
               layer6tbp6db = [layer6tbp6db dvpos];
            elseif ( (position >= whitematterpos(1)) & (position < whitematterpos(2)) )
            end

         else

            temp{1}{6} = 0;

            if ( (position >= layer1pos(1)) & (position < layer6pos(2)) )
               tlp6db = [tlp6db dvpos];
            end

            if ( (position >= layer1pos(1)) & (position < layer1pos(2)) ) % layer 1
               layer1tlp6db = [layer1tlp6db dvpos];
            elseif ( (position >= layer2pos(1)) & (position < layer2pos(2)) ) % layer 2
               layer2tlp6db = [layer2tlp6db dvpos];
            elseif ( (position >= layer3pos(1)) & (position < layer3pos(2)) ) % layer 3
               layer3tlp6db = [layer3tlp6db dvpos];
            elseif ( (position >= layer4pos(1)) & (position < layer4pos(2)) ) % layer 4
               layer4tlp6db = [layer4tlp6db dvpos];
            elseif ( (position >= layer5pos(1)) & (position < layer5pos(2)) ) % layer 5
               layer5tlp6db = [layer5tlp6db dvpos];
            elseif ( (position >= layer6pos(1)) & (position < layer6pos(2)) ) % layer 6
               layer6tlp6db = [layer6tlp6db dvpos];
            elseif ( (position >= whitematterpos(1)) & (position < whitematterpos(2)) )
            end

         end

         % ********************************************************
         % Process spectral modulation transfer function parameters
         % ********************************************************
         temp{1}{7} = rtf_params(j).best_xmf;

         xfit = min(xmf):0.001:max(xmf);
         mtf = interp1(xmf, xmtf, xfit, 'spline');

         indmax = find(mtf == max(mtf));

         ind3dblow = find(mtf(1:indmax) >= 0.7);
         ind3dbhi = find(mtf(indmax:end) >= 0.7);

         ind6dblow = find(mtf(1:indmax) >= 0.5);
         ind6dbhi = find(mtf(indmax:end) >= 0.5);
   
         temp{1}{8} = xfit(max(ind3dbhi+indmax-1)) - xfit(min(ind3dblow));
         temp{1}{9} = xfit(max(ind6dbhi+indmax-1)) - xfit(min(ind6dblow));

         if ( mtf(min(ind3dblow)) < 0.72 & mtf(max(ind3dbhi+indmax-1)) < 0.72 )

            temp{1}{10} = 1;

            if ( (position >= layer1pos(1)) & (position < layer6pos(2)) )
               xbp3db = [xbp3db dvpos];
            end

            if ( (position >= layer1pos(1)) & (position < layer1pos(2)) ) % layer 1
               layer1xbp3db = [layer1xbp3db dvpos];
            elseif ( (position >= layer2pos(1)) & (position < layer2pos(2)) ) % layer 2
               layer2xbp3db = [layer2xbp3db dvpos];
            elseif ( (position >= layer3pos(1)) & (position < layer3pos(2)) ) % layer 3
               layer3xbp3db = [layer3xbp3db dvpos];
            elseif ( (position >= layer4pos(1)) & (position < layer4pos(2)) ) % layer 4
               layer4xbp3db = [layer4xbp3db dvpos];
            elseif ( (position >= layer5pos(1)) & (position < layer5pos(2)) ) % layer 5
               layer5xbp3db = [layer5xbp3db dvpos];
            elseif ( (position >= layer6pos(1)) & (position < layer6pos(2)) ) % layer 6
               layer6xbp3db = [layer6xbp3db dvpos];
            elseif ( (position >= whitematterpos(1)) & (position < whitematterpos(2)) )
            end

         else

            temp{1}{10} = 0;

            if ( (position >= layer1pos(1)) & (position < layer6pos(2)) )
               xlp3db = [xlp3db dvpos];
            end

            if ( (position >= layer1pos(1)) & (position < layer1pos(2)) ) % layer 1
               layer1xlp3db = [layer1xlp3db dvpos];
            elseif ( (position >= layer2pos(1)) & (position < layer2pos(2)) ) % layer 2
               layer2xlp3db = [layer2xlp3db dvpos];
            elseif ( (position >= layer3pos(1)) & (position < layer3pos(2)) ) % layer 3
               layer3xlp3db = [layer3xlp3db dvpos];
            elseif ( (position >= layer4pos(1)) & (position < layer4pos(2)) ) % layer 4
               layer4xlp3db = [layer4xlp3db dvpos];
            elseif ( (position >= layer5pos(1)) & (position < layer5pos(2)) ) % layer 5
               layer5xlp3db = [layer5xlp3db dvpos];
            elseif ( (position >= layer6pos(1)) & (position < layer6pos(2)) ) % layer 6
               layer6xlp3db = [layer6xlp3db dvpos];
            elseif ( (position >= whitematterpos(1)) & (position < whitematterpos(2)) )
            end

         end

         if ( mtf(min(ind6dblow)) < 0.52 & mtf(max(ind6dbhi+indmax-1)) < 0.52 )
 
            temp{1}{11} = 1;

            if ( (position >= layer1pos(1)) & (position < layer6pos(2)) )
               xbp6db = [xbp6db dvpos];
            end

            if ( (position >= layer1pos(1)) & (position < layer1pos(2)) ) % layer 1
               layer1xbp6db = [layer1xbp6db dvpos];
            elseif ( (position >= layer2pos(1)) & (position < layer2pos(2)) ) % layer 2
               layer2xbp6db = [layer2xbp6db dvpos];
            elseif ( (position >= layer3pos(1)) & (position < layer3pos(2)) ) % layer 3
               layer3xbp6db = [layer3xbp6db dvpos];
            elseif ( (position >= layer4pos(1)) & (position < layer4pos(2)) ) % layer 4
               layer4xbp6db = [layer4xbp6db dvpos];
            elseif ( (position >= layer5pos(1)) & (position < layer5pos(2)) ) % layer 5
               layer5xbp6db = [layer5xbp6db dvpos];
            elseif ( (position >= layer6pos(1)) & (position < layer6pos(2)) ) % layer 6
               layer6xbp6db = [layer6xbp6db dvpos];
            elseif ( (position >= whitematterpos(1)) & (position < whitematterpos(2)) )
            end

         else

            temp{1}{11} = 0;

            if ( (position >= layer1pos(1)) & (position < layer6pos(2)) )
               xlp6db = [xlp6db dvpos];
            end

            if ( (position >= layer1pos(1)) & (position < layer1pos(2)) ) % layer 1
               layer1xlp6db = [layer1xlp6db dvpos];
            elseif ( (position >= layer2pos(1)) & (position < layer2pos(2)) ) % layer 2
               layer2xlp6db = [layer2xlp6db dvpos];
            elseif ( (position >= layer3pos(1)) & (position < layer3pos(2)) ) % layer 3
               layer3xlp6db = [layer3xlp6db dvpos];
            elseif ( (position >= layer4pos(1)) & (position < layer4pos(2)) ) % layer 4
               layer4xlp6db = [layer4xlp6db dvpos];
            elseif ( (position >= layer5pos(1)) & (position < layer5pos(2)) ) % layer 5
               layer5xlp6db = [layer5xlp6db dvpos];
            elseif ( (position >= layer6pos(1)) & (position < layer6pos(2)) ) % layer 6
               layer6xlp6db = [layer6xlp6db dvpos];
            elseif ( (position >= whitematterpos(1)) & (position < whitematterpos(2)) )
            end

         end



         % ********************************************************
         % Save data in layer specific variables
         % ********************************************************
         if ( (position >= layer1pos(1)) & (position < layer1pos(2)) )
            layer1 = [layer1 temp];
            clear('temp');

         elseif ( (position >= layer2pos(1)) & (position < layer2pos(2)) )
            layer2 = [layer2 temp];
            clear('temp');

         elseif ( (position >= layer3pos(1)) & (position < layer3pos(2)) )
            layer3 = [layer3 temp];
            clear('temp');

         elseif ( (position >= layer4pos(1)) & (position < layer4pos(2)) )
            layer4 = [layer4 temp];
            clear('temp');

         elseif ( (position >= layer5pos(1)) & (position < layer5pos(2)) )
            layer5 = [layer5 temp];
            clear('temp');

         elseif ( (position >= layer6pos(1)) & (position < layer6pos(2)) )
            layer6 = [layer6 temp];
            clear('temp');

         elseif ( (position >= whitematterpos(1)) & (position < whitematterpos(2)) )
            whitematter = [whitematter temp];
            clear('temp');

         else
            clear('temp');
         end

      end % (for j)

   end % (if)

end % (for i)




% ***************** Temporal Modulation Analysis *****************

% The following layer variables will be cell arrays 
% with the following fields:
%
% (1) dorsal-ventral postion
% (2) best TMF
% (3) tMTF width 3dB
% (4) tMTF width 6dB
% (5) tMTF 3dB shape - 1 = bp, 0 = lp
% (6) tMTF 6dB shape - 1 = bp, 0 = lp
% (7) best SMF
% (8) sMTF width 3dB
% (9) sMTF width 6dB
% (10) sMTF 3dB shape - 1 = bp, 0 = lp
% (11) sMTF 6dB shape - 1 = bp, 0 = lp

figure;

% ------------------------------------------------------------
%              Plot the best TMF layer results
% ------------------------------------------------------------

btmflayerdata = cell(1,6);

for i = 1:6

   subplot(6, 4, 4*i-3);
   hold on;
   temp1 = [];
   temp2 = [];
   dv = [];
   btmf = [];

   for j = 1:eval( ['length(layer' num2str(i) ')'] )
      eval(['temp1 = layer' num2str(i) '{' num2str(j) '}{1};']);
      eval(['temp2 = layer' num2str(i) '{' num2str(j) '}{2};']);
      dv = [dv temp1];
      btmf = [btmf temp2];
   end % (for j)

   [dv, index] = sort(dv);
   btmf = btmf(index);

   plot(dv, btmf, 'ko', 'markerfacecolor', 'k', 'markersize', 3);

   btmflayerdata{i}{1} = [];
   btmflayerdata{i}{2} = [];
   btmflayerdata{i}{3} = [];
   btmflayerdata{i}{4} = [];
   btmflayerdata{i}{5} = [];
   btmflayerdata{i}{6} = [];
   btmflayerdata{i}{7} = [];

   if ( ~isempty(dv) )

      udv = unique(dv);
      ubtmfmn = zeros(size(udv));
      ubtmfmnsm = zeros(size(udv));
      ubtmfmx = zeros(size(udv));
      ubtmfmxsm = zeros(size(udv));
      ubtmfmi = zeros(size(udv));
      ubtmfmism = zeros(size(udv));

      for k = 1:length(udv)
         kk = find(dv == udv(k));
         ubtmfmn(k) = mean(btmf(kk));
         ubtmfmx(k) = max(btmf(kk));
         ubtmfmi(k) = min(btmf(kk));
      end

      ubtmfmnsm = smooth(ubtmfmn,3);
      ubtmfmxsm = smooth(ubtmfmx,3);
      ubtmfmism = smooth(ubtmfmi,3);

      btmflayerdata{i}{1} = udv;
      btmflayerdata{i}{2} = ubtmfmn;
      btmflayerdata{i}{3} = ubtmfmnsm;
      btmflayerdata{i}{4} = ubtmfmx;
      btmflayerdata{i}{5} = ubtmfmxsm;
      btmflayerdata{i}{6} = ubtmfmi;
      btmflayerdata{i}{7} = ubtmfmism;
   end

   set(gca,'ytick', 0:10:40, 'yticklabel', 0:10:40);
   set(gca,'fontsize', 8);
   axis([-1 6 -5 45]);

   if ( i == 1 )
      title('bTMF');
   end

   if ( i < 6 )
      ylabel(sprintf('Layer%.0f\n[%.0f, %.0f]\nn=%.0f', i, ...
         eval(['layer' num2str(i) 'pos(1)']), ...
         eval(['layer' num2str(i) 'pos(2)']), ...
         eval(['length(layer' num2str(i) ')'])));
      set(gca, 'xticklabel', '');
   else
      xlabel('Ventral < -- > Dorsal (mm)');
      ylabel(sprintf('Layer%.0f\n[%.0f, %.0f]\nbTMF (Hz)\nn=%.0f', i, ...
         eval(['layer' num2str(i) 'pos(1)']), ...
         eval(['layer' num2str(i) 'pos(2)']), ...
         eval(['length(layer' num2str(i) ')'])));
   end

end % (for i)


% ------------------------------------------------------------
%           Plot the tMTF 3dB width layer results
% ------------------------------------------------------------

tw3dblayerdata = cell(1,6);

for i = 1:6

   subplot(6, 4, 4*i-2);
   hold on;
   temp1 = [];
   temp2 = [];
   temp3 = [];
   dv = [];
   tw3db = [];
   shape = [];

   for j = 1:eval( ['length(layer' num2str(i) ')'] )
      eval(['temp1 = layer' num2str(i) '{' num2str(j) '}{1};']);
      eval(['temp2 = layer' num2str(i) '{' num2str(j) '}{3};']);
      eval(['temp3 = layer' num2str(i) '{' num2str(j) '}{5};']);
      dv = [dv temp1];
      tw3db = [tw3db temp2];
      shape = [shape temp3];
   end % (for j)


   [dv, index] = sort(dv);
   tw3db = tw3db(index);
   shape = shape(index);

   for jj = 1:length(dv)
      if ( shape(jj) ) % it's bandpass
         plot(dv(jj), tw3db(jj), 'ro', 'markerfacecolor', 'r', 'markersize', 3);
      else
         plot(dv(jj), tw3db(jj), 'ko', 'markerfacecolor', 'k', 'markersize', 3);
      end
   end % (for jj)

   tw3dblayerdata{i}{1} = [];
   tw3dblayerdata{i}{2} = [];
   tw3dblayerdata{i}{3} = [];
   tw3dblayerdata{i}{4} = [];
   tw3dblayerdata{i}{5} = [];
   tw3dblayerdata{i}{6} = [];
   tw3dblayerdata{i}{7} = [];

   if ( ~isempty(dv) )

      udv = unique(dv);
      utw3dbmn = zeros(size(udv));
      utw3dbmnsm = zeros(size(udv));
      utw3dbmx = zeros(size(udv));
      utw3dbmxsm = zeros(size(udv));
      utw3dbmi = zeros(size(udv));
      utw3dbmism = zeros(size(udv));

      for k = 1:length(udv)
         kk = find(dv == udv(k));
         utw3dbmn(k) = mean(tw3db(kk));
         utw3dbmx(k) = max(tw3db(kk));
         utw3dbmi(k) = min(tw3db(kk));
      end

      utw3dbmnsm = smooth(utw3dbmn,3);
      utw3dbmxsm = smooth(utw3dbmx,3);
      utw3dbmism = smooth(utw3dbmi,3);

      tw3dblayerdata{i}{1} = udv;
      tw3dblayerdata{i}{2} = utw3dbmn;
      tw3dblayerdata{i}{3} = utw3dbmnsm;
      tw3dblayerdata{i}{4} = utw3dbmx;
      tw3dblayerdata{i}{5} = utw3dbmxsm;
      tw3dblayerdata{i}{6} = utw3dbmi;
      tw3dblayerdata{i}{7} = utw3dbmism;
   end

   set(gca,'ytick', 0:10:40, 'yticklabel', 0:10:40, 'tickdir', 'out');
   set(gca,'fontsize', 8);
   axis([-1 6 -5 45]);

   if ( i == 1 )
      plot(-10, -10, 'ro', 'markerfacecolor', 'r', 'markersize', 3);
      plot(-10, -10, 'ko', 'markerfacecolor', 'k', 'markersize', 3);
      title('Width - 3dB');
      legend('BP','LP',0);
   end

   if ( i == 6 )
      ylabel(sprintf('Width (Hz)'));
   else
      set(gca,'xticklabel', '');
   end

end % (for i)



% ------------------------------------------------------------
%           Plot the tMTF 6dB width layer results
% ------------------------------------------------------------

tw6dblayerdata = cell(1,6);

for i = 1:6

   subplot(6, 4, 4*i-1);
   hold on;
   temp1 = [];
   temp2 = [];
   temp3 = [];
   dv = [];
   tw6db = [];
   shape = [];

   for j = 1:eval( ['length(layer' num2str(i) ')'] )
      eval(['temp1 = layer' num2str(i) '{' num2str(j) '}{1};']);
      eval(['temp2 = layer' num2str(i) '{' num2str(j) '}{4};']);
      eval(['temp3 = layer' num2str(i) '{' num2str(j) '}{6};']);
      dv = [dv temp1];
      tw6db = [tw6db temp2];
      shape = [shape temp3];
   end % (for j)

   [dv, index] = sort(dv);
   tw6db = tw6db(index);
   shape = shape(index);

   for jj = 1:length(dv)
      if ( shape(jj) ) % it's bandpass
         plot(dv(jj), tw6db(jj), 'ro', 'markerfacecolor', 'r', 'markersize', 3);
      else
         plot(dv(jj), tw6db(jj), 'ko', 'markerfacecolor', 'k', 'markersize', 3);
      end
   end % (for jj)

   tw6dblayerdata{i}{1} = [];
   tw6dblayerdata{i}{2} = [];
   tw6dblayerdata{i}{3} = [];
   tw6dblayerdata{i}{4} = [];
   tw6dblayerdata{i}{5} = [];
   tw6dblayerdata{i}{6} = [];
   tw6dblayerdata{i}{7} = [];

   if ( ~isempty(dv) )

      udv = unique(dv);
      utw6dbmn = zeros(size(udv));
      utw6dbmnsm = zeros(size(udv));
      utw6dbmx = zeros(size(udv));
      utw6dbmxsm = zeros(size(udv));
      utw6dbmi = zeros(size(udv));
      utw6dbmism = zeros(size(udv));

      for k = 1:length(udv)
         kk = find(dv == udv(k));
         utw6dbmn(k) = mean(tw6db(kk));
         utw6dbmx(k) = max(tw6db(kk));
         utw6dbmi(k) = min(tw6db(kk));
      end

      utw6dbmnsm = smooth(utw6dbmn,3);
      utw6dbmxsm = smooth(utw6dbmx,3);
      utw6dbmism = smooth(utw6dbmi,3);

      tw6dblayerdata{i}{1} = udv;
      tw6dblayerdata{i}{2} = utw6dbmn;
      tw6dblayerdata{i}{3} = utw6dbmnsm;
      tw6dblayerdata{i}{4} = utw6dbmx;
      tw6dblayerdata{i}{5} = utw6dbmxsm;
      tw6dblayerdata{i}{6} = utw6dbmi;
      tw6dblayerdata{i}{7} = utw6dbmism;
   end

   set(gca,'ytick', 0:10:40, 'yticklabel', 0:10:40, 'tickdir', 'out');
   set(gca,'fontsize', 8);
   axis([-1 6 -5 45]);

   if ( i == 1 )
      title('Width - 6dB');
   end

   if ( i == 6 )
      ylabel(sprintf('Width (Hz)'));
   else
      set(gca,'xticklabel', '');
   end

end % (for i)


% ------------------------------------------------------------
%           Plot the filter shape layer results
% ------------------------------------------------------------

dvposbins = -1:0.5:6;

for i = 1:6

   subplot(6, 4, 4*i);
   hold on;

   if ( i == 1 ) % 3dB cutoff filter shapes histogram

      tnbp = histc([tbp3db -100], dvposbins);
      tnlp = histc([tlp3db -100], dvposbins);
      bar(dvposbins, [tnlp(:) tnbp(:)], 'stacked');
      temp = tnlp + tnbp;
      axis([-1 6 0 max([round(1.1*max(temp)) max(temp)+1])]);
      set(gca,'fontsize', 8)
      children = get(gca,'children');
      set(children(1),'facecolor', [0.75 0.75 0.75]);
      set(children(2),'facecolor', [0 0 0]);
      legend('LP', 'BP');
      title('Filter Shape - 3dB Cutoff');
      set(gca,'xticklabel', '', 'tickdir', 'out');

   elseif ( i == 2 ) % 6dB cutoff filter shapes histogram

      tnbp = histc([tbp6db -100], dvposbins);
      tnlp = histc([tlp6db -100], dvposbins);
      bar(dvposbins, [tnlp(:) tnbp(:)], 'stacked');
      temp = tnbp + tnlp;
      axis([-1 6 0 max([round(1.1*max(temp)) max(temp)+1])]);
      set(gca,'fontsize', 8)
      children = get(gca,'children');
      set(children(1),'facecolor', [0.75 0.75 0.75]);
      set(children(2),'facecolor', [0 0 0]);
      title('6dB Cutoff');
      set(gca,'xticklabel', '', 'tickdir', 'out');

   elseif ( i == 3 ) % 3dB cutoff bandpass filter shape histogram

      layer1tnbp = histc([layer1tbp3db -100], dvposbins);
      layer2tnbp = histc([layer2tbp3db -100], dvposbins);
      layer3tnbp = histc([layer3tbp3db -100], dvposbins);
      layer4tnbp = histc([layer4tbp3db -100], dvposbins);
      layer5tnbp = histc([layer5tbp3db -100], dvposbins);
      layer6tnbp = histc([layer6tbp3db -100], dvposbins);
      temp = layer1tnbp + layer2tnbp + layer3tnbp + ...
         layer4tnbp + layer5tnbp + layer6tnbp;
      bar(dvposbins, [layer6tnbp(:) layer5tnbp(:) layer4tnbp(:) ...
         layer3tnbp(:) layer2tnbp(:) layer1tnbp(:)], 'stacked');
      axis([-1 6 0 max([round(1.1*max(temp)) max(temp)+1])]);
      set(gca,'fontsize', 8);
      children = get(gca,'children');

      set(children(1),'facecolor', [1 1 1]); % layer 1
      set(children(2),'facecolor', [.85 .85 .85]); % layer 2
      set(children(3),'facecolor', [.7 .7 .7]); % layer 3
      set(children(4),'facecolor', [.55 .55 .55]); % layer 4
      set(children(5),'facecolor', [.4 .4 .4]); % layer 5
      set(children(6),'facecolor', [0 0 0]); % layer 6
      lh = legend('L1', 'L2', 'L3', 'L4', 'L5', 'L6',0);
      set(lh, 'ydir', 'reverse');
      title('3dB BandPass');
      set(gca,'xticklabel', '', 'tickdir', 'out');

   elseif ( i == 4 ) % 6dB cutoff bandpass filter shape histogram

      layer1tnbp = histc([layer1tbp6db -100], dvposbins);
      layer2tnbp = histc([layer2tbp6db -100], dvposbins);
      layer3tnbp = histc([layer3tbp6db -100], dvposbins);
      layer4tnbp = histc([layer4tbp6db -100], dvposbins);
      layer5tnbp = histc([layer5tbp6db -100], dvposbins);
      layer6tnbp = histc([layer6tbp6db -100], dvposbins);
      temp = layer1tnbp + layer2tnbp + layer3tnbp + ...
         layer4tnbp + layer5tnbp + layer6tnbp;
      bar(dvposbins, [layer6tnbp(:) layer5tnbp(:) layer4tnbp(:) ...
         layer3tnbp(:) layer2tnbp(:) layer1tnbp(:)], 'stacked');
      axis([-1 6 0 max([round(1.1*max(temp)) max(temp)+1])]);
      set(gca,'fontsize', 8);
      children = get(gca,'children');

      set(children(1),'facecolor', [1 1 1]); % layer 1
      set(children(2),'facecolor', [.85 .85 .85]); % layer 2
      set(children(3),'facecolor', [.7 .7 .7]); % layer 3
      set(children(4),'facecolor', [.55 .55 .55]); % layer 4
      set(children(5),'facecolor', [.4 .4 .4]); % layer 5
      set(children(6),'facecolor', [0 0 0]); % layer 6
      title('6dB BandPass');
      set(gca,'xticklabel', '', 'tickdir', 'out');

   elseif ( i == 5 ) % 3dB cutoff lowpass filter shape histogram

      layer1tnbp = histc([layer1tlp3db -100], dvposbins);
      layer2tnbp = histc([layer2tlp3db -100], dvposbins);
      layer3tnbp = histc([layer3tlp3db -100], dvposbins);
      layer4tnbp = histc([layer4tlp3db -100], dvposbins);
      layer5tnbp = histc([layer5tlp3db -100], dvposbins);
      layer6tnbp = histc([layer6tlp3db -100], dvposbins);
      temp = layer1tnbp + layer2tnbp + layer3tnbp + ...
         layer4tnbp + layer5tnbp + layer6tnbp;
      bar(dvposbins, [layer6tnbp(:) layer5tnbp(:) layer4tnbp(:) ...
         layer3tnbp(:) layer2tnbp(:) layer1tnbp(:)], 'stacked');
      axis([-1 6 0 max([round(1.1*max(temp)) max(temp)+1])]);
      set(gca,'fontsize', 8);
      children = get(gca,'children');

      set(children(1),'facecolor', [1 1 1]); % layer 1
      set(children(2),'facecolor', [.85 .85 .85]); % layer 2
      set(children(3),'facecolor', [.7 .7 .7]); % layer 3
      set(children(4),'facecolor', [.55 .55 .55]); % layer 4
      set(children(5),'facecolor', [.4 .4 .4]); % layer 5
      set(children(6),'facecolor', [0 0 0]); % layer 6
      title('3dB LowPass');
      set(gca,'xticklabel', '', 'tickdir', 'out');

   elseif ( i == 6 ) % 6dB cutoff lowpass filter shape histogram

      layer1tnbp = histc([layer1tlp6db -100], dvposbins);
      layer2tnbp = histc([layer2tlp6db -100], dvposbins);
      layer3tnbp = histc([layer3tlp6db -100], dvposbins);
      layer4tnbp = histc([layer4tlp6db -100], dvposbins);
      layer5tnbp = histc([layer5tlp6db -100], dvposbins);
      layer6tnbp = histc([layer6tlp6db -100], dvposbins);
      temp = layer1tnbp + layer2tnbp + layer3tnbp + ...
         layer4tnbp + layer5tnbp + layer6tnbp;
      bar(dvposbins, [layer6tnbp(:) layer5tnbp(:) layer4tnbp(:) ...
         layer3tnbp(:) layer2tnbp(:) layer1tnbp(:)], 'stacked');
      axis([-1 6 0 max([round(1.1*max(temp)) max(temp)+1])]);
      set(gca,'fontsize', 8);
      children = get(gca,'children');

      set(children(1),'facecolor', [1 1 1]); % layer 1
      set(children(2),'facecolor', [.85 .85 .85]); % layer 2
      set(children(3),'facecolor', [.7 .7 .7]); % layer 3
      set(children(4),'facecolor', [.55 .55 .55]); % layer 4
      set(children(5),'facecolor', [.4 .4 .4]); % layer 5
      set(children(6),'facecolor', [0 0 0]); % layer 6
      title('6dB LowPass');
      set(gca, 'tickdir', 'out');

   end

end % (for i)


if ( nargin == 2 )
   if ( isstr(exp) )
      suptitle([exp '  Temporal Modulation Parameter Topography']);
   end
end

orient landscape;


% ***************** Spectral Modulation Analysis *****************

% The following layer variables will be cell arrays 
% with the following fields:
%
% (1) dorsal-ventral postion
% (2) best TMF
% (3) tMTF width 3dB
% (4) tMTF width 6dB
% (5) tMTF 3dB shape - 1 = bp, 0 = lp
% (6) tMTF 6dB shape - 1 = bp, 0 = lp
% (7) best SMF
% (8) sMTF width 3dB
% (9) sMTF width 6dB
% (10) sMTF 3dB shape - 1 = bp, 0 = lp
% (11) sMTF 6dB shape - 1 = bp, 0 = lp

figure;

% ------------------------------------------------------------
%              Plot the best SMF layer results
% ------------------------------------------------------------

bsmflayerdata = cell(1,6);

for i = 1:6

   subplot(6, 4, 4*i-3);
   hold on;
   temp1 = [];
   temp2 = [];
   dv = [];
   bsmf = [];

   for j = 1:eval( ['length(layer' num2str(i) ')'] )
      eval(['temp1 = layer' num2str(i) '{' num2str(j) '}{1};']);
      eval(['temp2 = layer' num2str(i) '{' num2str(j) '}{7};']);
      dv = [dv temp1];
      bsmf = [bsmf temp2];
   end % (for j)

   [dv, index] = sort(dv);
   bsmf = bsmf(index);

   plot(dv, bsmf, 'ko', 'markerfacecolor', 'k', 'markersize', 3);

   bsmflayerdata{i}{1} = [];
   bsmflayerdata{i}{2} = [];
   bsmflayerdata{i}{3} = [];
   bsmflayerdata{i}{4} = [];
   bsmflayerdata{i}{5} = [];
   bsmflayerdata{i}{6} = [];
   bsmflayerdata{i}{7} = [];

   if ( ~isempty(dv) )

      udv = unique(dv);
      ubsmfmn = zeros(size(udv));
      ubsmfmnsm = zeros(size(udv));
      ubsmfmx = zeros(size(udv));
      ubsmfmxsm = zeros(size(udv));
      ubsmfmi = zeros(size(udv));
      ubsmfmism = zeros(size(udv));

      for k = 1:length(udv)
         kk = find(dv == udv(k));
         ubsmfmn(k) = mean(bsmf(kk));
         ubsmfmx(k) = max(bsmf(kk));
         ubsmfmi(k) = min(bsmf(kk));
      end

      ubsmfmnsm = smooth(ubsmfmn,3);
      ubsmfmxsm = smooth(ubsmfmx,3);
      ubsmfmism = smooth(ubsmfmi,3);

      bsmflayerdata{i}{1} = udv;
      bsmflayerdata{i}{2} = ubsmfmn;
      bsmflayerdata{i}{3} = ubsmfmnsm;
      bsmflayerdata{i}{4} = ubsmfmx;
      bsmflayerdata{i}{5} = ubsmfmxsm;
      bsmflayerdata{i}{6} = ubsmfmi;
      bsmflayerdata{i}{7} = ubsmfmism;
   end

   set(gca,'ytick', 0:1:4, 'yticklabel', 0:1:4, 'tickdir', 'out');
   set(gca,'fontsize', 8);
   axis([-1 6 -1 5]);

   if ( i == 1 )
      title('bSMF');
   end

   if ( i < 6 )
      ylabel(sprintf('Layer%.0f\n[%.0f, %.0f]\nn=%.0f', i, ...
         eval(['layer' num2str(i) 'pos(1)']), ...
         eval(['layer' num2str(i) 'pos(2)']), ...
         eval(['length(layer' num2str(i) ')'])));
      set(gca, 'xticklabel', '');
   else
      xlabel('Ventral < -- > Dorsal (mm)');
      ylabel(sprintf('Layer%.0f\n[%.0f, %.0f]\nbTMF (Hz)\nn=%.0f', i, ...
         eval(['layer' num2str(i) 'pos(1)']), ...
         eval(['layer' num2str(i) 'pos(2)']), ...
         eval(['length(layer' num2str(i) ')'])));
   end

end % (for i)



% ------------------------------------------------------------
%            Plot the sMTF 3dB width layer results
% ------------------------------------------------------------

sw3dblayerdata = cell(1,6);

for i = 1:6

   subplot(6, 4, 4*i-2);
   hold on;
   temp1 = [];
   temp2 = [];
   temp3 = [];
   dv = [];
   sw3db = [];
   shape = [];

   for j = 1:eval( ['length(layer' num2str(i) ')'] )
      eval(['temp1 = layer' num2str(i) '{' num2str(j) '}{1};']);
      eval(['temp2 = layer' num2str(i) '{' num2str(j) '}{8};']);
      eval(['temp3 = layer' num2str(i) '{' num2str(j) '}{10};']);
      dv = [dv temp1];
      sw3db = [sw3db temp2];
      shape = [shape temp3];
   end % (for j)

   [dv, index] = sort(dv);
   sw3db = sw3db(index);
   shape = shape(index);

   for jj = 1:length(dv)
      if ( shape(jj) ) % it's bandpass
         plot(dv(jj), sw3db(jj), 'ro', 'markerfacecolor', 'r', 'markersize', 3);
      else  % it's lowpass
         plot(dv(jj), sw3db(jj), 'ko', 'markerfacecolor', 'k', 'markersize', 3);
      end
   end % (for jj)

   sw3dblayerdata{i}{1} = [];
   sw3dblayerdata{i}{2} = [];
   sw3dblayerdata{i}{3} = [];
   sw3dblayerdata{i}{4} = [];
   sw3dblayerdata{i}{5} = [];
   sw3dblayerdata{i}{6} = [];
   sw3dblayerdata{i}{7} = [];

   if ( ~isempty(dv) )

      udv = unique(dv);
      usw3dbmn = zeros(size(udv));
      usw3dbmnsm = zeros(size(udv));
      usw3dbmx = zeros(size(udv));
      usw3dbmxsm = zeros(size(udv));
      usw3dbmi = zeros(size(udv));
      usw3dbmism = zeros(size(udv));

      for k = 1:length(udv)
         kk = find(dv == udv(k));
         usw3dbmn(k) = mean(sw3db(kk));
         usw3dbmx(k) = max(sw3db(kk));
         usw3dbmi(k) = min(sw3db(kk));
      end

      usw3dbmnsm = smooth(usw3dbmn,3);
      usw3dbmxsm = smooth(usw3dbmx,3);
      usw3dbmism = smooth(usw3dbmi,3);

      sw3dblayerdata{i}{1} = udv;
      sw3dblayerdata{i}{2} = usw3dbmn;
      sw3dblayerdata{i}{3} = usw3dbmnsm;
      sw3dblayerdata{i}{4} = usw3dbmx;
      sw3dblayerdata{i}{5} = usw3dbmxsm;
      sw3dblayerdata{i}{6} = usw3dbmi;
      sw3dblayerdata{i}{7} = usw3dbmism;
   end

   set(gca,'ytick', 0:1:4, 'yticklabel', 0:1:4, 'tickdir', 'out');
   set(gca,'fontsize', 8);
   axis([-1 6 -1 5]);

   if ( i == 1 )
      title('Width - 3dB');
      plot(-10, -10, 'ro', 'markerfacecolor', 'r', 'markersize', 3);
      plot(-10, -10, 'ko', 'markerfacecolor', 'r', 'markersize', 3);
      lh = legend('BP','LP',0);
      children = get(lh,'children');
      set(children(2),'tag', 'BP'); 
      set(children(4),'tag', 'LP'); 
   end

   if ( i == 6 )
      ylabel(sprintf('Width (cyc/oct)'));
   else
      set(gca,'xticklabel', '');
   end

end % (for i)



% ------------------------------------------------------------
%           Plot the sMTF 6dB width layer results
% ------------------------------------------------------------

sw6dblayerdata = cell(1,6);

for i = 1:6

   subplot(6, 4, 4*i-1);
   hold on;
   temp1 = [];
   temp2 = [];
   temp3 = [];
   dv = [];
   sw6db = [];
   shape = [];

   for j = 1:eval( ['length(layer' num2str(i) ')'] )
      eval(['temp1 = layer' num2str(i) '{' num2str(j) '}{1};']);
      eval(['temp2 = layer' num2str(i) '{' num2str(j) '}{9};']);
      eval(['temp3 = layer' num2str(i) '{' num2str(j) '}{11};']);
      dv = [dv temp1];
      sw6db = [sw6db temp2];
      shape = [shape temp3];
   end % (for j)

   [dv, index] = sort(dv);
   sw6db = sw6db(index);
   shape = shape(index);

   for jj = 1:length(dv)
      if ( shape(jj) ) % it's bandpass
         plot(dv(jj), sw6db(jj), 'ro', 'markerfacecolor', 'r', 'markersize', 3);
      else
         plot(dv(jj), sw6db(jj), 'ko', 'markerfacecolor', 'k', 'markersize', 3);
      end
   end % (for jj)

   sw6dblayerdata{i}{1} = [];
   sw6dblayerdata{i}{2} = [];
   sw6dblayerdata{i}{3} = [];
   sw6dblayerdata{i}{4} = [];
   sw6dblayerdata{i}{5} = [];
   sw6dblayerdata{i}{6} = [];
   sw6dblayerdata{i}{7} = [];

   if ( ~isempty(dv) )

      udv = unique(dv);
      usw6dbmn = zeros(size(udv));
      usw6dbmnsm = zeros(size(udv));
      usw6dbmx = zeros(size(udv));
      usw6dbmxsm = zeros(size(udv));
      usw6dbmi = zeros(size(udv));
      usw6dbmism = zeros(size(udv));

      for k = 1:length(udv)
         kk = find(dv == udv(k));
         usw6dbmn(k) = mean(sw6db(kk));
         usw6dbmx(k) = max(sw6db(kk));
         usw6dbmi(k) = min(sw6db(kk));
      end

      usw6dbmnsm = smooth(usw6dbmn,3);
      usw6dbmxsm = smooth(usw6dbmx,3);
      usw6dbmism = smooth(usw6dbmi,3);

      sw6dblayerdata{i}{1} = udv;
      sw6dblayerdata{i}{2} = usw6dbmn;
      sw6dblayerdata{i}{3} = usw6dbmnsm;
      sw6dblayerdata{i}{4} = usw6dbmx;
      sw6dblayerdata{i}{5} = usw6dbmxsm;
      sw6dblayerdata{i}{6} = usw6dbmi;
      sw6dblayerdata{i}{7} = usw6dbmism;
   end

   set(gca,'ytick', 0:1:4, 'yticklabel', 0:1:4, 'tickdir', 'out');
   set(gca,'fontsize', 8);
   axis([-1 6 -1 5]);

   if ( i == 1 )
      title('Width - 6dB');
   end

   if ( i == 6 )
      ylabel(sprintf('Width (cyc/oct)'));
   else
      set(gca,'xticklabel', '');
   end

end % (for i)


% Plot the filter shape layer results

dvposbins = -1:0.5:6;

for i = 1:6

   subplot(6, 4, 4*i);
   hold on;

   if ( i == 1 )

      xnbp = histc(xbp3db, dvposbins);
      xnlp = histc(xlp3db, dvposbins);
      bar(dvposbins, [xnlp(:) xnbp(:)], 'stacked');
      temp = xnlp + xnbp;
      axis([-1 6 0 max([round(1.1*max(temp)) max(temp)+1])]);
      set(gca,'fontsize', 8)
      children = get(gca,'children');
      set(children(1),'facecolor', [0.75 0.75 0.75]);
      set(children(2),'facecolor', [0 0 0]);
      legend('LP', 'BP');
      title('Filter Shape - 3dB Cutoff');
      set(gca,'xticklabel', '', 'tickdir', 'out');

   elseif ( i == 2 )

      xnbp = histc(xbp6db, dvposbins);
      xnlp = histc(xlp6db, dvposbins);
      bar(dvposbins, [xnlp(:) xnbp(:)], 'stacked');
      temp = xnbp + xnlp;
      axis([-1 6 0 max([round(1.1*max(temp)) max(temp)+1])]);
      set(gca,'fontsize', 8)
      children = get(gca,'children');
      set(children(1),'facecolor', [0.75 0.75 0.75]);
      set(children(2),'facecolor', [0 0 0]);
      title('6dB Cutoff');
      set(gca,'xticklabel', '', 'tickdir', 'out');

   elseif ( i == 3 )

      layer1xnbp = histc([layer1xbp3db -100], dvposbins);
      layer2xnbp = histc([layer2xbp3db -100], dvposbins);
      layer3xnbp = histc([layer3xbp3db -100], dvposbins);
      layer4xnbp = histc([layer4xbp3db -100], dvposbins);
      layer5xnbp = histc([layer5xbp3db -100], dvposbins);
      layer6xnbp = histc([layer6xbp3db -100], dvposbins);
      temp = layer1xnbp + layer2xnbp + layer3xnbp + layer4xnbp + layer5xnbp + layer6xnbp;
      bar(dvposbins, [layer6xnbp(:) layer5xnbp(:) layer4xnbp(:) layer3xnbp(:) layer2xnbp(:) layer1xnbp(:)], 'stacked');
      axis([-1 6 0 max([round(1.1*max(temp)) max(temp)+1])]);
      set(gca,'fontsize', 8);
      children = get(gca,'children');

      set(children(1),'facecolor', [1 1 1]); % layer 1
      set(children(2),'facecolor', [.85 .85 .85]); % layer 2
      set(children(3),'facecolor', [.7 .7 .7]); % layer 3
      set(children(4),'facecolor', [.55 .55 .55]); % layer 4
      set(children(5),'facecolor', [.4 .4 .4]); % layer 5
      set(children(6),'facecolor', [0 0 0]); % layer 6
      lh = legend('L1', 'L2', 'L3', 'L4', 'L5', 'L6',0);
      set(lh, 'ydir', 'reverse');
      title('3dB BandPass');
      set(gca,'xticklabel', '', 'tickdir', 'out');

   elseif ( i == 4 )

      layer1xnbp = histc([layer1xbp6db -100], dvposbins);
      layer2xnbp = histc([layer2xbp6db -100], dvposbins);
      layer3xnbp = histc([layer3xbp6db -100], dvposbins);
      layer4xnbp = histc([layer4xbp6db -100], dvposbins);
      layer5xnbp = histc([layer5xbp6db -100], dvposbins);
      layer6xnbp = histc([layer6xbp6db -100], dvposbins);
      temp = layer1xnbp + layer2xnbp + layer3xnbp + layer4xnbp + layer5xnbp + layer6xnbp;
      bar(dvposbins, [layer6xnbp(:) layer5xnbp(:) layer4xnbp(:) layer3xnbp(:) layer2xnbp(:) layer1xnbp(:)], 'stacked');
      axis([-1 6 0 max([round(1.1*max(temp)) max(temp)+1])]);
      set(gca,'fontsize', 8);
      children = get(gca,'children');

      set(children(1),'facecolor', [1 1 1]); % layer 1
      set(children(2),'facecolor', [.85 .85 .85]); % layer 2
      set(children(3),'facecolor', [.7 .7 .7]); % layer 3
      set(children(4),'facecolor', [.55 .55 .55]); % layer 4
      set(children(5),'facecolor', [.4 .4 .4]); % layer 5
      set(children(6),'facecolor', [0 0 0]); % layer 6
      title('6dB BandPass');
      set(gca,'xticklabel', '', 'tickdir', 'out');

   elseif ( i == 5 )

      layer1xnbp = histc([layer1xlp3db -100], dvposbins);
      layer2xnbp = histc([layer2xlp3db -100], dvposbins);
      layer3xnbp = histc([layer3xlp3db -100], dvposbins);
      layer4xnbp = histc([layer4xlp3db -100], dvposbins);
      layer5xnbp = histc([layer5xlp3db -100], dvposbins);
      layer6xnbp = histc([layer6xlp3db -100], dvposbins);
      temp = layer1xnbp + layer2xnbp + layer3xnbp + layer4xnbp + layer5xnbp + layer6xnbp;
      bar(dvposbins, [layer6xnbp(:) layer5xnbp(:) layer4xnbp(:) layer3xnbp(:) layer2xnbp(:) layer1xnbp(:)], 'stacked');
      axis([-1 6 0 max([round(1.1*max(temp)) max(temp)+1])]);
      set(gca,'fontsize', 8);
      children = get(gca,'children');

      set(children(1),'facecolor', [1 1 1]); % layer 1
      set(children(2),'facecolor', [.85 .85 .85]); % layer 2
      set(children(3),'facecolor', [.7 .7 .7]); % layer 3
      set(children(4),'facecolor', [.55 .55 .55]); % layer 4
      set(children(5),'facecolor', [.4 .4 .4]); % layer 5
      set(children(6),'facecolor', [0 0 0]); % layer 6
      title('3dB LowPass');
      set(gca,'xticklabel', '', 'tickdir', 'out');

   elseif ( i == 6 )

      layer1xnbp = histc([layer1xlp6db -100], dvposbins);
      layer2xnbp = histc([layer2xlp6db -100], dvposbins);
      layer3xnbp = histc([layer3xlp6db -100], dvposbins);
      layer4xnbp = histc([layer4xlp6db -100], dvposbins);
      layer5xnbp = histc([layer5xlp6db -100], dvposbins);
      layer6xnbp = histc([layer6xlp6db -100], dvposbins);
      temp = layer1xnbp + layer2xnbp + layer3xnbp + layer4xnbp + layer5xnbp + layer6xnbp;
      bar(dvposbins, [layer6xnbp(:) layer5xnbp(:) layer4xnbp(:) layer3xnbp(:) layer2xnbp(:) layer1xnbp(:)], 'stacked');
      axis([-1 6 0 max([round(1.1*max(temp)) max(temp)+1])]);
      set(gca,'fontsize', 8);
      children = get(gca,'children');

      set(children(1),'facecolor', [1 1 1]); % layer 1
      set(children(2),'facecolor', [.85 .85 .85]); % layer 2
      set(children(3),'facecolor', [.7 .7 .7]); % layer 3
      set(children(4),'facecolor', [.55 .55 .55]); % layer 4
      set(children(5),'facecolor', [.4 .4 .4]); % layer 5
      set(children(6),'facecolor', [0 0 0]); % layer 6
      title('6dB LowPass');
      set(gca, 'tickdir', 'out');

   end

end % (for i)


if ( nargin == 2 )
   if ( isstr(exp) )
      suptitle([exp '  Spectral Modulation Parameter Topography']);
   end
end

orient landscape;


% 
% % *********************************************************
% %           Plot the best TMF Layer Comparisons
% % *********************************************************
% 
% figure;
% 
% % ---------------------------------------------------
% %            Plot the raw layer data
% % ---------------------------------------------------
% 
% ylabeltext{1} = 'Mean filter';
% ylabeltext{2} = 'Mean filter smoothed';
% ylabeltext{3} = 'Max filter';
% ylabeltext{4} = 'Max filter smoothed';
% 
% colortype{1} = 'b-';
% colortype{2} = 'g-';
% colortype{3} = 'r-';
% colortype{4} = 'm:';
% colortype{5} = 'c-';
% colortype{6} = 'k-';
% 
% for i = 1:4  % 4 different smoothing procedures used
% 
%    subplot(4,3,3*i-2);
%    hold on;
% 
%    for j = 1:6  % layers 1 thru 6
% 
%       if ( ~isempty( eval(['btmflayerdata{' num2str(j) '}{1}']) ) )
%          plot(eval(['btmflayerdata{' num2str(j) '}{1}']), ...
%               eval(['btmflayerdata{' num2str(j) '}{' num2str(i+1) '}']), ...
%               colortype{j}, 'linewidth', 1);
%       else
%          plot(-2, 5, colortype{j}, 'linewidth', 1);
%       end
% 
%    end % (for j)
% 
%    ylabel(ylabeltext{i});
%    set(gca,'xtick', -1:6, 'xticklabel', -1:6);
%    set(gca,'ytick', 0:10:40, 'yticklabel', 0:10:40);
%    set(gca,'fontsize', 8);
%    axis([-1 6 -5 45]);
% 
%    if ( i == 4 )
%       xlabel('Ventral -- > Dorsal (mm)');
%       legend('L1', 'L2', 'L3', 'L4', 'L5', 'L6', 0);
%    elseif ( i == 1 )
%       title('Smoothing Procedures');
%       set(gca, 'xticklabel', '');
%    else
%       set(gca, 'xticklabel', '');
%    end
% 
% end % (for i)
% 
% 
% 
% % ---------------------------------------------------
% %     Plot the layer comparisons - normalization
% % ---------------------------------------------------
% 
% colortype{1} = 'b-';
% colortype{2} = 'g-';
% colortype{3} = 'r-';
% colortype{4} = 'c-';
% colortype{5} = 'k-';
% 
% layers = [1 2 3 5 6];
% 
% for i = 1:4  % 4 different smoothing procedures used
% 
%    subplot(4,3,3*i-1);
%    hold on;
% 
%    for j = 1:length(layers)  % layers 1 thru 6
% 
%       layerpos = eval(['btmflayerdata{' num2str(layers(j)) '}{1}']);
%       layer4pos = eval(['btmflayerdata{4}{1}']);
% 
%       [pos, ilayer, i4] = intersect(layerpos, layer4pos);
% 
%       if ( ~isempty( pos ) )
%          rlayer = eval(['btmflayerdata{' num2str(layers(j)) '}{' num2str(i+1) '}(ilayer)']);
%          r4 = eval(['btmflayerdata{' num2str(4) '}{' num2str(i+1) '}(i4)']);
%          plot(pos, rlayer ./ r4, colortype{j}, 'linewidth', 1);
%       else
%          plot(-10, -10, colortype{j}, 'linewidth', 1);
%       end
% 
%    end % (for j)
% 
%    set(gca,'fontsize', 8);
%    set(gca,'xtick', -1:6, 'xticklabel', -1:6);
%    axis([-1 6 0.5 1.75]);
%    
%    if ( i == 4 )
%       xlabel('Ventral -- > Dorsal (mm)');
%    elseif ( i == 3 )
%       lh = legend('1/4', '2/4', '3/4', '5/4', '6/4', 1);
%       set(gca, 'xticklabel', '');
%    elseif ( i == 1 )
%       title('Ratios relative to Layer 4');
%       set(gca, 'xticklabel', '');
%    else
%       set(gca, 'xticklabel', '');
%    end
% 
% end % (for i)
% 
% % ---------------------------------------------------
% %     Plot the layer comparisons - differences
% % ---------------------------------------------------
% 
% colortype{1} = 'b-';
% colortype{2} = 'g-';
% colortype{3} = 'r-';
% colortype{4} = 'c-';
% colortype{5} = 'k-';
% 
% layers = [1 2 3 5 6];
% 
% for i = 1:4  % 4 different smoothing procedures used
% 
%    subplot(4,3,3*i);
%    hold on;
% 
%    for j = 1:length(layers)  % layers 1 thru 6
% 
%       layerpos = eval(['btmflayerdata{' num2str(layers(j)) '}{1}']);
%       layer4pos = eval(['btmflayerdata{4}{1}']);
% 
%       [pos, ilayer, i4] = intersect(layerpos, layer4pos);
% 
%       if ( ~isempty( pos ) )
%          rlayer = eval(['btmflayerdata{' num2str(layers(j)) '}{' num2str(i+1) '}(ilayer)']);
%          r4 = eval(['btmflayerdata{' num2str(4) '}{' num2str(i+1) '}(i4)']);
%          plot(pos, abs(rlayer-r4) ./ r4 * 100, colortype{j}, 'linewidth', 1);
%       else
%          plot(-10, -10, colortype{j}, 'linewidth', 1);
%       end
% 
%    end % (for j)
% 
%    set(gca,'fontsize', 8);
%    set(gca,'xtick', -1:6, 'xticklabel', -1:6);
%    axis([-1 6 0 50]);
%    
%    if ( i == 4 )
%       xlabel('Ventral -- > Dorsal (mm)');
%    elseif ( i == 3 )
%       lh = legend('(1-4)/4', '(2-4)/4', '(3-4)/4', '(5-4)/4', '(6-4)/4', 1);
%       set(gca, 'xticklabel', '');
%    elseif ( i == 1 )
%       title('Percent difference relative to Layer 4');
%       set(gca, 'xticklabel', '');
%    else
%       set(gca, 'xticklabel', '');
%    end
% 
% end % (for i)
% 
% orient tall;
% 
% if ( nargin == 2 )
%    if ( isstr(exp) )
%       suptitle([exp '  Best TMF Layer Comparisons']);
%    end
% end
% 
% 
% % *********************************************************
% %       Plot the tMTF widths at 3dB Layer Comparisons
% % *********************************************************
% 
% figure;
% 
% % ---------------------------------------------------
% %            Plot the raw layer data
% % ---------------------------------------------------
% 
% ylabeltext{1} = 'Mean filter';
% ylabeltext{2} = 'Mean filter smoothed';
% ylabeltext{3} = 'Max filter';
% ylabeltext{4} = 'Max filter smoothed';
% 
% colortype{1} = 'b-';
% colortype{2} = 'g-';
% colortype{3} = 'r-';
% colortype{4} = 'm:';
% colortype{5} = 'c-';
% colortype{6} = 'k-';
% 
% for i = 1:4  % 4 different smoothing procedures used
% 
%    subplot(4,3,3*i-2);
%    hold on;
% 
%    for j = 1:6  % layers 1 thru 6
% 
%       if ( ~isempty( eval(['tw3dblayerdata{' num2str(j) '}{1}']) ) )
%          plot(eval(['tw3dblayerdata{' num2str(j) '}{1}']), ...
%               eval(['tw3dblayerdata{' num2str(j) '}{' num2str(i+1) '}']), ...
%               colortype{j}, 'linewidth', 1);
%       else
%          plot(-10, -10, colortype{j}, 'linewidth', 1);
%       end
% 
%    end % (for j)
% 
%    ylabel(ylabeltext{i});
%    set(gca,'xtick', -1:6, 'xticklabel', -1:6);
%    set(gca,'ytick', 0:10:40, 'yticklabel', 0:10:40);
%    set(gca,'fontsize', 8);
%    axis([-1 6 -5 45]);
% 
%    if ( i == 4 )
%       xlabel('Ventral -- > Dorsal (mm)');
%       legend('L1', 'L2', 'L3', 'L4', 'L5', 'L6', 0);
%    elseif ( i == 1 )
%       title('Smoothing Procedures');
%       set(gca, 'xticklabel', '');
%    else
%       set(gca, 'xticklabel', '');
%    end
% 
% end % (for i)
% 
% 
% 
% % ---------------------------------------------------
% %     Plot the layer comparisons - normalization
% % ---------------------------------------------------
% 
% colortype{1} = 'b-';
% colortype{2} = 'g-';
% colortype{3} = 'r-';
% colortype{4} = 'c-';
% colortype{5} = 'k-';
% 
% layers = [1 2 3 5 6];
% 
% for i = 1:4  % 4 different smoothing procedures used
% 
%    subplot(4,3,3*i-1);
%    hold on;
% 
%    for j = 1:length(layers)  % layers 1 thru 6
% 
%       layerpos = eval(['tw3dblayerdata{' num2str(layers(j)) '}{1}']);
%       layer4pos = eval(['tw3dblayerdata{4}{1}']);
% 
%       [pos, ilayer, i4] = intersect(layerpos, layer4pos);
% 
%       if ( ~isempty( pos ) )
%          rlayer = eval(['tw3dblayerdata{' num2str(layers(j)) '}{' num2str(i+1) '}(ilayer)']);
%          r4 = eval(['tw3dblayerdata{' num2str(4) '}{' num2str(i+1) '}(i4)']);
%          plot(pos, rlayer ./ r4, colortype{j}, 'linewidth', 1);
%       else
%          plot(-10, -10, colortype{j}, 'linewidth', 1);
%       end
% 
%    end % (for j)
% 
%    set(gca,'fontsize', 8);
%    set(gca,'xtick', -1:6, 'xticklabel', -1:6);
%    axis([-1 6 0.5 1.75]);
%    
%    if ( i == 4 )
%       xlabel('Ventral -- > Dorsal (mm)');
%    elseif ( i == 3 )
%       lh = legend('1/4', '2/4', '3/4', '5/4', '6/4', 1);
%       set(gca, 'xticklabel', '');
%    elseif ( i == 1 )
%       title('Ratios relative to Layer 4');
%       set(gca, 'xticklabel', '');
%    else
%       set(gca, 'xticklabel', '');
%    end
% 
% end % (for i)
% 
% % ---------------------------------------------------
% %     Plot the layer comparisons - differences
% % ---------------------------------------------------
% 
% colortype{1} = 'b-';
% colortype{2} = 'g-';
% colortype{3} = 'r-';
% colortype{4} = 'c-';
% colortype{5} = 'k-';
% 
% layers = [1 2 3 5 6];
% 
% for i = 1:4  % 4 different smoothing procedures used
% 
%    subplot(4,3,3*i);
%    hold on;
% 
%    for j = 1:length(layers)  % layers 1 thru 6
% 
%       layerpos = eval(['tw3dblayerdata{' num2str(layers(j)) '}{1}']);
%       layer4pos = eval(['tw3dblayerdata{4}{1}']);
% 
%       [pos, ilayer, i4] = intersect(layerpos, layer4pos);
% 
%       if ( ~isempty( pos ) )
%          rlayer = eval(['tw3dblayerdata{' num2str(layers(j)) '}{' num2str(i+1) '}(ilayer)']);
%          r4 = eval(['tw3dblayerdata{' num2str(4) '}{' num2str(i+1) '}(i4)']);
%          plot(pos, abs(rlayer-r4) ./ r4 * 100, colortype{j}, 'linewidth', 1);
%       else
%          plot(-10, -10, colortype{j}, 'linewidth', 1);
%       end
% 
%    end % (for j)
% 
%    set(gca,'fontsize', 8);
%    set(gca,'xtick', -1:6, 'xticklabel', -1:6);
%    axis([-1 6 0 50]);
%    
%    if ( i == 4 )
%       xlabel('Ventral -- > Dorsal (mm)');
%    elseif ( i == 3 )
%       lh = legend('(1-4)/4', '(2-4)/4', '(3-4)/4', '(5-4)/4', '(6-4)/4', 1);
%       set(gca, 'xticklabel', '');
%    elseif ( i == 1 )
%       title('Percent difference relative to Layer 4');
%       set(gca, 'xticklabel', '');
%    else
%       set(gca, 'xticklabel', '');
%    end
% 
% end % (for i)
% 
% orient tall;
% 
% if ( nargin == 2 )
%    if ( isstr(exp) )
%       suptitle([exp '  tMTF 3dB Width Layer Comparisons']);
%    end
% end
% 
% 
% 
% 
% 
% % *********************************************************
% %       Plot the tMTF widths at 6dB Layer Comparisons
% % *********************************************************
% 
% figure;
% 
% % ---------------------------------------------------
% %            Plot the raw layer data
% % ---------------------------------------------------
% 
% ylabeltext{1} = 'Mean filter';
% ylabeltext{2} = 'Mean filter smoothed';
% ylabeltext{3} = 'Max filter';
% ylabeltext{4} = 'Max filter smoothed';
% 
% colortype{1} = 'b-';
% colortype{2} = 'g-';
% colortype{3} = 'r-';
% colortype{4} = 'm:';
% colortype{5} = 'c-';
% colortype{6} = 'k-';
% 
% for i = 1:4  % 4 different smoothing procedures used
% 
%    subplot(4,3,3*i-2);
%    hold on;
% 
%    for j = 1:6  % layers 1 thru 6
% 
%       if ( ~isempty( eval(['tw6dblayerdata{' num2str(j) '}{1}']) ) )
%          plot(eval(['tw6dblayerdata{' num2str(j) '}{1}']), ...
%               eval(['tw6dblayerdata{' num2str(j) '}{' num2str(i+1) '}']), ...
%               colortype{j}, 'linewidth', 1);
%       else
%          plot(-10, -10, colortype{j}, 'linewidth', 1);
%       end
% 
%    end % (for j)
% 
%    ylabel(ylabeltext{i});
%    set(gca,'xtick', -1:6, 'xticklabel', -1:6);
%    set(gca,'ytick', 0:10:40, 'yticklabel', 0:10:40);
%    set(gca,'fontsize', 8);
%    axis([-1 6 -5 45]);
% 
%    if ( i == 4 )
%       xlabel('Ventral -- > Dorsal (mm)');
%       legend('L1', 'L2', 'L3', 'L4', 'L5', 'L6', 0);
%    elseif ( i == 1 )
%       title('Smoothing Procedures');
%       set(gca, 'xticklabel', '');
%    else
%       set(gca, 'xticklabel', '');
%    end
% 
% end % (for i)
% 
% 
% % ---------------------------------------------------
% %     Plot the layer comparisons - normalization
% % ---------------------------------------------------
% 
% colortype{1} = 'b-';
% colortype{2} = 'g-';
% colortype{3} = 'r-';
% colortype{4} = 'c-';
% colortype{5} = 'k-';
% 
% layers = [1 2 3 5 6];
% 
% for i = 1:4  % 4 different smoothing procedures used
% 
%    subplot(4,3,3*i-1);
%    hold on;
% 
%    for j = 1:length(layers)  % layers 1 thru 6
% 
%       layerpos = eval(['tw6dblayerdata{' num2str(layers(j)) '}{1}']);
%       layer4pos = eval(['tw6dblayerdata{4}{1}']);
% 
%       [pos, ilayer, i4] = intersect(layerpos, layer4pos);
% 
%       if ( ~isempty( pos ) )
%          rlayer = eval(['tw6dblayerdata{' num2str(layers(j)) '}{' num2str(i+1) '}(ilayer)']);
%          r4 = eval(['tw6dblayerdata{' num2str(4) '}{' num2str(i+1) '}(i4)']);
%          plot(pos, rlayer ./ r4, colortype{j}, 'linewidth', 1);
%       else
%          plot(-10, -10, colortype{j}, 'linewidth', 1);
%       end
% 
%    end % (for j)
% 
%    set(gca,'fontsize', 8);
%    set(gca,'xtick', -1:6, 'xticklabel', -1:6);
%    axis([-1 6 0.5 1.75]);
%    
%    if ( i == 4 )
%       xlabel('Ventral -- > Dorsal (mm)');
%    elseif ( i == 3 )
%       lh = legend('1/4', '2/4', '3/4', '5/4', '6/4', 1);
%       set(gca, 'xticklabel', '');
%    elseif ( i == 1 )
%       title('Ratios relative to Layer 4');
%       set(gca, 'xticklabel', '');
%    else
%       set(gca, 'xticklabel', '');
%    end
% 
% end % (for i)
% 
% 
% % ---------------------------------------------------
% %     Plot the layer comparisons - differences
% % ---------------------------------------------------
% 
% colortype{1} = 'b-';
% colortype{2} = 'g-';
% colortype{3} = 'r-';
% colortype{4} = 'c-';
% colortype{5} = 'k-';
% 
% layers = [1 2 3 5 6];
% 
% for i = 1:4  % 4 different smoothing procedures used
% 
%    subplot(4,3,3*i);
%    hold on;
% 
%    for j = 1:length(layers)  % layers 1 thru 6
% 
%       layerpos = eval(['tw6dblayerdata{' num2str(layers(j)) '}{1}']);
%       layer4pos = eval(['tw6dblayerdata{4}{1}']);
% 
%       [pos, ilayer, i4] = intersect(layerpos, layer4pos);
% 
%       if ( ~isempty( pos ) )
%          rlayer = eval(['tw6dblayerdata{' num2str(layers(j)) '}{' num2str(i+1) '}(ilayer)']);
%          r4 = eval(['tw6dblayerdata{' num2str(4) '}{' num2str(i+1) '}(i4)']);
%          plot(pos, abs(rlayer-r4) ./ r4 * 100, colortype{j}, 'linewidth', 1);
%       else
%          plot(-10, -10, colortype{j}, 'linewidth', 1);
%       end
% 
%    end % (for j)
% 
%    set(gca,'fontsize', 8);
%    set(gca,'xtick', -1:6, 'xticklabel', -1:6);
%    axis([-1 6 0 50]);
%    
%    if ( i == 4 )
%       xlabel('Ventral -- > Dorsal (mm)');
%    elseif ( i == 3 )
%       lh = legend('(1-4)/4', '(2-4)/4', '(3-4)/4', '(5-4)/4', '(6-4)/4', 1);
%       set(gca, 'xticklabel', '');
%    elseif ( i == 1 )
%       title('Percent difference relative to Layer 4');
%       set(gca, 'xticklabel', '');
%    else
%       set(gca, 'xticklabel', '');
%    end
% 
% end % (for i)
% 
% orient tall;
% 
% if ( nargin == 2 )
%    if ( isstr(exp) )
%       suptitle([exp '  tMTF 6dB Width Layer Comparisons']);
%    end
% end
% 
% 
% 
% 
% 
% % *********************************************************
% %           Plot the best SMF Layer Comparisons
% % *********************************************************
% 
% figure;
% 
% % ---------------------------------------------------
% %            Plot the raw layer data
% % ---------------------------------------------------
% 
% ylabeltext{1} = 'Mean filter';
% ylabeltext{2} = 'Mean filter smoothed';
% ylabeltext{3} = 'Max filter';
% ylabeltext{4} = 'Max filter smoothed';
% 
% colortype{1} = 'b-';
% colortype{2} = 'g-';
% colortype{3} = 'r-';
% colortype{4} = 'm:';
% colortype{5} = 'c-';
% colortype{6} = 'k-';
% 
% for i = 1:4  % 4 different smoothing procedures used
% 
%    subplot(4,3,3*i-2);
%    hold on;
% 
%    for j = 1:6  % layers 1 thru 6
% 
%       if ( ~isempty( eval(['bsmflayerdata{' num2str(j) '}{1}']) ) )
%          plot(eval(['bsmflayerdata{' num2str(j) '}{1}']), ...
%               eval(['bsmflayerdata{' num2str(j) '}{' num2str(i+1) '}']), ...
%               colortype{j}, 'linewidth', 1);
%       else
%          plot(-10, -10, colortype{j}, 'linewidth', 1);
%       end
% 
%    end % (for j)
% 
%    ylabel(ylabeltext{i});
%    set(gca,'xtick', -1:6, 'xticklabel', -1:6);
%    set(gca,'ytick', 0:1:4, 'yticklabel', 0:1:4, 'tickdir', 'out');
%    set(gca,'fontsize', 8);
%    axis([-1 6 -1 5]);
% 
%    if ( i == 4 )
%       xlabel('Ventral -- > Dorsal (mm)');
%       legend('L1', 'L2', 'L3', 'L4', 'L5', 'L6', 0);
%    elseif ( i == 1 )
%       title('Smoothing Procedures');
%       set(gca, 'xticklabel', '');
%    else
%       set(gca, 'xticklabel', '');
%    end
% 
% end % (for i)
% 
% 
% % ---------------------------------------------------
% %     Plot the layer comparisons - normalization
% % ---------------------------------------------------
% 
% colortype{1} = 'b-';
% colortype{2} = 'g-';
% colortype{3} = 'r-';
% colortype{4} = 'c-';
% colortype{5} = 'k-';
% 
% layers = [1 2 3 5 6];
% 
% for i = 1:4  % 4 different smoothing procedures used
% 
%    subplot(4,3,3*i-1);
%    hold on;
% 
%    for j = 1:length(layers)  % layers 1 thru 6
% 
%       layerpos = eval(['bsmflayerdata{' num2str(layers(j)) '}{1}']);
%       layer4pos = eval(['bsmflayerdata{4}{1}']);
% 
%       [pos, ilayer, i4] = intersect(layerpos, layer4pos);
% 
%       if ( ~isempty( pos ) )
%          rlayer = eval(['bsmflayerdata{' num2str(layers(j)) '}{' num2str(i+1) '}(ilayer)']);
%          r4 = eval(['bsmflayerdata{' num2str(4) '}{' num2str(i+1) '}(i4)']);
%          plot(pos, rlayer ./ r4, colortype{j}, 'linewidth', 1);
%       else
%          plot(-10, -10, colortype{j}, 'linewidth', 1);
%       end
% 
%    end % (for j)
% 
%    set(gca,'fontsize', 8);
%    set(gca,'xtick', -1:6, 'xticklabel', -1:6);
%    axis([-1 6 0.5 1.75]);
%    
%    if ( i == 4 )
%       xlabel('Ventral -- > Dorsal (mm)');
%    elseif ( i == 3 )
%       lh = legend('1/4', '2/4', '3/4', '5/4', '6/4', 1);
%       set(gca, 'xticklabel', '');
%    elseif ( i == 1 )
%       title('Ratios relative to Layer 4');
%       set(gca, 'xticklabel', '');
%    else
%       set(gca, 'xticklabel', '');
%    end
% 
% end % (for i)
% 
% 
% % ---------------------------------------------------
% %     Plot the layer comparisons - differences
% % ---------------------------------------------------
% 
% colortype{1} = 'b-';
% colortype{2} = 'g-';
% colortype{3} = 'r-';
% colortype{4} = 'c-';
% colortype{5} = 'k-';
% 
% layers = [1 2 3 5 6];
% 
% for i = 1:4  % 4 different smoothing procedures used
% 
%    subplot(4,3,3*i);
%    hold on;
% 
%    for j = 1:length(layers)  % layers 1 thru 6
% 
%       layerpos = eval(['bsmflayerdata{' num2str(layers(j)) '}{1}']);
%       layer4pos = eval(['bsmflayerdata{4}{1}']);
% 
%       [pos, ilayer, i4] = intersect(layerpos, layer4pos);
% 
%       if ( ~isempty( pos ) )
%          rlayer = eval(['bsmflayerdata{' num2str(layers(j)) '}{' num2str(i+1) '}(ilayer)']);
%          r4 = eval(['bsmflayerdata{' num2str(4) '}{' num2str(i+1) '}(i4)']);
%          plot(pos, abs(rlayer-r4) ./ r4 * 100, colortype{j}, 'linewidth', 1);
%       else
%          plot(-10, -10, colortype{j}, 'linewidth', 1);
%       end
% 
%    end % (for j)
% 
%    set(gca,'fontsize', 8);
%    set(gca,'xtick', -1:6, 'xticklabel', -1:6);
%    axis([-1 6 0 50]);
%    
%    if ( i == 4 )
%       xlabel('Ventral -- > Dorsal (mm)');
%    elseif ( i == 3 )
%       lh = legend('(1-4)/4', '(2-4)/4', '(3-4)/4', '(5-4)/4', '(6-4)/4', 1);
%       set(gca, 'xticklabel', '');
%    elseif ( i == 1 )
%       title('Percent difference relative to Layer 4');
%       set(gca, 'xticklabel', '');
%    else
%       set(gca, 'xticklabel', '');
%    end
% 
% end % (for i)
% 
% orient tall;
% 
% if ( nargin == 2 )
%    if ( isstr(exp) )
%       suptitle([exp '  Best SMF Layer Comparisons']);
%    end
% end
% 
% 
% % *********************************************************
% %       Plot the sMTF widths at 3dB Layer Comparisons
% % *********************************************************
% 
% figure;
% 
% % ---------------------------------------------------
% %            Plot the raw layer data
% % ---------------------------------------------------
% 
% ylabeltext{1} = 'Mean filter';
% ylabeltext{2} = 'Mean filter smoothed';
% ylabeltext{3} = 'Max filter';
% ylabeltext{4} = 'Max filter smoothed';
% 
% colortype{1} = 'b-';
% colortype{2} = 'g-';
% colortype{3} = 'r-';
% colortype{4} = 'm:';
% colortype{5} = 'c-';
% colortype{6} = 'k-';
% 
% for i = 1:4  % 4 different smoothing procedures used
% 
%    subplot(4,3,3*i-2);
%    hold on;
% 
%    for j = 1:6  % layers 1 thru 6
% 
%       if ( ~isempty( eval(['sw3dblayerdata{' num2str(j) '}{1}']) ) )
%          plot(eval(['sw3dblayerdata{' num2str(j) '}{1}']), ...
%               eval(['sw3dblayerdata{' num2str(j) '}{' num2str(i+1) '}']), ...
%               colortype{j}, 'linewidth', 1);
%       else
%          plot(-10, -10, colortype{j}, 'linewidth', 1);
%       end
% 
%    end % (for j)
% 
%    ylabel(ylabeltext{i});
%    set(gca,'xtick', -1:6, 'xticklabel', -1:6);
%    set(gca,'ytick', 0:1:4, 'yticklabel', 0:1:4, 'tickdir', 'out');
%    set(gca,'fontsize', 8);
%    axis([-1 6 -1 5]);
% 
%    if ( i == 4 )
%       xlabel('Ventral -- > Dorsal (mm)');
%       legend('L1', 'L2', 'L3', 'L4', 'L5', 'L6', 0);
%    elseif ( i == 1 )
%       title('Smoothing Procedures');
%       set(gca, 'xticklabel', '');
%    else
%       set(gca, 'xticklabel', '');
%    end
% 
% end % (for i)
% 
% 
% 
% % ---------------------------------------------------
% %     Plot the layer comparisons - normalization
% % ---------------------------------------------------
% 
% colortype{1} = 'b-';
% colortype{2} = 'g-';
% colortype{3} = 'r-';
% colortype{4} = 'c-';
% colortype{5} = 'k-';
% 
% layers = [1 2 3 5 6];
% 
% for i = 1:4  % 4 different smoothing procedures used
% 
%    subplot(4,3,3*i-1);
%    hold on;
% 
%    for j = 1:length(layers)  % layers 1 thru 6
% 
%       layerpos = eval(['sw3dblayerdata{' num2str(layers(j)) '}{1}']);
%       layer4pos = eval(['sw3dblayerdata{4}{1}']);
% 
%       [pos, ilayer, i4] = intersect(layerpos, layer4pos);
% 
%       if ( ~isempty( pos ) )
%          rlayer = eval(['sw3dblayerdata{' num2str(layers(j)) '}{' num2str(i+1) '}(ilayer)']);
%          r4 = eval(['sw3dblayerdata{' num2str(4) '}{' num2str(i+1) '}(i4)']);
%          plot(pos, rlayer ./ r4, colortype{j}, 'linewidth', 1);
%       else
%          plot(-10, -10, colortype{j}, 'linewidth', 1);
%       end
% 
%    end % (for j)
% 
%    set(gca,'fontsize', 8);
%    set(gca,'xtick', -1:6, 'xticklabel', -1:6);
%    axis([-1 6 0.5 1.75]);
%    
%    if ( i == 4 )
%       xlabel('Ventral -- > Dorsal (mm)');
%    elseif ( i == 3 )
%       lh = legend('1/4', '2/4', '3/4', '5/4', '6/4', 1);
%       set(gca, 'xticklabel', '');
%    elseif ( i == 1 )
%       title('Ratios relative to Layer 4');
%       set(gca, 'xticklabel', '');
%    else
%       set(gca, 'xticklabel', '');
%    end
% 
% end % (for i)
% 
% % ---------------------------------------------------
% %     Plot the layer comparisons - differences
% % ---------------------------------------------------
% 
% colortype{1} = 'b-';
% colortype{2} = 'g-';
% colortype{3} = 'r-';
% colortype{4} = 'c-';
% colortype{5} = 'k-';
% 
% layers = [1 2 3 5 6];
% 
% for i = 1:4  % 4 different smoothing procedures used
% 
%    subplot(4,3,3*i);
%    hold on;
% 
%    for j = 1:length(layers)  % layers 1 thru 6
% 
%       layerpos = eval(['sw3dblayerdata{' num2str(layers(j)) '}{1}']);
%       layer4pos = eval(['sw3dblayerdata{4}{1}']);
% 
%       [pos, ilayer, i4] = intersect(layerpos, layer4pos);
% 
%       if ( ~isempty( pos ) )
%          rlayer = eval(['sw3dblayerdata{' num2str(layers(j)) '}{' num2str(i+1) '}(ilayer)']);
%          r4 = eval(['sw3dblayerdata{' num2str(4) '}{' num2str(i+1) '}(i4)']);
%          plot(pos, abs(rlayer-r4) ./ r4 * 100, colortype{j}, 'linewidth', 1);
%       else
%          plot(-10, -10, colortype{j}, 'linewidth', 1);
%       end
% 
%    end % (for j)
% 
%    set(gca,'fontsize', 8);
%    set(gca,'xtick', -1:6, 'xticklabel', -1:6);
%    axis([-1 6 0 50]);
%    
%    if ( i == 4 )
%       xlabel('Ventral -- > Dorsal (mm)');
%    elseif ( i == 3 )
%       lh = legend('(1-4)/4', '(2-4)/4', '(3-4)/4', '(5-4)/4', '(6-4)/4', 1);
%       set(gca, 'xticklabel', '');
%    elseif ( i == 1 )
%       title('Percent difference relative to Layer 4');
%       set(gca, 'xticklabel', '');
%    else
%       set(gca, 'xticklabel', '');
%    end
% 
% end % (for i)
% 
% orient tall;
% 
% if ( nargin == 2 )
%    if ( isstr(exp) )
%       suptitle([exp '  sMTF 3dB Width Layer Comparisons']);
%    end
% end
% 
% 
% 
% 
% 
% % *********************************************************
% %       Plot the sMTF widths at 6dB Layer Comparisons
% % *********************************************************
% 
% figure;
% 
% % ---------------------------------------------------
% %            Plot the raw layer data
% % ---------------------------------------------------
% 
% ylabeltext{1} = 'Mean filter';
% ylabeltext{2} = 'Mean filter smoothed';
% ylabeltext{3} = 'Max filter';
% ylabeltext{4} = 'Max filter smoothed';
% 
% colortype{1} = 'b-';
% colortype{2} = 'g-';
% colortype{3} = 'r-';
% colortype{4} = 'm:';
% colortype{5} = 'c-';
% colortype{6} = 'k-';
% 
% for i = 1:4  % 4 different smoothing procedures used
% 
%    subplot(4,3,3*i-2);
%    hold on;
% 
%    for j = 1:6  % layers 1 thru 6
% 
%       if ( ~isempty( eval(['sw6dblayerdata{' num2str(j) '}{1}']) ) )
%          plot(eval(['sw6dblayerdata{' num2str(j) '}{1}']), ...
%               eval(['sw6dblayerdata{' num2str(j) '}{' num2str(i+1) '}']), ...
%               colortype{j}, 'linewidth', 1);
%       else
%          plot(-10, -10, colortype{j}, 'linewidth', 1);
%       end
% 
%    end % (for j)
% 
%    ylabel(ylabeltext{i});
%    set(gca,'xtick', -1:6, 'xticklabel', -1:6);
%    set(gca,'ytick', 0:1:4, 'yticklabel', 0:1:4, 'tickdir', 'out');
%    set(gca,'fontsize', 8);
%    axis([-1 6 -1 5]);
% 
%    if ( i == 4 )
%       xlabel('Ventral -- > Dorsal (mm)');
%       legend('L1', 'L2', 'L3', 'L4', 'L5', 'L6', 0);
%    elseif ( i == 1 )
%       title('Smoothing Procedures');
%       set(gca, 'xticklabel', '');
%    else
%       set(gca, 'xticklabel', '');
%    end
% 
% end % (for i)
% 
% 
% 
% % ---------------------------------------------------
% %     Plot the layer comparisons - normalization
% % ---------------------------------------------------
% 
% colortype{1} = 'b-';
% colortype{2} = 'g-';
% colortype{3} = 'r-';
% colortype{4} = 'c-';
% colortype{5} = 'k-';
% 
% layers = [1 2 3 5 6];
% 
% for i = 1:4  % 4 different smoothing procedures used
% 
%    subplot(4,3,3*i-1);
%    hold on;
% 
%    for j = 1:length(layers)  % layers 1 thru 6
% 
%       layerpos = eval(['sw6dblayerdata{' num2str(layers(j)) '}{1}']);
%       layer4pos = eval(['sw6dblayerdata{4}{1}']);
% 
%       [pos, ilayer, i4] = intersect(layerpos, layer4pos);
% 
%       if ( ~isempty( pos ) )
%          rlayer = eval(['sw6dblayerdata{' num2str(layers(j)) '}{' num2str(i+1) '}(ilayer)']);
%          r4 = eval(['sw6dblayerdata{' num2str(4) '}{' num2str(i+1) '}(i4)']);
%          plot(pos, rlayer ./ r4, colortype{j}, 'linewidth', 1);
%       else
%          plot(-10, -10, colortype{j}, 'linewidth', 1);
%       end
% 
%    end % (for j)
% 
%    set(gca,'fontsize', 8);
%    set(gca,'xtick', -1:6, 'xticklabel', -1:6);
%    axis([-1 6 0.5 1.75]);
%    
%    if ( i == 4 )
%       xlabel('Ventral -- > Dorsal (mm)');
%    elseif ( i == 3 )
%       lh = legend('1/4', '2/4', '3/4', '5/4', '6/4', 1);
%       set(gca, 'xticklabel', '');
%    elseif ( i == 1 )
%       title('Ratios relative to Layer 4');
%       set(gca, 'xticklabel', '');
%    else
%       set(gca, 'xticklabel', '');
%    end
% 
% end % (for i)
% 
% % ---------------------------------------------------
% %     Plot the layer comparisons - differences
% % ---------------------------------------------------
% 
% colortype{1} = 'b-';
% colortype{2} = 'g-';
% colortype{3} = 'r-';
% colortype{4} = 'c-';
% colortype{5} = 'k-';
% 
% layers = [1 2 3 5 6];
% 
% for i = 1:4  % 4 different smoothing procedures used
% 
%    subplot(4,3,3*i);
%    hold on;
% 
%    for j = 1:length(layers)  % layers 1 thru 6
% 
%       layerpos = eval(['sw6dblayerdata{' num2str(layers(j)) '}{1}']);
%       layer4pos = eval(['sw6dblayerdata{4}{1}']);
% 
%       [pos, ilayer, i4] = intersect(layerpos, layer4pos);
% 
%       if ( ~isempty( pos ) )
%          rlayer = eval(['sw6dblayerdata{' num2str(layers(j)) '}{' num2str(i+1) '}(ilayer)']);
%          r4 = eval(['sw6dblayerdata{' num2str(4) '}{' num2str(i+1) '}(i4)']);
%          plot(pos, abs(rlayer-r4) ./ r4 * 100, colortype{j}, 'linewidth', 1);
%       else
%          plot(-10, -10, colortype{j}, 'linewidth', 1);
%       end
% 
%    end % (for j)
% 
%    set(gca,'fontsize', 8);
%    set(gca,'xtick', -1:6, 'xticklabel', -1:6);
%    axis([-1 6 0 50]);
% 
%    if ( i == 4 )
%       xlabel('Ventral -- > Dorsal (mm)');
%    elseif ( i == 3 )
%       lh = legend('(1-4)/4', '(2-4)/4', '(3-4)/4', '(5-4)/4', '(6-4)/4', 1);
%       set(gca, 'xticklabel', '');
%    elseif ( i == 1 )
%       title('Percent difference relative to Layer 4');
%       set(gca, 'xticklabel', '');
%    else
%       set(gca, 'xticklabel', '');
%    end
% 
% end % (for i)
% 
% orient tall;
% 
% if ( nargin == 2 )
%    if ( isstr(exp) )
%       suptitle([exp '  sMTF 6dB Width Layer Comparisons']);
%    end
% end



btmflayerdata{end+1} = 'TBMF';
tw3dblayerdata{end+1} = 'TW3dB';
tw6dblayerdata{end+1} = 'TW6dB';
bsmflayerdata{end+1} = 'SBMF';
sw3dblayerdata{end+1} = 'SW3dB';
sw6dblayerdata{end+1} = 'SW6dB';


plot_the_mtf_data_layer_comparisons(btmflayerdata, exp);
plot_the_mtf_data_layer_comparisons(tw3dblayerdata, exp);
plot_the_mtf_data_layer_comparisons(tw6dblayerdata, exp);
plot_the_mtf_data_layer_comparisons(bsmflayerdata, exp);
plot_the_mtf_data_layer_comparisons(sw3dblayerdata, exp);
plot_the_mtf_data_layer_comparisons(sw6dblayerdata, exp);


function plot_the_mtf_data_layer_comparisons(data, exp)

figure;

% ---------------------------------------------------
%            Plot the raw layer data
% ---------------------------------------------------

ylabeltext{1} = 'Mean filter';
ylabeltext{2} = 'Mean filter smoothed';
ylabeltext{3} = 'Max filter';
ylabeltext{4} = 'Max filter smoothed';
ylabeltext{5} = 'Min filter';
ylabeltext{6} = 'Min filter smoothed';

colortype{1} = 'b';
colortype{2} = 'g';
colortype{3} = 'r';
colortype{4} = 'm';
colortype{5} = 'c';
colortype{6} = 'k';

markerlinetype{1} = 's-';
markerlinetype{2} = 's-';
markerlinetype{3} = 's-';
markerlinetype{4} = 's:';
markerlinetype{5} = 's-';
markerlinetype{6} = 's-';


for i = 1:6  % 4 different smoothing procedures used

   subplot(6,3,3*i-2);
   hold on;

   for j = 1:6  % layers 1 thru 6

      if ( ~isempty( eval(['data{' num2str(j) '}{1}']) ) )
         plot(eval(['data{' num2str(j) '}{1}']), ...
              eval(['data{' num2str(j) '}{' num2str(i+1) '}']), ...
              [colortype{j} markerlinetype{j}], ...
              'markerfacecolor', colortype{j}, 'markersize', 2, 'linewidth', 1);
      else
         plot(-10, -100, [colortype{j} markerlinetype{j}], ...
              'markerfacecolor', colortype{j}, 'markersize', 2, 'linewidth', 1);
      end

   end % (for j)

   ylabel(ylabeltext{i});
   set(gca,'xtick', -1:6, 'xticklabel', -1:6);
   set(gca,'fontsize', 8);

   if ( strcmp(data{end},'TBMF') | strcmp(data{end},'TW3dB') | strcmp(data{end},'TW6dB'))
      set(gca,'ytick', 0:10:40, 'yticklabel', 0:10:40);
      axis([-1 6 -5 45]);
   elseif ( strcmp(data{end},'SBMF') | strcmp(data{end},'SW3dB') | strcmp(data{end},'SW6dB'))
      set(gca,'ytick', 0:1:4, 'yticklabel', 0:1:4, 'tickdir', 'out');
      axis([-1 6 -1 5]);
   else
      error('Bad data given to plot_the_mtf_data_layer_comparisons.m');
   end

   if ( i == 6 )
      xlabel('Ventral -- > Dorsal (mm)');
      legend('L1', 'L2', 'L3', 'L4', 'L5', 'L6', 0);
   elseif ( i == 1 )
      title('Smoothing Procedures');
      set(gca, 'xticklabel', '');
   else
      set(gca, 'xticklabel', '');
   end

end % (for i)


% ---------------------------------------------------
%   Plot the layer comparisons - normalization
% ---------------------------------------------------

colortype{1} = 'b';
colortype{2} = 'g';
colortype{3} = 'r';
colortype{4} = 'c';
colortype{5} = 'k';

markerlinetype{1} = 's-';
markerlinetype{2} = 's-';
markerlinetype{3} = 's-';
markerlinetype{4} = 's-';
markerlinetype{5} = 's-';

layers = [1 2 3 5 6];

for i = 1:6  % 4 different smoothing procedures used

   subplot(6,3,3*i-1);
   hold on;

   for j = 1:length(layers)  % layers 1 thru 6

      layerpos = eval(['data{' num2str(layers(j)) '}{1}']);
      layer4pos = eval(['data{4}{1}']);

      [pos, ilayer, i4] = intersect(layerpos, layer4pos);

      if ( ~isempty( pos ) )
         resplayer = eval(['data{' num2str(layers(j)) '}{' num2str(i+1) '}(ilayer)']);
         resp4 = eval(['data{' num2str(4) '}{' num2str(i+1) '}(i4)']);
         plot(pos, resplayer ./ (resp4+eps), ...
              [colortype{j} markerlinetype{j}], ...
              'markerfacecolor', colortype{j}, 'markersize', 2, 'linewidth', 1);
      else
         plot(-10, -100, [colortype{j} markerlinetype{j}], ...
              'markerfacecolor', colortype{j}, 'markersize', 2, 'linewidth', 1);
      end

   end % (for j)

   set(gca,'xtick', -1:6, 'xticklabel', -1:6);
   set(gca,'fontsize', 8);
   axis([-1 6 0 3]);

   if ( i == 3 )
      lh = legend('1/4', '2/4', '3/4', '5/4', '6/4', 1);
      set(gca, 'xticklabel', '');
   elseif ( i == 1 )
      title('Ratios relative to Layer 4');
      set(gca, 'xticklabel', '');
   elseif ( i ~= 6 )
      set(gca, 'xticklabel', '');
   end

end % (for i)


% ---------------------------------------------------
%   Plot the layer comparisons - differences
% ---------------------------------------------------

colortype{1} = 'b';
colortype{2} = 'g';
colortype{3} = 'r';
colortype{4} = 'c';
colortype{5} = 'k';

markerlinetype{1} = 's-';
markerlinetype{2} = 's-';
markerlinetype{3} = 's-';
markerlinetype{4} = 's-';
markerlinetype{5} = 's-';

layers = [1 2 3 5 6];

for i = 1:6  % 4 different smoothing procedures used

   subplot(6,3,3*i);
   hold on;

   for j = 1:length(layers)  % layers 1 thru 6

      layerpos = eval(['data{' num2str(layers(j)) '}{1}']);
      layer4pos = eval(['data{4}{1}']);

      [pos, ilayer, i4] = intersect(layerpos, layer4pos);

      if ( ~isempty( pos ) )
         resplayer = eval(['data{' num2str(layers(j)) '}{' num2str(i+1) '}(ilayer)']);
         resp4 = eval(['data{' num2str(4) '}{' num2str(i+1) '}(i4)']);
         plot(pos, abs(resplayer-resp4) ./ (resp4+eps) * 100, ...
              [colortype{j} markerlinetype{j}], ...
              'markerfacecolor', colortype{j}, 'markersize', 2, 'linewidth', 1);

      else
         plot(-10, -100, [colortype{j} markerlinetype{j}], ...
              'markerfacecolor', colortype{j}, 'markersize', 2, 'linewidth', 1);
      end

   end % (for j)

   set(gca,'xtick', -1:6, 'xticklabel', -1:6);
   set(gca,'fontsize', 8);
   axis([-1 6 -10 150]);
   
   if ( i == 3 )
      lh = legend('(1-4)/4', '(2-4)/4', '(3-4)/4', '(5-4)/4', '(6-4)/4', 1);
      set(gca, 'xticklabel', '');
   elseif ( i == 1 )
      title('Percent difference relative to Layer 4');
      set(gca, 'xticklabel', '');
   elseif ( i ~= 6 )
      set(gca, 'xticklabel', '');
   end

end % (for i)


orient tall;

if ( nargin == 2 )
   if ( isstr(exp) )
      suptitle([exp '  MTF ' data{end} ' Layer Comparisons']);
   end
else
   suptitle(['MTF ' data{end} ' Layer Comparisons']);
end










