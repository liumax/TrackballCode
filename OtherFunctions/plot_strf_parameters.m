function plot_strf_parameters(params)
%plot_strf_parameters - plots best temporal, spectral modulation values
%   along with asymmetry indices and separability indices.
%
%   caa 10/28/02

exp = params(1).exp;
site = params(1).site;
stim = params(1).stim;

layer2 = 200;
layer3 = 400;
layer4 = 800;
layer5 = 1100;
layer6 = 1500;
whitematter = 2000;

layer1label = 100;
layer2label = 300;
layer3label = 600;
layer4label = 950;
layer5label = 1300;
layer6label = 1750;
whitematterlabel = 2200;

% initialize the variables
position = zeros(1,length(params));
energy = zeros(1,length(params));
pli = zeros(1,length(params));
tci5 = zeros(1,length(params));
tci10 = zeros(1,length(params));
tci15 = zeros(1,length(params));
scip25 = zeros(1,length(params));
scin25 = zeros(1,length(params));
sepindex = zeros(1,length(params));
sepindex95 = zeros(1,length(params));

tmf = params(1).tmf;
tmf = tmf(:)';
ind0 = find(tmf==0);
tmf_right = tmf(ind0:end);

xmf = params(1).xmf;
xmf = xmf(:)';


% First get all the separability indices

for i = 1:length(params)

   position(i) = params(i).position;

   energy(i) = params(i).energy(end,2);
   pli(i) = params(i).pli;
   tci5(i) = params(i).tci(1,2);
   tci10(i) = params(i).tci(2,2);
   tci15(i) = params(i).tci(3,2);
   scip25(i) = params(i).sci(8,2);
   scin25(i) = params(i).sci(9,2);

   eigvals = params(i).eigvals;
   sepindex(i) = 1 - eigvals(1) / (sum(eigvals)+eps);

   % Can use the following lines to have a power cutoff for eigvals
   percent_cutoff = 95;
   eigpercent = 100 * cumsum(eigvals) / sum(eigvals);
   i95percent = find(eigpercent < percent_cutoff);

   sepindex95(i) = 1 - eigvals(1) / (sum(eigvals(i95percent))+eps);

   tmtf = sum(params(i).rtf,1);
   tmtf = tmtf(:)'
   tmtf_left = tmtf(1:ind0-1)
   tmtf_right = tmtf(ind0:end)
   tmtf = tmtf_right + [0 fliplr(tmtf_left)];

   [flipud([tmtf_left 0]') tmtf_right' tmtf']

   xmtf = sum(params(i).rtf,2)';
   xmtf = xmtf(:)';

   f1 = figure;
   subplot(3,1,1);
   imagesc(params(i).rtf);
   axis xy;
   subplot(3,1,2);
   plot(tmf_right, tmtf, 'ko-');
   subplot(3,1,3);
   plot(xmf, xmtf, 'ko-');
   pause;

end

figure;
suptitle(sprintf('%s site%.0f %s', exp, site, stim));

% First plot the separability index vs. depth
subplot(2,2,1);
hold on;
plot(sepindex, position, 'ko', 'markerfacecolor','k');
% plot(sepindex95, position, 'kv', 'markerfacecolor','k');
box on;
axis([0 1 0 2400]);
plot([0 1],[layer2 layer2],'k:');
plot([0 1],[layer3 layer3],'k:');
plot([0 1],[layer4 layer4],'k:');
plot([0 1],[layer5 layer5],'k:');
plot([0 1],[layer6 layer6],'k:');
plot([0 1],[whitematter whitematter],'k-');
text([0.9],[layer1label],'L I');
text([0.9],[layer2label],'L II');
text([0.9],[layer3label],'L III');
text([0.9],[layer4label],'L IV');
text([0.9],[layer5label],'L V');
text([0.9],[layer6label],'L VI');
text([0.9],[whitematterlabel],'WM');

set(gca,'ytick',[0 layer2 layer3 layer4 layer5 layer6 whitematter]);
xlabel('SI');
ylabel('Depth (um)');
title('Separability Index');
set(gca,'ydir','rev');


% plot the strf energy vs. depth
subplot(2,2,2);
hold on;
plot(energy, position, 'ko', 'markerfacecolor','k');
box on;
axis([0 1 0 2400]);
plot([0 1],[layer2 layer2],'k:');
plot([0 1],[layer3 layer3],'k:');
plot([0 1],[layer4 layer4],'k:');
plot([0 1],[layer5 layer5],'k:');
plot([0 1],[layer6 layer6],'k:');
plot([0 1],[whitematter whitematter],'k-');
text([0.9],[layer1label],'L I');
text([0.9],[layer2label],'L II');
text([0.9],[layer3label],'L III');
text([0.9],[layer4label],'L IV');
text([0.9],[layer5label],'L V');
text([0.9],[layer6label],'L VI');
text([0.9],[whitematterlabel],'WM');

set(gca,'ytick',[0 layer2 layer3 layer4 layer5 layer6 whitematter]);
xlabel('Energy');
ylabel('Depth (um)');
title('STRF Energy');
set(gca,'ydir','rev');



% plot phase-locking index vs. depth
subplot(2,2,3);
hold on;
plot(pli, position, 'ko', 'markerfacecolor','k');
box on;
axis([0 0.5 0 2400]);
plot([0 1],[layer2 layer2],'k:');
plot([0 1],[layer3 layer3],'k:');
plot([0 1],[layer4 layer4],'k:');
plot([0 1],[layer5 layer5],'k:');
plot([0 1],[layer6 layer6],'k:');
plot([0 1],[whitematter whitematter],'k-');
text([0.45],[layer1label],'L I');
text([0.45],[layer2label],'L II');
text([0.45],[layer3label],'L III');
text([0.45],[layer4label],'L IV');
text([0.45],[layer5label],'L V');
text([0.45],[layer6label],'L VI');
text([0.45],[whitematterlabel],'WM');

set(gca,'ytick',[0 layer2 layer3 layer4 layer5 layer6 whitematter]);
xlabel('PLI');
ylabel('Depth (um)');
title('Phase-locking index');
set(gca,'ydir','rev');



% plot temporal correlation index vs. depth
subplot(2,2,4);
hold on;
plot(tci5, position, 'ko', 'markerfacecolor','k');
plot(tci10, position, 'rs', 'markerfacecolor','r');
plot(tci15, position, 'gv', 'markerfacecolor','g');
legend('5ms', '10ms', '15ms', 0);
box on;
axis([0 1 0 2400]);
plot([0 1],[layer2 layer2],'k:');
plot([0 1],[layer3 layer3],'k:');
plot([0 1],[layer4 layer4],'k:');
plot([0 1],[layer5 layer5],'k:');
plot([0 1],[layer6 layer6],'k:');
plot([0 1],[whitematter whitematter],'k-');
text([0.9],[layer1label],'L I');
text([0.9],[layer2label],'L II');
text([0.9],[layer3label],'L III');
text([0.9],[layer4label],'L IV');
text([0.9],[layer5label],'L V');
text([0.9],[layer6label],'L VI');
text([0.9],[whitematterlabel],'WM');
% plot([0.75], [2175], 'rs', 'markerfacecolor', 'r');
% text([0.78],[2175],'10ms');

set(gca,'ytick',[0 layer2 layer3 layer4 layer5 layer6 whitematter]);
xlabel('TCI');
ylabel('Depth (um)');
title('Temporal correlation index');
set(gca,'ydir','rev');

orient landscape

% print_mfilename(mfilename);



if ( 0 )



suptitle(sprintf('%s site%.0f %s', exp, site, stim));
print_mfilename(mfilename);
legend('All Evals', '95% Evals', 0);

% We also want to plot the 
% (1) phase locking index
% (2) temporal correlation index
% (3) strf energy
% (4) mtfs
% (5) best temporal modulation frequency
% (6) best spectral modulation frequency


pos = [rfstats.position];
bsm = [rfstats.bestsm];
btm = [rfstats.besttm];
bvel = [rfstats.bestvel];
ai = [rfstats.ai];
pval = [rfstats.aipval];
si = [rfstats.si];


subplot(3,2,5);
hold on;
plot(si, pos, 'ko', 'markerfacecolor','k');
box on;
axis([0 1 0 2500]);
plot([0 1],[layer2 layer2],'k:');
plot([0 1],[layer3 layer3],'k:');
plot([0 1],[layer4 layer4],'k:');
plot([0 1],[layer5 layer5],'k:');
plot([0 1],[layer6 layer6],'k:');
xlabel('SI');
ylabel('Depth (um)');
title('Separability Index');
set(gca,'ydir','rev');


subplot(3,2,1);
hold on;
plot(bsm,pos,'ko','markerfacecolor','k');
xlim = get(gca,'xlim');
axis([xlim(1) xlim(2) 0 2500]);
plot([xlim(1) xlim(2)],[layer2 layer2],'k:');
plot([xlim(1) xlim(2)],[layer3 layer3],'k:');
plot([xlim(1) xlim(2)],[layer4 layer4],'k:');
plot([xlim(1) xlim(2)],[layer5 layer5],'k:');
plot([xlim(1) xlim(2)],[layer6 layer6],'k:');
box on;
set(gca,'ydir','rev');
xlabel('Spec Mod (cyc/oct)');
ylabel('Depth (um)');
title('Best Spec Mod');

subplot(3,2,2);
hold on;
plot(btm,pos,'ko','markerfacecolor','k');
xlim = get(gca,'xlim');
temp = max(abs(xlim));
axis([-temp temp 0 2500]);
plot([-temp temp],[layer2 layer2],'k:');
plot([-temp temp],[layer3 layer3],'k:');
plot([-temp temp],[layer4 layer4],'k:');
plot([-temp temp],[layer5 layer5],'k:');
plot([-temp temp],[layer6 layer6],'k:');
box on;
set(gca,'ydir','rev');
xlabel('Temp Mod (cyc/s)');
title('Best Temp Mod');


subplot(3,2,3);
hold on;
plot(bvel,pos,'ko','markerfacecolor','k');
box on;
xlim = get(gca,'xlim');
axis([xlim(1) xlim(2) 0 2500]);
plot([xlim(1) xlim(2)],[layer2 layer2],'k:');
plot([xlim(1) xlim(2)],[layer3 layer3],'k:');
plot([xlim(1) xlim(2)],[layer4 layer4],'k:');
plot([xlim(1) xlim(2)],[layer5 layer5],'k:');
plot([xlim(1) xlim(2)],[layer6 layer6],'k:');
xlabel('Velocity (oct/s)');
ylabel('Depth (um)');
title('Best Velocity');
set(gca,'ydir','rev');

subplot(3,2,4);
hold on;
plot(aisig, possig, 'ko', 'markerfacecolor','k');
plot(ainotsig, posnotsig, 'ko', 'markerfacecolor','none');
box on;
maxai = max([abs(min(ai)) max(ai)]);
xlo = -maxai-0.05;
xhi = maxai+0.05;
axis([xlo xhi 0 2500]);
plot([xlo xhi],[layer2 layer2],'k:');
plot([xlo xhi],[layer3 layer3],'k:');
plot([xlo xhi],[layer4 layer4],'k:');
plot([xlo xhi],[layer5 layer5],'k:');
plot([xlo xhi],[layer6 layer6],'k:');
legend('p<=0.05','p>0.05');
xlabel('AI');
title('Asymmetry Index');
set(gca,'ydir','rev');

subplot(3,2,5);
hold on;
plot(si, pos, 'ko', 'markerfacecolor','k');
box on;
axis([0 1 0 2500]);
plot([0 1],[layer2 layer2],'k:');
plot([0 1],[layer3 layer3],'k:');
plot([0 1],[layer4 layer4],'k:');
plot([0 1],[layer5 layer5],'k:');
plot([0 1],[layer6 layer6],'k:');
xlabel('SI');
ylabel('Depth (um)');
title('Separability Index');
set(gca,'ydir','rev');

suptitle(sprintf('%s site%.0f %s', ...
         rfstats(1).exp, rfstats(1).site, rfstats(1).stim));

print_mfilename(mfilename);
orient tall;

end