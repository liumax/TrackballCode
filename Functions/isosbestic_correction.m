function photoSig = isosbestic_correction(data)
% performs isosbestic correction and returns photoSig with dF/F

uncorrSig = double(data.streams.x70G.data'); % uncorrected signal
isoSig = double(data.streams.x05G.data'); % isosbestic signal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REMOVE OUTLIERS TO GET BEST FIT: (SO I DON'T FIT SHIT THAT HAPPENS WHEN
% TDT IS TURNED ON/OFF:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z_sig = zscore(uncorrSig);
z_iso = zscore(isoSig);
a = find(z_sig<-10);
b = find(z_iso<-10);
rmIdx = union(a,b);
caSig_rm = uncorrSig; caSig_rm(rmIdx) = []; % Removing signal >10std below mean
iso_rm = isoSig; iso_rm(rmIdx) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIT SIGNAL AND GET dF/F
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b = glmfit(iso_rm,caSig_rm);
iso_scaled = b(1) + isoSig*b(2);
photoSig = (uncorrSig-iso_scaled)./iso_scaled;
photoSig(rmIdx) = 0; % Replace outliers with 0

