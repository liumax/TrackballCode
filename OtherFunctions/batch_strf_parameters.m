function batch_strf_parameters(sitenum)
% batch_strf_parameters - compute all strf parameters,
% such as energy, phase-locking index, svd
% singular values, in one big process.
%
% batch_strf_parameters(sitenum)
%
% sitenum : optional. If specified it tells which site(s)
%   to process. With no input all sites in a given experimental
%   directory are processed.
%
% caa 9/16/03

if ( nargin == 0 )
   sitenum = [];
end


if ( isempty(sitenum) )
   dsite = dir('site*');
else
   for i = 1:length(sitenum)
      dsite(i).name = ['site' num2str(sitenum(i))];
   end
end


stim{1} = 'dmr1';
stim{2} = 'dmr2';

for i = 1:length(dsite)

   for k = 1:length(stim)

      site = dsite(i).name;
      filename = [site '\*-' stim{k} '-*-strfcmb.mat'];
      dfile = dir(filename);

      if ( ~isempty(dfile) )

         infile = dfile.name;
         indstrfcmb = findstr(infile, 'strfcmb');

         if ( strcmp(stim{k}, 'dmr1') | strcmp(stim{k}, 'dmr2') | strcmp(stim{k}, 'rn1') )

            outfile = [dfile.name(1:indstrfcmb-1) 'strf-params2'];
            fprintf('%s\n', outfile);

            load([site '\' infile]);

            if ( exist('strf')==1 & exist('trigger')==1 & exist('cmb')==1 )
 
               params = strf_parameters(strf, trigger);

               save([site '\' outfile], 'params');

            end % (if)

            clear params trigger strf cmb  % get rid of workspace variables

         end % (if)

      end % (if)

   end % (for k)

end % (for i)





