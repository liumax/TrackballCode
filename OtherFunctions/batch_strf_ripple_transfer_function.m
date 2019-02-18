function batch_strf_ripple_transfer_function(sitenum)
% batch_strf_ripple_transfer_function - compute all strf parameters,
% such as energy, phase-locking index, svd
% singular values, in one big process.
%
% batch_strf_ripple_transfer_function(sitenum)
%
% sitenum : optional. If specified it tells which site(s)
%   to process. This is done because in some cases
%   different tuning curve parameters were used for
%   different sites. Should have the form [1 4 8 12].
%
% caa 9/16/03

if ( nargin == 0 )
   sitenum = [];
end


if ( isempty(sitenum) )
   dsite = dir('site*');
else
   for i = 1:length(sitenum)
      dsite(i).name = ['site' num2str(sitenum)];
   end
end


stim{1} = 'dmr1';
stim{2} = 'dmr2';
stim{3} = 'rn1';


for i = 1:length(dsite)

   for k = 1:length(stim)

      site = dsite(i).name;
      filename = [site '\*-' stim{k} '-*-strf-params2.mat'];
      dfile = dir(filename);

      if ( ~isempty(dfile) )

         infile = dfile.name;
         indstrfparams = findstr(infile, 'strf-params2');

         if ( strcmp(stim{k}, 'dmr1') | strcmp(stim{k}, 'dmr2') | strcmp(stim{k}, 'rn1') )

            outfile = [dfile.name(1:indstrfparams-1) 'rtf-params2'];
            fprintf('%s\n', outfile);

            load([site '\' infile]);

            if ( exist('params')==1 )
 
               rtf_params = strf_ripple_transfer_function(params);

               save([site '\' outfile], 'rtf_params');

            end % (if)

            clear params rtf_params  % get rid of workspace variables

         end % (if)

      end % (if)

   end % (for k)

end % (for i)



