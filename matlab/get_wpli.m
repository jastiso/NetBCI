function [ wpli ] = get_wpli( data, st, en, freqs )
%% Get wpli

% wpli for average trials

%@author JStiso jeni.stiso@gmail.com

% Change Log
% August 28th, 2018 - created script


%% Get only relevant time

cfg = [];
cfg.toilim    = [st en];
data = ft_redefinetrial(cfg,data);

%% Fourier spectra
% repeats are over trials

cfg = [];
cfg.method = 'mtmfft';
cfg.taper = 'dpss';
cfg.tapsmofrq = 4;
cfg.keeptrials = 'yes';
cfg.output = 'powandcsd';
cfg.foilim = freqs;
csd = ft_freqanalysis(cfg,data);
% average over freqs
% csd.crsspctrm = mean(csd.crsspctrm,3);


%% combine gradiometers
% not sure what step this should happen in....

%psd_cmb = combine_gradiometers(psd,0); % 0 flags this as PSD


%% Get debiased wPLI
% with the exception of the abs at the end, this is the same as the
% fieldtrip function
% The debiased WPLI-square estimator is computed by (i) computing the 
% imaginary components of the cross-spectra (Fig. 7), 
x = imag(csd.crsspctrm);
% (ii) computing the average imaginary component of the cross-spectra 
% (Fig. 8, top), and 
num = zeros(size(x,2),size(x,3));
for i = 1:size(x,1)
    for j = (i+1):size(x,1)
        num = num + squeeze(x(i,:,:)).*squeeze(x(j,:,:));
    end
end
% (iii) normalizing by the computed average over the magnitudes of the 
% imaginary component of the cross-spectra (Fig. 8, bottom).
denom = zeros(size(num));
for i = 1:size(x,1)
    for j = (i+1):size(x,1)
        denom = denom + abs(squeeze(x(i,:,:)).*squeeze(x(j,:,:)));
    end
end
wpli_v = num./denom;

%average over freqs
wpli_v = abs(mean(wpli_v,2));

% reshape into matrix
wpli = zeros(numel(csd.label));
for i = 1:numel(csd.label)
    row = csd.label{i};
    for j = (i+1):numel(csd.label)
        col = csd.label{j};
        row_idx = find(strcmp(csd.labelcmb(:,2),row));
        col_idx = find(strcmp(csd.labelcmb(:,1),col));
        idx = intersect(row_idx,col_idx);
        wpli(i,j) = wpli_v(idx);
    end
end

wpli = wpli+wpli';

% cfg = [];
% cfg.dojack = 1;
% cfg.method = 'wpli_debiased';
% [wpli_ft] = ft_connectivityanalysis(cfg,csd);
%  
 % average over frequencies
% wpli_ft.wpli_debiasedspctrm = mean(wpli_ft.wpli_debiasedspctrm,2);
% wpli. freq = mean(csd.freq);
% wpli.dimord = 'chancmb';


end

