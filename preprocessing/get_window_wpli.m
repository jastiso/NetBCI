function [ wpli ] = get_window_wpli( data, srate, st, en, freqs, t, grad)
%% Get wpli

% wpli for single trials
% ovlp should be a decimal, 0-1
% st and en are in seconds
% winLength should be ~5*your longest wavelength, in ms

%@author JStiso jeni.stiso@gmail.com

% Change Log
% August 27th, 2018 - created script
% August 29th - updated to calculate wpli over a sliding window
% August 30th - updated to use wavelets


%% get cross spectral density
% The product of the Fourier transform of one signal with the conjugate of the Fourier Transform of the signal gives the Cross-Spectral density (CSD).

% this will be estimated with wavelet coherence, doesn't assume
% stationarity, and gets you lots of time points
% defining constants

cfg = [];
cfg.method     = 'wavelet';
cfg.width      = 6;
cfg.output     = 'powandcsd';
cfg.foi        = freqs;
cfg.toi        = st:1/srate:en;
cfg.pad        = 'nextpow2';
cfg.trials     = t;
wave = ft_freqanalysis(cfg, data);

% housekeeping for formatting things later (basically to move between matrix and vector representations). There might be a better way to
% do this
nSens = size(data.trial{1},1);
nEdge = (nSens^2 - nSens)/2;
cnt = 0;
label = zeros(nEdge,2);
for i = 1:nSens
    for j = (i+1):nSens
        cnt = cnt +1;
        label(cnt,:) = [i,j];
    end
end



%% Get wPLI
% this is code for getting the debiased wPLI - but with 750 timepoints I
% dont think we need to debias it (confirmed this with Javi)
% with the exception of the abs at the end, this is the same as the
% fieldtrip function
% The debiased WPLI-square estimator is computed by (i) computing the
% imaginary components of the cross-spectra (Fig. 7),
% X = imag(wave.crsspctrm);
% % (ii) computing the average imaginary component of the cross-spectra
% % (Fig. 8, top), and
% num = zeros(size(X,1),size(X,2));
% for i = 1:size(X,1)
%     for j = (i+1):size(X,1)
%         num = num + squeeze(X(:,:,i)).*squeeze(X(:,:,i));
%     end
% end
% % (iii) normalizing by the computed average over the magnitudes of the
% % imaginary component of the cross-spectra (Fig. 8, bottom).
% denom = zeros(size(num));
% for i = 1:size(X,1)
%     for j = (i+1):size(X,1)
%         denom = denom + abs(squeeze(X(i,:,:)).*squeeze(X(j,:,:)));
%     end
% end
% wpli_v = num./denom;

% wpli - not debiased
% before abs, this gives the same output as fieldtrip
X = imag(wave.crsspctrm);
wpli_v = abs(nanmean(X,3))./nanmean(abs(X),3);

%average over freqs
wpli_v = abs(mean(wpli_v,2));

% put back into matrix
wpli = zeros(nSens);
for i = 1:nSens
    for j = (i+1):nSens
        row_idx = find(label(:,1)==i);
        col_idx = find(label(:,2)==j);
        idx = intersect(row_idx,col_idx);
        wpli(i,j) = wpli_v(idx);
    end
end

wpli = wpli+wpli';

%% Testing some wavelets

% cfg = [];
% cfg.method     = 'wavelet';
% cfg.width      = 6;
% cfg.output     = 'powandcsd';
% cfg.foi        = 7:1:14;
% cfg.toi        = 3:1/srate:6;
% wave = ft_freqanalysis(cfg, data_grad);

% Check that output is the same as fieldtrip - it is
% cfg = [];
% cfg.method = 'wpli';
% wave.dimord = 'rpt_chancmb_freq';
% wave.crsspctrm = permute(wave.crsspctrm, [3,1,2]);
% wave.powspctrm = permute(wave.powspctrm, [3,1,2]);
% ft = ft_connectivityanalysis(cfg,wave);

end

