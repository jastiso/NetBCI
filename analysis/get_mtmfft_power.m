function [ pow ] = get_mtmfft_power( data, srate, st, en, freqs,t)
%% Get power

% getting power, using mtmfft

%@author JStiso jeni.stiso@gmail.com

% Change Log
% April 8th - made script


%% get cross spectral density
% The product of the Fourier transform of one signal with the conjugate of the Fourier Transform of the signal gives the Cross-Spectral density (CSD).

% this will be estimated with wavelet coherence, doesn't assume
% stationarity, and gets you lots of time points
% defining constants

cfg = [];
cfg.method     = 'mtmfft';
cfg.tapsmofrq  = 0.5;
cfg.output     = 'pow';
cfg.foi        = freqs;
cfg.toi        = st:1/srate:en;
cfg.pad        = 'nextpow2';
cfg.trials     = t;
wave = ft_freqanalysis(cfg, data);

% avg, and log transform
pow = squeeze(mean(log(wave.powspctrm),2));

end

