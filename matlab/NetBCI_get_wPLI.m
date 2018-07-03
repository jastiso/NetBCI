function [ conn ] = NetBCI_get_wPLI( data, trials, sens )
%Get wPLIT from fieldtrip data
% needs data, and which trials to use, and whether or not these are
% gradiometers which need to be combined

% field trip mtmfft method
cfg.output     = 'fourier';
cfg.method     = 'mtmfft';
cfg.foilim     = freqs(f,:);
cfg.tapsmofrq  = 5;
cfg.keeptrials = 'yes';
cfg.trials = trials;
freq    = ft_freqanalysis(cfg, data);
freq.freq = mean(freq.freq);
freq.fourierspctrm = mean(freq.fourierspctrm,3);

if strcmp(lower(sens), 'grad')
    freq = combine_gradiometers(freq);
end

cfg = [];
cfg.method = 'wpli_debiased';
cfg.keeptrials = 'yes';
%cfg.complex = 'wpli';
wpli = ft_connectivityanalysis(cfg, freq);
conn = wpli;


end

