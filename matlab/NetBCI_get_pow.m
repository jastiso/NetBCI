function [ freq_mt, freq_st ] = NetBCI_get_pow( data, st, en, bl_st, bl_en)
% get multitaper power for HF, and single tapers for LF
% Change log
% June 21, 018 - Changed to incorporate baseline corrections

% multitaper
cfg = [];
cfg.output     = 'pow';
cfg.method     = 'mtmconvol';
cfg.keeptrials = 'yes';
cfg.foi        = 31:2:70;
cfg.t_ftimwin  = 5./cfg.foi;
cfg.tapsmofrq  = 0.4 *cfg.foi;
cfg.toi        = 0:.004:en;
freq = ft_freqanalysis(cfg, data);
% baseline correct
bl = freq;
bl_idx = find(freq.time > bl_st & freq.time <= bl_en);
bl.time = freq.time(bl_idx);
bl.powspctrm = freq.powspctrm(:,:,:,bl_idx);
freq_mt = freq;
idx = find(freq.time > st & freq.time <= en);
freq_mt.time = freq.time(idx);
freq_mt.powspctrm = freq.powspctrm(:,:,:,idx);
freq_mt = baseline_correct(freq_mt, bl);
% single tapers
cfg              = [];
cfg.output       = 'pow';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.keeptrials = 'yes';
cfg.foi          = 2:2:30;
cfg.t_ftimwin    = 7./cfg.foi;  % 7 cycles per time window
cfg.toi          = 0:.004:en;
freq = ft_freqanalysis(cfg, data);
% baseline correct
bl = freq;
bl_idx = find(freq.time > bl_st & freq.time <= bl_en);
bl.time = freq.time(bl_idx);
bl.powspctrm = freq.powspctrm(:,:,:,bl_idx);
freq_st = freq;
idx = find(freq_st.time > st & freq_st.time <= en);
freq_st.time = freq.time(idx);
freq_st.powspctrm = freq.powspctrm(:,:,:,idx);
freq_st = baseline_correct(freq_st, bl);
end

