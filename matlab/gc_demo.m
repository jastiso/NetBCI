%% Demo for generating covariance based GC FC matrices on single trials

% parameters
st = 3; % task start
en = 6; % in seconds, the feedback period: 3-6s
bl_st = 0.5; % intertrial interval - baseline start
bl_en = 2; % baseline end
% GC params
dt = 100; % window
lag = 10; % lag
t0 = dt+lag;

% split up data into baseline and task datasets
bl = freq; % freq is the output of ft_freqanalysis
bl_idx = find(freq.time > bl_st & freq.time <= bl_en);
bl.time = freq.time(bl_idx);
bl.powspctrm = freq.powspctrm(:,:,:,bl_idx);
freq_grad = freq;
idx = find(freq.time > st & freq.time <= en);
freq_grad.time = freq.time(idx);
freq_grad.powspctrm = freq.powspctrm(:,:,:,idx);
freq_grad = baseline_correct(freq_grad, bl);

% now get FC
gc_grad = zeros((nNode^2-nNode)/2, nTrial);
for t = 1:nTrial
    tmp = zeros(nNode);
    %freq_grad is the baseline corrected output of ft_freqanalysis, with
    %powspctrm as a field
    curr = squeeze(freq_grad.powspctrm(t,:,:));
    [GC, pairs] = cov_GC(curr,dt,lag,t0); % function from andrea
    GC = sum(GC,2);
    % plot
    tmp(logical(tril(ones(nNode),-1))) = GC;
    tmp = tmp + tmp';
    figure(5); clf
    imagesc(tmp); colorbar; pause(0.001)

    % save
    gc_grad(:,t) = GC;
end