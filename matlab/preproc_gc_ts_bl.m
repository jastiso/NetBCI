%% Preprocessing

% preprocess MEG data, and get into connectivity matrix
% @author JStiso jeni.stio@gmail.com Nov 2017



%% Define Global Variables

addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830')
addpath(genpath('/Users/stiso/Documents/MATLAB/eeglab_current'))

top_dir = '/Users/stiso/Documents/MATLAB/NetBCI/';
data_dir = '/Users/stiso/Resilio Sync/NETBCI.RAW/DataBase/1_Signals/2_Segmentation/2_MEG/';
raw_dir = '/Users/stiso/Resilio Sync/NETBCI.RAW/DataBase/1_Signals/0_RawData/2_MEG/';
sessions = [{'Session1'}, {'Session2'}, {'Session3'}, {'Session4'}];
condition = [{'test01'}, {'test02'}, {'test03'}, {'test04'}, {'test05'}, {'test06'}]; % not including rest
subjs = [13];
nNode = 102;
nChunks = 1; %number of sets to divide trials into
freqs = [7,14;15,30;31,45;55,70];
bands = [{'alpha'},{'beta'},{'low_gamma'}, {'gamma'}];
st = 3;
en = 6; % in seconds, the feedback period: 3-6s
num_neg_grad = [];
num_neg_mag = [];
cnte = 1;
% GC params
dt = 100;
lag = 10;
t0 = dt+lag;
% time shift stats
ts_mu = 0;
ts_sigma = 100;

%% Start Loop

for i = 1:numel(sessions)
    session = sessions{i};
    for j = 1:numel(condition)
        cond = condition{j};
        for k = subjs
            subj = k;
            ext = sprintf('%03d',subj);
            % load data
            d = dir([data_dir, session, '/*', cond, '*', '/Seg_MEG_Subj', ext, '_Ses', session(end), '_', cond, '.mat']);
            d_raw = dir([raw_dir, session, '/*', cond, '*', '/RawData_MEG_Subj', ext, '_Ses', session(end), '_', cond, '.mat']);
            try
                load([d.folder, '/', d.name]);
                img_dir = [top_dir, session, '/', cond, '/', ext, '/diagnostics/ts/'];
                save_dir = [top_dir, session, '/', cond, '/', ext, '/FCmatrices/ts/'];
                try
                    load([d_raw.folder, '/', d_raw.name]);
                catch
                    d_raw = dir([raw_dir, 'Session1/*', cond, '*', '/RawData_MEG_Subj', ext, '_Ses1_', cond, '.mat']);
                    load([d_raw.folder, '/', d_raw.name]);
                end
                try
                    eval(['data_bl = ', d.name(1:end-4), '.baseline;']);
                    eval(['data = ', d.name(1:end-4), '.MI;']);
                    eval(['clear ', d.name(1:end-4)]);
                    
                catch
                    eval(['data_bl = ', d.name(1:13), '16', d.name(16:(end-4)), '.baseline;']);
                    eval(['data = ', d.name(1:13), '16', d.name(16:(end-4)), '.MI;']);
                    eval(['clear ', d.name(1:13), '16', d.name(16:(end-4))]);
                end
                try
                    eval(['data_raw = ', d_raw.name(1:end-4), ';']);
                    eval(['clear ', d_raw.name(1:end-4)]);
                catch
                    eval(['data_raw = ', d_raw.name(1:17), '16', d_raw.name(20:(end-4)), ';']);
                    eval(['clear ', d_raw.name(1:17), '16', d_raw.name(20:(end-4))]);
                end
                
                % make directories
                if ~exist([img_dir, '/grad/'], 'dir')
                    mkdir([img_dir, '/grad/'])
                end
                % make directories
                if ~exist([img_dir, '/mag/'], 'dir')
                    mkdir([img_dir, '/mag/'])
                end
                % make directories
                if ~exist([save_dir, '/grad/'], 'dir')
                    mkdir([save_dir, '/grad/'])
                end
                % make directories
                if ~exist([save_dir, '/mag/'], 'dir')
                    mkdir([save_dir, '/mag/'])
                end
                
                 % remove trials with NaNs - this is for subj 13
                if numel(data.trial) == numel(data_bl.trial)
                    for m = 1:numel(data.trial)
                        if any(any(isnan(data.trial{m}))) || any(any(isnan(data_bl.trial{m})))
                            data.trial = {data.trial{[1:(m-1), (m+1):end]}};
                            data_bl.trial = {data_bl.trial{[1:(m-1), (m+1):end]}};
                            
                            data.time = {data.time{[1:(m-1), (m+1):end]}};
                            data_bl.time = {data_bl.time{[1:(m-1), (m+1):end]}};
                            
                            data.sampleinfo = data.sampleinfo([1:(m-1), (m+1):end],:);
                            data_bl.sampleinfo = data_bl.sampleinfo([1:(m-1), (m+1):end],:);
                        end
                    end
                end
                
                
                % view the data (optional)
                cfg = [];
                cfg.viewmode = 'vertical';
                %ft_databrowser(cfg, data)
                
                % grab some more variables
                trl_len = numel(data.trial{1}(1,:));
                srate = data.fsample;
                nTrial_o = numel(data.trial);
                
                % check fft - will fix later if there is noise
                figure(1); clf
                spectopo(data.trial{1}, 0, round(data.fsample));
                pause(.1)
                saveas(gca, [img_dir, 'fft_before.png'], 'png')
                
                % separate into magnetometers and gradiometers
                data.senstype = 'neuromag306';
                mag_idx = strcmp(data_raw.meg_sensors.chantype,'megmag');
                grad_idx = ~mag_idx;
                data_mag = data;
                data_grad = data;
                % grad
                data_grad.label = data_grad.label(grad_idx);
                for m = 1:numel(data_grad.trial)
                    data_grad.trial{m} = data_grad.trial{m}(grad_idx,:);
                end
                
                % now mag
                data_mag.label = data_mag.label(mag_idx);
                for m = 1:numel(data_mag.trial)
                    data_mag.trial{m} = data_mag.trial{m}(mag_idx,:);
                end
                
                
                % timeshift (only for task)
                for t = 1:nTrial_o
                    data_mag1.trial{t} = time_shift(data_mag.trial{t}, ...
                        (ts_mu*srate)/1000, round((ts_sigma*srate)/1000), 'mag', data_mag.label);
                    data_grad1.trial{t} = time_shift(data_grad.trial{t}, ...
                        (ts_mu*srate)/1000, round((ts_sigma*srate)/1000), 'grad', data_grad.label);
                end
                
                
                
                % now for BL
                data_mag_bl = data_bl;
                data_grad_bl = data_bl;
                data_grad_bl.label = data_grad_bl.label(grad_idx);
                % grad
                for m = 1:numel(data_grad_bl.trial)
                    data_grad_bl.trial{m} = data_grad_bl.trial{m}(grad_idx,:);
                end
                
                % mag
                data_mag_bl.label = data_mag_bl.label(mag_idx);
                for m = 1:numel(data_mag_bl.trial)
                    data_mag_bl.trial{m} = data_mag_bl.trial{m}(mag_idx,:);
                end
                
                %% Make connectivity matrices
                
                % get power in task and BL
                % use multitapers for gamma, and single tapers for low
                % freq
                
                % magnetometers
                [freq_mag_mt, freq_mag_st] = NetBCI_get_pow(data_mag, st, en);
                [freq_mag_bl_mt, freq_mag_bl_st] = NetBCI_get_pow(data_mag_bl, st, en);
                %gradiometers
                [freq_grad_mt, freq_grad_st] = NetBCI_get_pow(data_grad, st, en);
                [freq_grad_bl_mt, freq_grad_bl_st] = NetBCI_get_pow(data_grad_bl, st, en);
                
                
                % plot
                figure(1); clf
                imagesc(squeeze(freq_mag_mt.powspctrm(end,94,:,:)))
                colorbar
                pause(0.01)
                figure(2); clf
                imagesc(squeeze(freq_mag_bl_mt.powspctrm(end,94,:,:)))
                colorbar
                pause(0.01)
                figure(3); clf
                imagesc(squeeze(freq_mag_st.powspctrm(end,94,:,:)))
                colorbar
                pause(0.01)
                figure(4); clf
                imagesc(squeeze(freq_mag_bl_st.powspctrm(end,94,:,:)))
                colorbar
                pause(0.01)
                
                for f = 1:numel(bands)
                    f_range = freqs(f,:);
                    curr_band = bands{f};
                    % average to get current band
                    if contains(curr_band, 'gamma')
                        % magnetometers
                        % initialize
                        freq_mag = freq_mag_mt;
                        freq_mag_bl = freq_mag_bl_mt;
                        % get freq range
                        idx = freq_mag_mt.freq >= f_range(1) & freq_mag_mt.freq < f_range(end);
                        % average
                        pow = log10(freq_mag_mt.powspctrm(:,:,idx,:));
                        pow = squeeze(nanmean(pow,3));
                        freq_mag.powspctrm = pow;
                        freq_mag.freq = squeeze(mean(freq_mag.freq(idx)));
                        %average bl
                        pow = log10(freq_mag_bl_mt.powspctrm(:,:,idx,:));
                        pow = squeeze(nanmean(pow,3));
                        freq_mag_bl.powspctrm = pow;
                        freq_mag_bl.freq = squeeze(mean(freq_mag_bl.freq(idx)));
                        
                        % gradiometers
                        % initialize
                        freq_grad = freq_grad_mt;
                        freq_grad_bl = freq_grad_bl_mt;
                        % get freq range
                        idx = freq_grad_mt.freq >= f_range(1) & freq_grad_mt.freq < f_range(end);
                        % average
                        pow = log10(freq_grad_mt.powspctrm(:,:,idx,:));
                        pow = squeeze(nanmean(pow,3));
                        freq_grad.powspctrm = pow;
                        freq_grad.freq = squeeze(mean(freq_grad.freq(idx)));
                        %average bl
                        pow = log10(freq_grad_bl_mt.powspctrm(:,:,idx,:));
                        pow = squeeze(nanmean(pow,3));
                        freq_grad_bl.powspctrm = pow;
                        freq_grad_bl.freq = squeeze(mean(freq_grad_bl.freq(idx)));
                    else
                        % magnetometers
                        % initialize
                        freq_mag = freq_mag_st;
                        freq_mag_bl = freq_mag_bl_st;
                        % get freq range
                        idx = freq_mag_st.freq >= f_range(1) & freq_mag_st.freq < f_range(end);
                        % average
                        pow = log10(freq_mag_st.powspctrm(:,:,idx,:));
                        pow = squeeze(nanmean(pow,3));
                        freq_mag.powspctrm = pow;
                        freq_mag.freq = squeeze(mean(freq_mag.freq(idx)));
                        %average bl
                        pow = log10(freq_mag_bl_st.powspctrm(:,:,idx,:));
                        pow = squeeze(nanmean(pow,3));
                        freq_mag_bl.powspctrm = pow;
                        freq_mag_bl.freq = squeeze(mean(freq_mag_bl.freq(idx)));
                        
                        % gradiometers
                        % initialize
                        freq_grad = freq_grad_st;
                        freq_grad_bl = freq_grad_bl_st;
                        % get freq range
                        idx = freq_grad_st.freq >= f_range(1) & freq_grad_st.freq < f_range(end);
                        % average
                        pow = log10(freq_grad_st.powspctrm(:,:,idx,:));
                        pow = squeeze(nanmean(pow,3));
                        freq_grad.powspctrm = pow;
                        freq_grad.freq = squeeze(mean(freq_grad.freq(idx)));
                        %average bl
                        pow = log10(freq_grad_bl_st.powspctrm(:,:,idx,:));
                        pow = squeeze(nanmean(pow,3));
                        freq_grad_bl.powspctrm = pow;
                        freq_grad_bl.freq = squeeze(mean(freq_grad_bl.freq(idx)));
                    end
                    
                    % combine gradiometers
                    freq_grad = combine_gradiometers(freq_grad);
                    freq_grad_bl = combine_gradiometers(freq_grad_bl);

                    % plot
                    figure(1); clf
                    imagesc(squeeze(freq_grad.powspctrm(end,:,:))); colorbar
                    saveas(gca, [img_dir, bands{f}, 'pow.png'], 'png')
                    
                    figure(2); clf
                    imagesc(squeeze(freq_grad_bl.powspctrm(end,:,:))); colorbar
                    saveas(gca, [img_dir, bands{f}, 'pow_bl.png'], 'png')
                    
                    
                    % baseline correct power (z-score)
                    % mag
                    freq_mag = baseline_correct(freq_mag, freq_mag_bl);
                    % grad
                    freq_grad = baseline_correct(freq_grad, freq_grad_bl);
                    nTrial = size(freq_mag.powspctrm,1);

                                        
                    figure(3); clf
                    imagesc(squeeze(freq_grad.powspctrm(end,:,:))); colorbar
                    caxis([-2 5])
                    saveas(gca, [img_dir, bands{f}, 'pow_z.png'], 'png')
                    
                    
                    
                    % get cov based GC
                    gc_mag = zeros((nNode^2-nNode)/2, nTrial);
                    gc_grad = zeros((nNode^2-nNode)/2, nTrial);
                    for t = 1:nTrial
                        tmp = zeros(nNode);
                        % magnetometers
                        curr = squeeze(freq_mag.powspctrm(t,:,:));
                        [GC, pairs] = cov_GC(curr,dt,lag,t0);
                        GC = sum(GC,2);
                        % plot
                        tmp(logical(tril(ones(nNode),-1))) = GC;
                        tmp = tmp + tmp';
                        figure(5); clf
                        imagesc(tmp); colorbar; pause(0.001)
                        if t == nTrial
                            saveas(gca, [img_dir, bands{f}, 'gc_mag.png'], 'png')
                        end
                        % save
                        gc_mag(:,t) = GC;
                        
                        %gradiometers
                        tmp = zeros(nNode);
                        curr = squeeze(freq_grad.powspctrm(t,:,:));
                        [GC, pairs] = cov_GC(curr,dt,lag,t0);
                        GC = sum(GC,2);
                        % plot
                        tmp(logical(tril(ones(nNode),-1))) = GC;
                        tmp = tmp + tmp';
                        figure(5); clf
                        imagesc(tmp); colorbar; pause(0.001)
                        if t == nTrial
                            saveas(gca, [img_dir, bands{f}, 'gc_grad.png'], 'png')
                        end
                        % save
                        gc_grad(:,t) = GC;
                    end
                    save([save_dir, 'NMF_', curr_band, '_mag_gc.mat'], 'gc_mag')
                    save([save_dir, 'NMF_', curr_band, '_grad_gc.mat'], 'gc_grad')
                end
            catch
                errors{cnte,1} = [d.folder, '/', d.name];
                cnte = cnte + 1;
            end
        end
    end
end







