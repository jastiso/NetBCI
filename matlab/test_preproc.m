%% Preprocessing

% preprocess MEG data, and get into connectivity matrix
% @author JStiso jeni.stio@gmail.com Nov 2017



%% Define Global Variables

addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830')
addpath(genpath('/Users/stiso/Documents/MATLAB/eeglab_current'))

top_dir = '/Users/stiso/Documents/MATLAB/NetBCI/';
data_dir = '/Users/stiso/Resilio Sync/NETBCI.RAW/DataBase/1_Signals/2_Segmentation/2_MEG/';
raw_dir = '/Users/stiso/Resilio Sync/NETBCI.RAW/DataBase/1_Signals/0_RawData/2_MEG/';
sessions = [{'Session2'}];
condition = [{'test03'},{'test05'}]; % not including rest
nSubj = 20;
nChunks = 2; %number of sets to divide trials into
freqs = [3,6;7,14;15,30;31,45;55,70];
bands = [{'theta'},{'alpha'},{'beta'},{'low_gamma'}, {'gamma'}];
dur = 3.5; % in seconds, the length of time you want going into each window
num_neg_grad = [];
num_neg_mag = [];
cnte = 1;

%% Start Loop

for i = 1:numel(sessions)
    session = sessions{i};
    for j = 1:numel(condition)
        cond = condition{j};
        for k = 12
            ext = sprintf('%03d',k);
            % load data
            d = dir([data_dir, session, '/*', cond, '*', '/Seg_MEG_Subj', ext, '_Ses', session(end), '_', cond, '.mat']);
            d_raw = dir([raw_dir, session, '/*', cond, '*', '/RawData_MEG_Subj', ext, '_Ses', session(end), '_', cond, '.mat']);
            try
                load([d.folder, '/', d.name]);
                img_dir = [top_dir, session, '/', cond, '/', ext, '/diagnostics/'];
                save_dir = [top_dir, session, '/', cond, '/', ext, '/FCmatrices/'];
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
                    eval(['data_raw = ', d_raw.name(1:end-4), ';']);
                    eval(['clear ', d_raw.name(1:end-4)]);
                catch
                    eval(['data_bl = ', d.name(1:13), '16', d.name(16:(end-4)), '.baseline;']);
                    eval(['data = ', d.name(1:13), '16', d.name(16:(end-4)), '.MI;']);
                    eval(['clear ', d.name(1:13), '16', d.name(16:(end-4))]);
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
                
                % view the data (optional)
                cfg = [];
                cfg.viewmode = 'vertical';
                %ft_databrowser(cfg, data)
                
                % grab some more variables
                trl_len = numel(data.trial{1}(1,:));
                srate = data.fsample;
                nTrial = numel(data.trial);
                
                % check fft
                figure(1); clf
                spectopo(data.trial{1}, 0, round(data.fsample));
                pause(.1)
                saveas(gca, [img_dir, 'fft_before.png'], 'png')
                
                % filter
                % high pass filter
                [a,b] = butter(4, 1/(data.fsample/2), 'high');
                data.trial{1} = filtfilt(a,b,data.trial{1}')';
                [a,b] = butter(4, 1/(data.fsample/2), 'high');
                data_bl.trial{1} = filtfilt(a,b,data_bl.trial{1}')';
                
                
                % low pass filter
                %                 [a,b] = butter(4, 160/(data.fsample/2), 'low');
                %                 data.trial{1} = filtfilt(a,b,data.trial{1}')';
                %                 [a,b] = butter(4, 160/(data.fsample/2), 'low');
                %                 data_bl.trial{1} = filtfilt(a,b,data_bl.trial{1}')';
                
                
                % image again to see what changed
                figure(1); clf
                spectopo(data.trial{1}, 0, floor(data.fsample));
                pause(.01)
                saveas(gca, [img_dir, 'fft_after.png'], 'png')
                
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
                % now for BL
                data_mag_bl = data;
                data_grad_bl = data;
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
                for f = 1:numel(bands)
                    cnt = 0;
                    conn_tv_mag = zeros((size(data_mag.trial{1},1)^2-size(data_mag.trial{1},1))/2, numel(1:dur*srate:trl_len));
                    conn_tv_grad = zeros((size(data_mag.trial{1},1)^2-size(data_mag.trial{1},1))/2, numel(1:dur*srate:trl_len)); % after combining, grad and nag will be the same size
                    conn_tv_mag_bl = zeros(size(conn_tv_mag));
                    conn_tv_grad_bl = zeros(size(conn_tv_grad));
                    st = 1;
                    for n = 1:nChunks
                        for t  = 1:dur*srate:trl_len
                            cnt = cnt+1;
                            % first magnetometer
                            % filter
                            cfg            = [];
                            cfg.toilim    = [t/srate (t+srate-1)/srate];
                            chunk_mag = ft_redefinetrial(cfg,data_mag);
                            chunk_grad = ft_redefinetrial(cfg,data_grad);
                            chunk_mag_bl = ft_redefinetrial(cfg,data_mag_bl);
                            chunk_grad_bl = ft_redefinetrial(cfg,data_grad_bl);
                            
                            cfg.output     = 'fourier';
                            cfg.method     = 'mtmfft';
                            cfg.foilim     = freqs(f,:);
                            cfg.tapsmofrq  = 5;
                            cfg.keeptrials = 'yes';
                            if n ~= nChunks
                                cfg.trials = st:floor(st+nTrial/nChunks)-1;
                            else
                                cfg.trials = st:nTrial;
                            end
                            freq    = ft_freqanalysis(cfg, chunk_mag);
                            freq.freq = mean(freq.freq);
                            freq.fourierspctrm = mean(freq.fourierspctrm,3);
                            % for bl
                            freq_bl    = ft_freqanalysis(cfg, chunk_mag_bl);
                            freq_bl.freq = mean(freq_bl.freq);
                            freq_bl.fourierspctrm = mean(freq_bl.fourierspctrm,3);
                            
                            
                            cfg = [];
                            cfg.method = 'wpli_debiased';
                            cfg.keeptrials = 'yes';
                            %cfg.complex = 'wpli';
                            wpli = ft_connectivityanalysis(cfg, freq);
                            conn_tv_mag(:,cnt) = reshape(wpli.wpli_debiasedspctrm(logical(triu(ones(size(wpli.wpli_debiasedspctrm)),1))), 1, []);
                            save([save_dir, 'mag/', bands{f}, '_wpli'], 'wpli')
                            wpli_bl = ft_connectivityanalysis(cfg, freq_bl);
                            conn_tv_mag_bl(:,cnt) = reshape(wpli_bl.wpli_debiasedspctrm(logical(triu(ones(size(wpli_bl.wpli_debiasedspctrm)),1))), 1, []);
                            save([save_dir, 'mag/', bands{f}, '_wpli_bl'], 'wpli_bl')
                            
                            %plot
                            figure(1); clf
                            imagesc(wpli.wpli_debiasedspctrm);
                            colorbar
                            saveas(gca, [img_dir, 'mag/', bands{f}, '_wpli.png'], 'png')
                            figure(2); clf
                            imagesc(wpli_bl.wpli_debiasedspctrm);
                            colorbar
                            
                            
                            
                            % now gradiometer
                            % filter
                            cfg            = [];
                            cfg.output     = 'fourier';
                            cfg.method     = 'mtmfft';
                            cfg.foilim     = freqs(f,:);
                            cfg.tapsmofrq  = 5;
                            cfg.keeptrials = 'yes';
                            if n ~= nChunks
                                cfg.trials = st:floor(st+nTrial/nChunks)-1;
                            else
                                cfg.trials = st:nTrial;
                            end
                            freq    = ft_freqanalysis(cfg, chunk_grad);
                            freq.freq = mean(freq.freq);
                            freq.fourierspctrm = mean(freq.fourierspctrm,3);
                            freq_bl    = ft_freqanalysis(cfg, chunk_grad_bl);
                            freq_bl.freq = mean(freq_bl.freq);
                            freq_bl.fourierspctrm = mean(freq_bl.fourierspctrm,3);
                            
                            % combine gradiometers
                            freq = combine_gradiometers(freq);
                            freq_bl = combine_gradiometers(freq_bl);
                            
                            
                            cfg = [];
                            cfg.method = 'wpli_debiased';
                            %cfg.complex = 'wpli';
                            wpli = ft_connectivityanalysis(cfg, freq);
                            conn_tv_grad(:,cnt) = reshape(wpli.wpli_debiasedspctrm(logical(triu(ones(size(wpli.wpli_debiasedspctrm)),1))), 1, []);
                            save([save_dir, 'grad/', bands{f}, '_wpli'], 'wpli')
                            wpli_bl = ft_connectivityanalysis(cfg, freq_bl);
                            conn_tv_grad_bl(:,cnt) = reshape(wpli_bl.wpli_debiasedspctrm(logical(triu(ones(size(wpli_bl.wpli_debiasedspctrm)),1))), 1, []);
                            save([save_dir, 'grad/', bands{f}, '_wpli_bl'], 'wpli_bl')
                            
                            %plot
                            figure(3); clf
                            imagesc(wpli.wpli_debiasedspctrm);
                            colorbar
                            saveas(gca, [img_dir, 'grad/', bands{f}, '_wpli.png'], 'png')
                            figure(4); clf
                            imagesc(wpli_bl.wpli_debiasedspctrm);
                            colorbar
                        end
                        st = st + floor(st+nTrial/nChunks)-1;
                    end
                    num_neg_grad = [num_neg_grad, sum(conn_tv_grad < 0)/numel(conn_tv_grad)];
                    num_neg_mag = [num_neg_mag, sum(conn_tv_mag < 0)/numel(conn_tv_mag)];
                    conn_tv_grad(conn_tv_grad < 0) = 0;
                    conn_tv_mag(conn_tv_mag < 0) = 0;
                    conn_tv_grad_bl(conn_tv_grad_bl < 0) = 0;
                    conn_tv_mag_bl(conn_tv_mag_bl < 0) = 0;
                    save([save_dir, 'NMF_', bands{f}, '_grad'], 'conn_tv_grad');
                    save([save_dir, 'NMF_', bands{f}, '_mag'], 'conn_tv_mag');
                    save([save_dir, 'NMF_', bands{f}, '_grad_bl'], 'conn_tv_grad_bl');
                    save([save_dir, 'NMF_', bands{f}, '_mag_bl'], 'conn_tv_mag_bl');
                end
            catch
                errors{cnte} = [d.folder, '/', d.name];
                cnte = cnte + 1;
            end
        end
    end
end






