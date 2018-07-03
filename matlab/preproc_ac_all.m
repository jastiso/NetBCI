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
subjs = [1,17:20];
nNode = 102;
nChunks = 1; %number of sets to divide trials into

% get log spaced intervals
freqs = 3:3:45;
freqs = freqs(~(freqs <= 52 & freqs >= 48)); % don't look at 50 Hz noise
freqs(2,:) = freqs + 2;
freqs = freqs';

st = 3;
en = 6; % in seconds, the feedback period: 3-6s
num_neg_grad = [];
num_neg_mag = [];
cnte = 1;

%% Start Loop

for i = 1:numel(sessions)
    session = sessions{i};
    for j = 1:numel(condition)
        cond = condition{j};
       for k = 1:numel(subjs)
            subj = subjs(k);
            ext = sprintf('%03d',subj);
            
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
                
                % view the data (optional)
                cfg = [];
                cfg.viewmode = 'vertical';
                %ft_databrowser(cfg, data)
                
                % grab some more variables
                trl_len = numel(data.trial{1}(1,:));
                srate = data.fsample;
                nTrial = numel(data.trial);
                
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
                t = (st*srate):(en*srate);
                ac_mag = zeros((size(data_mag.trial{1},1)^2-size(data_mag.trial{1},1))/2, size(freqs,1));
                ac_grad = zeros((size(data_mag.trial{1},1)^2-size(data_mag.trial{1},1))/2, size(freqs,1)); % after combining, grad and nag will be the same size
                ac_mag_bl = zeros(size(ac_mag));
                ac_grad_bl = zeros(size(ac_grad));
                % first magnetometer
                % filter
                cfg            = [];
                cfg.toilim    = [t(1)/srate t(end)/srate];
                chunk_mag = ft_redefinetrial(cfg,data_mag);
                chunk_grad = ft_redefinetrial(cfg,data_grad);
                chunk_mag_bl = ft_redefinetrial(cfg,data_mag_bl);
                chunk_grad_bl = ft_redefinetrial(cfg,data_grad_bl);
                
                % prewhiten
                cfg = [];
                cfg.derivative = 'yes';
                chunk_mag = ft_preprocessing(cfg, chunk_mag);
                chunk_grad = ft_preprocessing(cfg,chunk_grad);
                chunk_mag_bl = ft_preprocessing(cfg,chunk_mag_bl);
                chunk_grad_bl = ft_preprocessing(cfg,chunk_grad_bl);
                
                %                         % get power spectra via pwelch method
                %                         spectra = zeros(nTrial,numel(f_range), nNode);
                %                         for m = 1:nTrial
                %                             [spectra(m,:,:), fr] = pwelch(chunk_mag.trial{m}', [], [] ,f_range, srate, 'power');
                %                         end
                %                         spectra = squeeze(mean(spectra));
                %
                %  % icov
                %                         cnt = 1;
                %                         S = cov(h');
                %                         lambda = logspace(-3,2,20);
                %                         sparsity = zeros(1,numel(lambda));
                %                         for l = lambda
                %                             A = L1precisionBCD(S,l);
                %                             sparsity(cnt) = sum(sum(abs(A == 0)))/numel(A);
                %                             cnt = cnt + 1;
                %                         end
                
                % get amplitude envelope
                % bandpass filter and get analytical signal
                % magnitude
                for f = 1:size(freqs,1)
                    f_range = freqs(f,:);
                    amp_corr = NetBCI_get_amp_corr(chunk_mag, t, srate, freqs(f,:));
                    
                    % now for baseline
                    amp_corr_bl = NetBCI_get_amp_corr(chunk_mag_bl, t, srate, freqs(f,:));
                    
                    
                    ac_mag(:,f) = reshape(amp_corr(logical(triu(ones(size(amp_corr)),1))), 1, []);
                    ac_mag_bl(:,f) = reshape(amp_corr_bl(logical(triu(ones(size(amp_corr_bl)),1))), 1, []);
                    
                    %plot
                    figure(1); clf
                    imagesc(amp_corr);
                    colorbar; pause(0.01)
                    figure(2); clf
                    imagesc(amp_corr_bl);
                    colorbar; pause(0.01)
                    
                    
                    % now gradiometer
                    
                    % pick 1 gradiometer
                    cfg = [];
                    cfg.channel = data.label(strcmp(data_raw.meg_sensors.chanunit, 'T/cm'));
                    cfg.channel = cfg.channel(1:2:numel(cfg.channel));
                    chunk_grad = ft_preprocessing(cfg,chunk_grad);
                    chunk_grad_bl = ft_preprocessing(cfg,chunk_grad_bl);
                    
                    
                    % get amplitude envelope
                    % bandpass filter and get analytical signal
                    % magnitude
                    amp_corr = NetBCI_get_amp_corr(chunk_grad, t, srate, freqs(f,:));
                    
                    % now for baseline
                    amp_corr_bl = NetBCI_get_amp_corr(chunk_grad_bl, t, srate, freqs(f,:));
                    
                    
                    ac_grad(:,f) = reshape(amp_corr(logical(triu(ones(size(amp_corr)),1))), 1, []);
                    ac_grad_bl(:,f) = reshape(amp_corr_bl(logical(triu(ones(size(amp_corr_bl)),1))), 1, []);
                    
                    %plot
                    figure(1); clf
                    imagesc(amp_corr);
                    colorbar ; pause(0.01)
                    figure(2); clf
                    imagesc(amp_corr_bl);
                    colorbar; pause(0.01)
                    
                    
                    
                end
                figure(3); clf
                imagesc(ac_grad); colorbar;
                saveas(gca, [img_dir, 'grad/_all.png'], 'png');
                figure(4); clf
                imagesc(ac_mag); colorbar;
                saveas(gca, [img_dir, 'mag/_all.png'], 'png');
                
                
                % reshape to vector
                ac_grad = reshape(ac_grad, [], 1);
                ac_grad_bl = reshape(ac_grad_bl, [], 1);
                ac_mag = reshape(ac_mag, [], 1);
                ac_mag_bl = reshape(ac_mag_bl, [], 1);
                
                num_neg_grad = [num_neg_grad, sum(ac_grad < 0)/numel(ac_grad)];
                num_neg_mag = [num_neg_mag, sum(ac_mag < 0)/numel(ac_mag)];
                ac_grad(ac_grad < 0) = 0;
                ac_mag(ac_mag < 0) = 0;
                ac_grad_bl(ac_grad_bl < 0) = 0;
                ac_mag_bl(ac_mag_bl < 0) = 0;
                save([save_dir, 'NMF_all_grad_ac'], 'ac_grad');
                save([save_dir, 'NMF_all_mag_ac'], 'ac_mag');
                save([save_dir, 'NMF_all_grad_bl_ac'], 'ac_grad_bl');
                save([save_dir, 'NMF_all_mag_bl_ac'], 'ac_mag_bl');
                
            catch
                errors{cnte} = [d.folder, '/', d.name];
                cnte = cnte + 1;
            end
        end
    end
end







