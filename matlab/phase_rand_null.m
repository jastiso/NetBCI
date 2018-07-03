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
nSubj = 20;
nChunks = 1; %number of sets to divide trials into
freqs = [3,6;7,14;15,30;31,45;55,70];
bands = [{'theta'},{'alpha'},{'beta'},{'low_gamma'}, {'gamma'}];
dur = 3.5; % in seconds, the length of time you want going into each window
num_neg_grad = [];
num_neg_mag = [];
cnte = 0;

%% Start Loop

for i = 1:numel(sessions)
    session = sessions{i};
    for j = 1:numel(condition)
        cond = condition{j};
        for k = 1:nSubj
            ext = sprintf('%03d',k);
            % load data
            d = dir([data_dir, session, '/*', cond, '*', '/Seg_MEG_Subj', ext, '_Ses', session(end), '_', cond, '.mat']);
            d_raw = dir([raw_dir, session, '/*', cond, '*', '/RawData_MEG_Subj', ext, '_Ses', session(end), '_', cond, '.mat']);
            try
                load([d.folder, '/', d.name]);
                img_dir = [top_dir, session, '/', cond, '/', ext, '/diagnostics/pr/'];
                save_dir = [top_dir, session, '/', cond, '/', ext, '/FCmatrices/pr/'];
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
                
                % create phase randomized surrogate
                for t = 1:nTrial
                    data.trial{t} = linsurr(data.trial);
                    data_bl.trial{t} = linsurr(data_bl.trial);
                end
                
                % check fft
                figure(1); clf
                spectopo(data.trial{1}, 0, round(data.fsample));
                pause(.1)
                saveas(gca, [img_dir, 'fft_before.png'], 'png')
                
                % filter
                data = NetBCI_filtering(data);
                data_bl = NetBCI_filtering(data_bl);
                
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
                
                
                %% Make connectivity matrices
                for f = 1:numel(bands)
                    cnt = 0;
                    conn_tv_mag = zeros((size(data_mag.trial{1},1)^2-size(data_mag.trial{1},1))/2, numel(1:dur*srate:trl_len));
                    conn_tv_grad = zeros((size(data_mag.trial{1},1)^2-size(data_mag.trial{1},1))/2, numel(1:dur*srate:trl_len)); % after combining, grad and nag will be the same size
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
                            
                            % get power spectra via pwelch method
                            %[spectra, fr] = pwelch(chunk_mag.trial{1}', [], [] ,freqs(f,:), srate, 'power')
                            
                            % get wPLI using fieldtrip
                            if n ~= nChunks
                                trials = st:floor(st+nTrial/nChunks)-1;
                            else
                                trials = st:nTrial;
                            end
                            wpli = NetBCI_get_wPLI(chunk_mag, trials, 'mag');
                            wpli_bl = NetBCI_get_wPLI(chunk_mag_bl, trials, 'mag');
                            % reshape and save
                            conn_tv_mag(:,cnt) = reshape(wpli.wpli_debiasedspctrm(logical(triu(ones(size(wpli.wpli_debiasedspctrm)),1))), 1, []);
                            save([save_dir, 'mag/', bands{f}, '_wpli'], 'wpli')
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
                            wpli = NetBCI_get_wPLI(chunk_grad, trials, 'grad');
                            wpli_bl = NetBCI_get_wPLI(chunk_grad_bl, trials, 'grad');
                            % reshape and save
                            conn_tv_mag(:,cnt) = reshape(wpli.wpli_debiasedspctrm(logical(triu(ones(size(wpli.wpli_debiasedspctrm)),1))), 1, []);
                            save([save_dir, 'mag/', bands{f}, '_wpli'], 'wpli')
                            conn_tv_mag_bl(:,cnt) = reshape(wpli_bl.wpli_debiasedspctrm(logical(triu(ones(size(wpli_bl.wpli_debiasedspctrm)),1))), 1, []);
                            save([save_dir, 'mag/', bands{f}, '_wpli_bl'], 'wpli_bl')
                            
                            
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
                    save([save_dir, 'NMF_', bands{f}, '_grad'], 'conn_tv_grad');
                    save([save_dir, 'NMF_', bands{f}, '_mag'], 'conn_tv_mag');
                    
                end
            catch
                errors{cnte} = [data_dir, session, '/*', cond, '*', '/Seg_MEG_Subj', ext, '_Ses', session(end), '_', cond, '.mat'];
                cnte = cnte + 1;
            end
        end
    end
end





