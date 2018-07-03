%% Analyze resting state data
%% Define Global Variables

addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830')
addpath(genpath('/Users/stiso/Documents/MATLAB/eeglab_current'))

top_dir = '/Users/stiso/Documents/MATLAB/NetBCI/';
data_dir = [top_dir, 'DataBase/1_Signals/1_PreProcessing/2_MEG/'];
raw_dir = [top_dir, 'DataBase/1_Signals/0_RawData/2_MEG/'];
sessions = [{'Session1'}, {'Session2'}, {'Session3'}, {'Session4'}];
condition = [{'RS1_BEGIN'}, {'RS2_END'}]; % not including rest
nSubj = 20;
nChunks = 2; %number of sets to divide trials into
freqs = [3,6;7,14;15,30;40,69];
bands = [{'theta'},{'alpha'},{'beta'},{'gamma'}];
dur = 2; % in seconds, the length of time you want going into each window
num_neg_grad = [];
num_neg_mag = [];

%% Start Loop

for i = 1:numel(sessions)
    session = sessions{i};
    for j = 1:numel(condition)
        cond = condition{j};
        c_ext = strsplit(cond, '_');
        c_ext = c_ext{1};
        for k = 1:nSubj
            ext = sprintf('%03d',k);
            % load data
            d = dir([data_dir, session, '/*', cond, '*', '/PreProc_MEG_Subj', ext, '_Ses', num2str(i), '_', c_ext, '.mat']);
            d_raw = dir([raw_dir, session, '/*', cond, '*', '/RawData_MEG_Subj', ext, '_Ses', num2str(i), '_', c_ext, '.mat']);
            if numel(d) > 0
                load([d.folder, '/', d.name]);
                img_dir = [top_dir, session, '/', cond, '/', ext, '/diagnostics/'];
                save_dir = [top_dir, session, '/', cond, '/', ext, '/FCmatrices/'];
                load([d_raw.folder, '/', d_raw.name]);
                try
                    eval(['data = ', d.name(1:end-4), ';']);
                    eval(['clear ', d.name(1:end-4)]);
                    eval(['data_raw = ', d_raw.name(1:end-4), ';']);
                    eval(['clear ', d_raw.name(1:end-4)]);
                catch
                    eval(['data = ', d.name(1:13), '16', d.name(16:(end-4)), ';']);
                    eval(['clear ', d.name(1:13), '16', d.name(16:(end-4))]);
                    eval(['data_raw = ', d_raw.name(1:17), '16', d_raw.name(20:(end-4)), ';']);
                    eval(['clear ', d_raw.name(1:17), '16', d_raw.name(20:(end-4))]);
                end
                srate = data.fsample;
                rs_len = numel(data.trial{1,1});
                
                % make directories
                if ~exist([img_dir, '/grad/RS/'], 'dir')
                    mkdir([img_dir, '/grad/RS/'])
                end
                % make directories
                if ~exist([img_dir, '/mag/RS/'], 'dir')
                    mkdir([img_dir, '/mag/RS/'])
                end
                % make directories
                if ~exist([save_dir, '/grad/RS/'], 'dir')
                    mkdir([save_dir, '/grad/RS/'])
                end
                % make directories
                if ~exist([save_dir, '/mag/RS/'], 'dir')
                    mkdir([save_dir, '/mag/RS/'])
                end
                
                 % filter
                % now filter: 50 and 100
                % butter(order, [frequencies], type); stop means sharp cutoff
                [a,b] = butter(4, [48.5/(data.fsample/2), 51.5/(data.fsample/2)], 'stop');
                data.trial{1} = filtfilt(a,b,data.trial{1}')'; % important to use filtfilt and not filt here, to preserve phase relationship
                
                % high pass filter
                [a,b] = butter(4, 1/(data.fsample/2), 'high');
                data.trial{1} = filtfilt(a,b,data.trial{1}')';
                
                
                % low pass filter
                [a,b] = butter(4, 70/(data.fsample/2), 'low');
                data.trial{1} = filtfilt(a,b,data.trial{1}')';
                
                
                % image again to see what changed
                figure(1); clf
                spectopo(data.trial{1}, 0, floor(data.fsample));
                pause(.01)
                saveas(gca, [img_dir, 'fft_after.png'], 'png')
                
                % separate into magnetometers and gradiometers
                mag_idx = strcmp(data_raw.meg_sensors.chantype,'megmag');
                grad_idx = ~mag_idx;
                data_mag = data;
                data_grad = data;
                data_grad.label = data_grad.label(grad_idx);
                for m = 1:numel(data_grad.trial)
                    data_grad.trial{m} = data_grad.trial{m}(grad_idx,:);
                end
                data_mag.label = data_mag.label(mag_idx);
                for m = 1:numel(data_mag.trial)
                    data_mag.trial{m} = data_mag.trial{m}(mag_idx,:);
                end
                

                %% Make connectivity matrices
                for f = 1:numel(bands)
                    cnt = 0;
                    conn_tv_mag = zeros((size(data_mag.trial{1},1)^2-size(data_mag.trial{1},1))/2, numel(1:dur*srate:rs_len));
                    conn_tv_grad = zeros((size(data_grad.trial{1},1)^2-size(data_grad.trial{1},1))/2, numel(1:dur*srate:rs_len));

                    st = 1;
                        %for t  = 1:dur*srate:rs_len
                            cnt = cnt+1;
                            % first magnetometer
                            % filter
%                             cfg            = [];
%                             cfg.toilim    = [t/srate (t+srate-1)/srate];
%                             chunk_mag = ft_redefinetrial(cfg,data_mag);
%                             chunk_grad = ft_redefinetrial(cfg,data_grad);
                            
                            
                            cfg.output     = 'pow';
                            cfg.method     = 'mtmconvol';
                            cfg.foi         = freqs(f,1):freqs(f,2);
                            cfg.t_ftimwin  = 5./cfg.foi;
                            cfg.tapsmofrq  = 0.4 *cfg.foi;
                            cfg.toi        = 0:0.05:30;
                            freq    = ft_freqanalysis(cfg, data_mag);
                            freq.freq = mean(freq.freq);
                            freq.powspctrm = squeeze(nanmean(freq.powspctrm,2));
                            freq.powspctrm(isnan(freq.powspctrm)) = [];
%                             cfg = [];
%                             cfg.method = 'wpli_debiased';
%                             cfg.keeptrials = 'yes';
%                             %cfg.complex = 'wpli';
%                             wpli = ft_connectivityanalysis(cfg, freq);
%                             conn_tv_mag(:,cnt) = reshape(wpli.wpli_debiasedspctrm(logical(triu(ones(size(wpli.wpli_debiasedspctrm)),1))), 1, []);
%                             save([save_dir, 'mag/', bands{f}, '_wpli'], 'wpli')
                            
                            %plot
                            figure(1); clf
                            imagesc(wpli.wpli_debiasedspctrm);
                            colorbar
                            saveas(gca, [img_dir, 'mag/', bands{f}, '_wpli.png'], 'png')
                            
                            
                            
                            % now gradiometer
                            % filter
                            cfg            = [];
                            cfg.output     = 'fourier';
                            cfg.method     = 'mtmfft';
                            cfg.foilim     = freqs(f,:);
                            cfg.tapsmofrq  = 5;
                            freq    = ft_freqanalysis(cfg, chunk_grad);
                            freq.freq = mean(freq.freq);
                            freq.fourierspctrm = mean(freq.fourierspctrm,3);
                            
                            cfg = [];
                            cfg.method = 'wpli_debiased';
                            %cfg.complex = 'wpli';
                            wpli = ft_connectivityanalysis(cfg, freq);
                            conn_tv_grad(:,cnt) = reshape(wpli.wpli_debiasedspctrm(logical(triu(ones(size(wpli.wpli_debiasedspctrm)),1))), 1, []);
                            save([save_dir, 'grad/', bands{f}, '_wpli'], 'wpli')

                            
                            %plot
                            figure(3); clf
                            imagesc(wpli.wpli_debiasedspctrm);
                            colorbar
                            saveas(gca, [img_dir, 'grad/', bands{f}, '_wpli.png'], 'png')
                        %end
                        st = st + floor(st+nTrial/nChunks)-1;
                    num_neg_grad = [num_neg_grad, sum(conn_tv_grad < 0)/numel(conn_tv_grad)];
                    num_neg_mag = [num_neg_mag, sum(conn_tv_mag < 0)/numel(conn_tv_mag)];
                    conn_tv_grad(conn_tv_grad < 0) = 0;
                    conn_tv_mag(conn_tv_mag < 0) = 0;

                    save([save_dir, 'NMF_', bands{f}, '_grad'], 'conn_tv_grad');
                    save([save_dir, 'NMF_', bands{f}, '_mag'], 'conn_tv_mag');

                end
            end
        end
    end
end