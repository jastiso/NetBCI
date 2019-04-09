%% Get "brain states" for network control validation
% Get power centered at each timepoint across sensors

% Change Log
% April 5, 2019 - created

%% Define Global Variables

addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830')
addpath(genpath('/Users/stiso/Documents/MATLAB/eeglab_current'))
addpath(genpath('/Users/stiso/Documents/Code/NetBCI/matlab/'))

top_dir = '/Users/stiso/Documents/MATLAB/NetBCI/';
data_dir = '/Users/stiso/Resilio Sync/NETBCI.RAW/DataBase/1_Signals/2_Segmentation/2_MEG/';
raw_dir = '/Users/stiso/Resilio Sync/NETBCI.RAW/DataBase/1_Signals/0_RawData/2_MEG/';
sessions = [{'Session1'}, {'Session2'}, {'Session3'}, {'Session4'}];
condition = [{'test01'}, {'test02'}, {'test03'}, {'test04'}, {'test05'}, {'test06'}]; % not including rest
subjs = [1:20];
nNode = 102;

freq_range = [7,14;15,30;31,45];
bands = [{'alpha'},{'beta'},{'low_gamma'}];

st = 3;
en = 6; % in seconds, the feedback period: 3-6s

%% Loop through subjects

for k = subjs
    ext = sprintf('%03d',k);
    states = cell(numel(bands),1);
    for i = 1:numel(sessions)
        session = sessions{i};
        for j = 1:numel(condition)
            cond = condition{j};
            
            % load data
            d = dir([data_dir, session, '/*', cond, '*', '/Seg_MEG_Subj', ext, '_Ses', session(end), '_', cond, '.mat']);
            d_raw = dir([raw_dir, session, '/*', cond, '*', '/RawData_MEG_Subj', ext, '_Ses', session(end), '_', cond, '.mat']);
            load([d.folder, '/', d.name]);
            img_dir = [top_dir, session, '/', cond, '/', ext, '/diagnostics/'];
            save_dir = [top_dir, session, '/', cond, '/', ext, '/states/'];
            try
                load([d_raw.folder, '/', d_raw.name]);
            catch
                d_raw = dir([raw_dir, 'Session1/*', cond, '*', '/RawData_MEG_Subj', ext, '_Ses1_', cond, '.mat']);
                load([d_raw.folder, '/', d_raw.name]);
            end
            try
                eval(['data = ', d.name(1:end-4), '.MI;']);
                eval(['clear ', d.name(1:end-4)]);
                
            catch
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
            if ~exist([img_dir], 'dir')
                mkdir([img_dir])
            end
            % make directories
            if ~exist([save_dir], 'dir')
                mkdir([save_dir])
            end
            
            % view the data (optional)
            cfg = [];
            cfg.viewmode = 'vertical';
            %ft_databrowser(cfg, data)
            
            
            % check fft - will fix later if there is noise
            figure(1); clf
            spectopo(data.trial{1}, 0, round(data.fsample));
            pause(.1)
            saveas(gca, [img_dir, 'fft_before.png'], 'png')
            
            % cut out trials that have NaN values at the end...namely subj
            % 13
            for m = 1:numel(data.trial)
                if any(any(isnan(data.trial{m})))
                    data.trial = data.trial([1:(m-1),(m+1):end]);
                    data.time = data.time([1:(m-1),(m+1):end]);
                    data.sampleinfo = data.sampleinfo([1:(m-1),(m+1):end],:);
                end
            end
            
            % grab some more variables
            trl_len = numel(data.trial{1}(1,:));
            srate = data.fsample;
            nTrial = numel(data.trial);
            
            
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
            
            
            %% Make state matrices
            
            % combine gradiometers
            cfg = [];
            cfg.method = 'sum';
            data_grad = ft_combineplanar(cfg,data_grad);
            
            
            for f = 1:numel(bands)
                f_range = freq_range(f,1):freq_range(f,2);
                curr_band = bands{f};
                
                %gradiometers
                for t = 1:nTrial
                    pow = get_wavelet_power(data_grad, srate, st, en, f_range, t);
                    
                    % concat
                    states{f} = [states{f}, pow];
                end
            end
        end
    end
    
    for f = 1:numel(bands)
        states{f} = zscore(states{f},0,2);
        figure(f+1); clf
        imagesc(states{f}); colorbar;
        caxis([-3,5]); title(bands{f});
        saveas(gca, [img_dir, bands{f}, '_wpli_pow.png'], 'png')
    end
    save([save_dir, 'states.mat'], 'states')
end
