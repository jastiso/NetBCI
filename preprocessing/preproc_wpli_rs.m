%% Resting State MEG data
% get connectivity matrices for pre task MEG
% @author JStiso March 2020


%% Define Global Variables
clear
clc

addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830')
addpath(genpath('/Users/stiso/Documents/MATLAB/eeglab_current'))
addpath(genpath('/Users/stiso/Documents/Code/netBCI/matlab/'))

top_dir = '/Users/stiso/Documents/MATLAB/NetBCI/';
data_dir = '/Users/stiso/Resilio Sync/NETBCI.RAW/DataBase/1_Signals/1_PreProcessing/2_MEG/';
raw_dir = '/Users/stiso/Resilio Sync/NETBCI.RAW/DataBase/1_Signals/0_RawData/2_MEG/';
cnte = 1;

sessions = [{'Session4'}];
condition = [{'RS2_END'}]; % not including rest
subjs = [1:20];
nNode = 102;

freqs = [7,14;15,30;31,45];
bands = [{'alpha'},{'beta'},{'low_gamma'}];

for s = subjs
    subj = s;
    for i = 1:numel(sessions)
        session = sessions{i};
        for j = 1:numel(condition)
            cond = condition{j};
            ext = sprintf('%03d',subj);
            % load data
            d = dir([data_dir, session, '/*', cond, '*', '/PreProc_MEG_Subj', ext, '_Ses', session(end), '_RS2.mat']);
            d_raw = dir([raw_dir, session, '/*', cond, '*', '/RawData_MEG_Subj', ext, '_Ses', session(end), '_RS2.mat']);
            %try
                load([d.folder, '/', d.name]);
                img_dir = [top_dir, session, '/', cond, '/', ext, '/diagnostics/'];
                save_dir = [top_dir, session, '/', cond, '/', ext, '/FCmatrices/'];
                if ~exist(img_dir, 'dir')
                    mkdir(img_dir);
                end
                if ~exist(save_dir, 'dir')
                    mkdir(save_dir);
                end
                
                try
                    load([d_raw.folder, '/', d_raw.name]);
                catch
                    d_raw = dir([raw_dir, 'Session1/*', cond, '*', '/RawData_MEG_Subj', ext, '_Ses4_RS2.mat']);
                    load([d_raw.folder, '/', d_raw.name]);
                end
                try
                    eval(['data = ', d.name(1:end-4), ';']);
                    eval(['clear ', d.name(1:end-4), ';']);
                    
                catch
                    eval(['data = ', d.name(1:17), '16', d.name(20:(end-4)), ';']);
                    eval(['clear ', d.name(1:17), '16', d.name(20:(end-4))]);
                end
                try
                    eval(['data_raw = ', d_raw.name(1:end-4), ';']);
                    eval(['clear ', d_raw.name(1:end-4)]);
                catch
                    eval(['data_raw = ', d_raw.name(1:17), '16', d_raw.name(20:(end-4)), ';']);
                    eval(['clear ', d_raw.name(1:17), '16', d_raw.name(20:(end-4))]);
                end
                
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
                dur = size(data.trial{1},2)/srate;
                
                % separate into magnetometers and gradiometers
                data.senstype = 'neuromag306';
                mag_idx = strcmp(data_raw.meg_sensors.chantype,'megmag');
                grad_idx = ~mag_idx;
                data_grad = data;
                % grad
                data_grad.label = data_grad.label(grad_idx);
                for m = 1:numel(data_grad.trial)
                    data_grad.trial{m} = data_grad.trial{m}(grad_idx,:);
                end
                
                
                %% Make connectivity matrices
                
                % combine gradiometers
                cfg = [];
                cfg.method = 'sum';
                data_grad = ft_combineplanar(cfg,data_grad);
                
                for f = 1:numel(bands)
                    f_range = freqs(f,1):freqs(f,2);
                    curr_band = bands{f};
                    
                    % get wPLI
                    %gradiometers
                    cfg = [];
                    cfg.method     = 'wavelet';
                    cfg.width      = 7;
                    cfg.output     = 'powandcsd';
                    cfg.foi        = f_range;
                    cfg.toi        = 0:0.5:dur;
                    cfg.pad        = 'nextpow2';
                    wave = ft_freqanalysis(cfg, data_grad);
                    
                    % housekeeping for formatting things later (basically to move between matrix and vector representations). There might be a better way to
                    % do this
                    nEdge = (nNode^2 - nNode)/2;
                    cnt = 0;
                    label = zeros(nEdge,2);
                    for n = 1:nNode
                        for m = (n+1):nNode
                            cnt = cnt +1;
                            label(cnt,:) = [n,m];
                        end
                    end
                    
                    X = imag(wave.crsspctrm);
                    wpli_v = abs(nanmean(X,3))./nanmean(abs(X),3);
                    
                    %average over freqs
                    wpli_v = abs(mean(wpli_v,2));
                    
                    % put back into matrix
                    wpli = zeros(nNode);
                    for n = 1:nNode
                        for m = (n+1):nNode
                            row_idx = find(label(:,1)==n);
                            col_idx = find(label(:,2)==m);
                            idx = intersect(row_idx,col_idx);
                            wpli(n,m) = wpli_v(idx);
                        end
                    end
                    wpli = wpli+wpli';
                    
                    % plot
                    figure(1); clf
                    imagesc(wpli); colorbar;
                    
                    % save
                    if any(any(isnan(wpli)))
                        error('This matrix contains nans')
                    end
                    save([save_dir, curr_band, '_rs_wpli.mat'], 'wpli')
                end
%             catch
%                 errors{cnte,1} = [d.folder, '/', d.name];
%                 cnte = cnte + 1;
%             end
        end
    end
end