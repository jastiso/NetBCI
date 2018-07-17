%% Test the variability of null models
% to save time, you can sometimes use one null model per trial, rather than
% a dist of 1000 models for every trial. however, you need to be sure that
% you variability across instantiations is small, compared to your
% variability across trials. This script tests that.

% @author JStiso jeni.stiso@gmail.com
% Change Log:
% July 9, 2018 - created script
% July 12 - changed baseline to 0.2-0.8

%% Define Global Variables

addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830')
addpath(genpath('/Users/stiso/Documents/MATLAB/eeglab_current'))

top_dir = '/Users/stiso/Documents/MATLAB/NetBCI/';
data_dir = '/Users/stiso/Resilio Sync/NETBCI.RAW/DataBase/1_Signals/2_Segmentation/2_MEG/';
raw_dir = '/Users/stiso/Resilio Sync/NETBCI.RAW/DataBase/1_Signals/0_RawData/2_MEG/';
sessions = [{'Session1'}, {'Session2'}, {'Session3'}, {'Session4'}];
condition = [{'test01'}, {'test02'}, {'test03'}, {'test04'}, {'test05'}, {'test06'}]; % not including rest
nNode = 102;
nChunks = 1; %number of sets to divide trials into
freqs = [7,14;15,30;31,45;55,70];
bands = [{'alpha'},{'beta'},{'low_gamma'}, {'gamma'}];
st = 3;
en = 6; % in seconds, the feedback period: 3-6s
bl_st = 0.2; % half o fthe intertrial interval
bl_en = 0.8;
num_neg_grad = [];
num_neg_mag = [];
cnte = 1;
% GC params
dt = 100;
lag = 10;
t0 = dt+lag;


%% Get nulls networks
% initialize final struct
nulls = struct();
nulls.trial_dist = cell(1,numel(bands));
nulls.null_dist = cell(1,numel(bands));

for i = 1:numel(sessions)
    session = sessions{i};
    for j = 1:numel(condition)
        cond = condition{j};
        subj = 16;
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
            data_grad = data;
            % grad
            data_grad.label = data_grad.label(grad_idx);
            for m = 1:numel(data_grad.trial)
                data_grad.trial{m} = data_grad.trial{m}(grad_idx,:);
            end
            % phase randomize
            for t = 1:nTrial
                data_grad.trial{t} = linsurr(data_grad.trial{t});
            end
            
            
            
            
            % Make connectivity matrices
            
            % get power in task and BL correct
            % use multitapers for gamma, and single tapers for low
            % freq
            
            %gradiometers
            [freq_grad_mt, freq_grad_st] = NetBCI_get_pow(data_grad, st, en, bl_st, bl_en);
            
            
            % plot
            figure(1); clf
            imagesc(squeeze(freq_grad_mt.powspctrm(end,94,:,:)))
            colorbar
            pause(0.01)
            figure(3); clf
            imagesc(squeeze(freq_grad_st.powspctrm(end,94,:,:)))
            colorbar
            pause(0.01)
            
            
            for f = 1:numel(bands)
                f_range = freqs(f,:);
                curr_band = bands{f};
                % average to get current band
                if contains(curr_band, 'gamma')
                    
                    % gradiometers
                    % initialize
                    freq_grad = freq_grad_mt;
                    % get freq range
                    idx = freq_grad_mt.freq >= f_range(1) & freq_grad_mt.freq < f_range(end);
                    % average
                    pow = (freq_grad_mt.powspctrm(:,:,idx,:));
                    pow = squeeze(nanmean(pow,3));
                    freq_grad.powspctrm = pow;
                    freq_grad.freq = squeeze(mean(freq_grad.freq(idx)));
                    
                else
                    % gradiometers
                    % initialize
                    freq_grad = freq_grad_st;
                    % get freq range
                    idx = freq_grad_st.freq >= f_range(1) & freq_grad_st.freq < f_range(end);
                    % average
                    pow = (freq_grad_st.powspctrm(:,:,idx,:));
                    pow = squeeze(nanmean(pow,3));
                    freq_grad.powspctrm = pow;
                    freq_grad.freq = squeeze(mean(freq_grad.freq(idx)));
                    
                end
                
                % combine gradiometers
                freq_grad = combine_gradiometers(freq_grad);
                
                % plot
                figure(1); clf
                imagesc(squeeze(freq_grad.powspctrm(end,:,:))); colorbar
                
                
                % get cov based GC
                gc_grad = zeros((nNode^2-nNode)/2, nTrial);
                for t = 1:nTrial
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
                    
                    % save
                    gc_grad(:,t) = GC;
                end
                nulls.trial_dist{f} = [nulls.trial_dist{f}, gc_grad];
            end
        catch
            errors{cnte,1} = [d.folder, '/', d.name];
            cnte = cnte + 1;
        end
    end
end
%% now get dist across nulls

nNull = size(nulls.trial_dist{1},2);
% get distribution of nulls
null_inst = cell(1,nNull);
for m = 1:nNull
    null_inst{m} = data_grad;
    for t = 1:nTrial
        null_inst{m}.trial{t} = linsurr(data_grad.trial{t});
    end
end




% Make connectivity matrices

% get power in task and BL correct
% use multitapers for gamma, and single tapers for low
% freq

for m = 1:nNull
    curr_data = null_inst{m};
    %gradiometers
    [freq_grad_mt, freq_grad_st] = NetBCI_get_pow(curr_data, st, en, bl_st, bl_en);
    
    
    % plot
    figure(1); clf
    imagesc(squeeze(freq_grad_mt.powspctrm(end,94,:,:)))
    colorbar
    pause(0.01)
    figure(3); clf
    imagesc(squeeze(freq_grad_st.powspctrm(end,94,:,:)))
    colorbar
    pause(0.01)
    
    
    for f = 1:numel(bands)
        f_range = freqs(f,:);
        curr_band = bands{f};
        % average to get current band
        if contains(curr_band, 'gamma')
            
            % gradiometers
            % initialize
            freq_grad = freq_grad_mt;
            % get freq range
            idx = freq_grad_mt.freq >= f_range(1) & freq_grad_mt.freq < f_range(end);
            % average
            pow = (freq_grad_mt.powspctrm(:,:,idx,:));
            pow = squeeze(nanmean(pow,3));
            freq_grad.powspctrm = pow;
            freq_grad.freq = squeeze(mean(freq_grad.freq(idx)));
            
        else
            % gradiometers
            % initialize
            freq_grad = freq_grad_st;
            % get freq range
            idx = freq_grad_st.freq >= f_range(1) & freq_grad_st.freq < f_range(end);
            % average
            pow = (freq_grad_st.powspctrm(:,:,idx,:));
            pow = squeeze(nanmean(pow,3));
            freq_grad.powspctrm = pow;
            freq_grad.freq = squeeze(mean(freq_grad.freq(idx)));
            
        end
        
        % combine gradiometers
        freq_grad = combine_gradiometers(freq_grad);
        
        % plot
        figure(1); clf
        imagesc(squeeze(freq_grad.powspctrm(end,:,:))); colorbar
        
        
       
        %gradiometers
        tmp = zeros(nNode);
        curr = squeeze(freq_grad.powspctrm(1,:,:));
        [GC, pairs] = cov_GC(curr,dt,lag,t0);
        GC = sum(GC,2);
        % plot
        tmp(logical(tril(ones(nNode),-1))) = GC;
        tmp = tmp + tmp';
        figure(5); clf
        imagesc(tmp); colorbar; pause(0.001)
        
        % save
        nulls.null_dist{f} = [nulls.null_dist{f}, GC];
    end
end


save([top_dir, 'nulls.mat'], 'nulls', '-v7.3');

%% Plot

load([top_dir, 'nulls.mat']);
img_dir = [top_dir, 'null_dist/'];

% for every band, plot histograms of summary statistics
for i = 1:numel(bands)
   f = bands{i};
   trial_dist = nulls.trial_dist{i};
   null_dist = nulls.null_dist{i};
   
   % mean
   [N, X] = hist(mean(null_dist,1),50);
   [N2, X2] = hist(mean(trial_dist,1),50);
   figure(1); clf
   bar(X, N, 0.8, 'facecolor', [230,27,68]./255, 'edgecolor', [230,27,68]./255,'facealpha', .5); hold on
   bar(X2, N2, 0.8, 'facecolor', [0,121,194]./255, 'edgecolor', [0,121,194]./255,'facealpha', .5);
   legend([{'null'}, {'trial'}])
   hold off
   saveas(gca, [img_dir, 'mean_', f, '.png'])
   
   % all values
   [N, X] = hist(reshape(null_dist,1,[]),50);
   [N2, X2] = hist(reshape(trial_dist,1,[]),50);
   figure(2); clf
   bar(X, N, 0.8, 'facecolor', [230,27,68]./255, 'edgecolor', [230,27,68]./255,'facealpha', .5); hold on
   bar(X2, N2, 0.8, 'facecolor', [0,121,194]./255, 'edgecolor', [0,121,194]./255,'facealpha', .5);
   legend([{'null'}, {'trial'}])
   hold off
   saveas(gca, [img_dir, 'all_', f, '.png'])
   
   % var
   [N, X] = hist(std(null_dist,0,1),50);
   [N2, X2] = hist(std(trial_dist,0,1),50);
   figure(3); clf
   bar(X, N, 0.8, 'facecolor', [230,27,68]./255, 'edgecolor', [230,27,68]./255,'facealpha', .5); hold on
   bar(X2, N2, 0.8, 'facecolor', [0,121,194]./255, 'edgecolor', [0,121,194]./255,'facealpha', .5);
   legend([{'null'}, {'trial'}])
   hold off
   saveas(gca, [img_dir, 'sd_', f, '.png'])
   
   
end
