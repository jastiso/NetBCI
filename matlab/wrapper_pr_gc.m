function [ errors ] = wrapper_pr_gc( session, condition, subjs, data_dir, raw_dir, top_dir, bands, freqs, bl_st, bl_en, dt, lag, t0, st, en, nNode )
% helper function for paralelizing pipeline
errors = {};
cnte = 1;

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
            
            % phase randomize
            for t = 1:nTrial
                data_mag.trial{t} = linsurr(data_mag.trial{t});
                data_grad.trial{t} = linsurr(data_grad.trial{t});
            end
            
            
            
            %% Make connectivity matrices
            
            % get power in task and BL
            % use multitapers for gamma, and single tapers for low
            % freq
            
            % magnetometers
            [freq_mag_mt, freq_mag_st] = NetBCI_get_pow(data_mag, st, en, bl_st, bl_en);
            %gradiometers
            [freq_grad_mt, freq_grad_st] = NetBCI_get_pow(data_grad, st, en, bl_st, bl_en);
            
            
            % plot
            figure(1); clf
            imagesc(squeeze(freq_mag_mt.powspctrm(end,94,:,:)))
            colorbar
            pause(0.01)
            figure(3); clf
            imagesc(squeeze(freq_mag_st.powspctrm(end,94,:,:)))
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
                    % get freq range
                    idx = freq_mag_mt.freq >= f_range(1) & freq_mag_mt.freq < f_range(end);
                    % average
                    pow = (freq_mag_mt.powspctrm(:,:,idx,:));
                    pow = squeeze(nanmean(pow,3));
                    freq_mag.powspctrm = pow;
                    freq_mag.freq = squeeze(mean(freq_mag.freq(idx)));
                    
                    
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
                    % magnetometers
                    % initialize
                    freq_mag = freq_mag_st;
                    % get freq range
                    idx = freq_mag_st.freq >= f_range(1) & freq_mag_st.freq < f_range(end);
                    % average
                    pow = (freq_mag_st.powspctrm(:,:,idx,:));
                    pow = squeeze(nanmean(pow,3));
                    freq_mag.powspctrm = pow;
                    freq_mag.freq = squeeze(mean(freq_mag.freq(idx)));
                    
                    
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
                saveas(gca, [img_dir, bands{f}, 'pow_pr.png'], 'png')
                
                
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
                    %figure(5); clf
                    %imagesc(tmp); colorbar; pause(0.001)
                    %if t == nTrial
                    %    saveas(gca, [img_dir, bands{f}, 'pr_gc_mag.png'], 'png')
                    %end
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
                        saveas(gca, [img_dir, bands{f}, 'pr_gc_grad.png'], 'png')
                    end
                    % save
                    gc_grad(:,t) = GC;
                end
                 if any(any(isnan(gc_mag)))
                    error('This matrix contains nans')
                end
                save([save_dir, 'NMF_', curr_band, '_mag_pr.mat'], 'gc_mag')
                save([save_dir, 'NMF_', curr_band, '_grad_pr.mat'], 'gc_grad')
            end
        catch
            errors{cnte,1} = [d.folder, '/', d.name];
            cnte = cnte + 1;
        end
    end
end
end

