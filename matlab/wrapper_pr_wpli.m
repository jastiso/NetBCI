function [ errors ] = wrapper_pr_wpli( sessions, condition, subj, data_dir, raw_dir, top_dir, bands, freqs, st, en, nNode)
% helper function for paralelizing pipeline
errors = {};
cnte = 1;

for i = 1:numel(sessions)
    session = sessions{i};
    for j = 1:numel(condition)
        cond = condition{j};
        
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
            
            % phase randomize
            for t = 1:nTrial
                data_mag.trial{t} = linsurr(data_mag.trial{t});
                data_grad.trial{t} = linsurr(data_grad.trial{t});
            end
            
            %% Make connectivity matrices
            
            % combine gradiometers
            % until I figure out how to do this, we are just selecting
            % one
            idx = zeros(size(data_grad.label));
            for i = 1:numel(data_grad.label)
                idx(i) = data_grad.label{i}(end) == '2';
            end
            idx = logical(idx);
            data_grad.label = data_grad.label(idx);
            for m = 1:numel(data_grad.trial)
                data_grad.trial{m} = data_grad.trial{m}(idx,:);
            end
            
            for f = 1:numel(bands)
                f_range = freqs(f,1):freqs(f,2);
                curr_band = bands{f};
                
                % get wPLI
                %wpli_mag = zeros((nNode^2-nNode)/2, nTrial);
                wpli_grad = zeros((nNode^2-nNode)/2, nTrial);
                for t = 1:nTrial
                    %                    tmp = zeros(nNode);
                    %                     % magnetometers
                    %                     curr = squeeze(freq_mag.powspctrm(t,:,:));
                    %                     [GC, pairs] = cov_GC(curr,dt,lag,t0);
                    %                     if any(any(isnan(GC)))
                    %                         warning('There are Nans')
                    %                     end
                    %                     GC = sum(GC,2);
                    %                     % plot
                    %                     tmp(logical(tril(ones(nNode),-1))) = GC;
                    %                     tmp = tmp + tmp';
                    %                     figure(5); clf
                    %                     imagesc(tmp); colorbar; pause(0.001)
                    %                     if t == nTrial
                    %                         saveas(gca, [img_dir, bands{f}, 'wpli_mag.png'], 'png')
                    %                     end
                    %                     % save
                    %                     wpli_mag(:,t) = GC;
                    %
                    %gradiometers
                    wpli = get_window_wpli( data_grad, srate, st, en, f_range, t);
                    % plot
                    figure(5); clf
                    imagesc(wpli); colorbar; pause(0.001)
                    if t == nTrial
                        saveas(gca, [img_dir, bands{f}, '_wpli_pr_grad.png'], 'png')
                    end
                    % save
                    v_wpli = reshape(wpli(logical(triu(ones(nNode),1))),1,[]);
                    wpli_grad(:,t) = v_wpli;
                end
                if any(any(isnan(wpli_grad)))
                    error('This matrix contains nans')
                end
                %save([save_dir, 'NMF_', curr_band, '_mag_wpli_pr.mat'], 'wpli_mag')
                save([save_dir, 'NMF_', curr_band, '_grad_wpli_pr.mat'], 'wpli_grad')
            end
        catch
            errors{cnte,1} = [d.folder, '/', d.name];
            cnte = cnte + 1;
        end
    end
end
end

