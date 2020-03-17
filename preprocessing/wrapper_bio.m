function [ errors ] = wrapper_bio( sessions, condition, subj, bio_dir, data_dir, raw_dir, top_dir, bands, freqs, st, en, nNode)
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
        d_bio = dir([bio_dir, session, '/*', cond, '*', '/Seg_BIO_Subj', ext, '_Ses', session(end), '_', cond, '.mat']);
        %try
            load([d.folder, '/', d.name]);
            load([d_bio.folder, '/', d_bio.name]);
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
            try
                eval(['data_bio = ', d_bio.name(1:end-4), ';']);
                eval(['clear ', d_bio.name(1:end-4)]);
            catch
                eval(['data_bio = ', d_bio.name(1:13), '16', d_bio.name(16:(end-4)), ';']);
                eval(['clear ', d_bio.name(1:13), '16', d_bio.name(16:(end-4))]);
            end
            bio = data_bio.MI;   
            
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
            nBio = 2; %number of EOG channels
            
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
            
            % combine bio with brain data
            for t = 1:nTrial
                data_grad.trial{t} = [data_grad.trial{t}; bio.trial{t}(4:5,:)]; % EOG channels are the last two
            end
            data_grad.label = [data_grad.label(:)', bio.label(4:5,:)'];
            bio_idx = zeros(nNode+nBio,1);
            bio_idx(nNode+1:end) = 1;
            sens_idx = ~bio_idx;
            
            for f = 1:numel(bands)
                f_range = freqs(f,1):freqs(f,2);
                curr_band = bands{f};
                
                % get wPLI
                wpli_bio = zeros((nNode*nBio), nTrial);
                for t = 1:nTrial
                    %gradiometers
                    wpli = get_window_wpli( data_grad, srate, st, en, f_range, t);
                    
                    % plot
%                     plot_data.powspctrm = wpli(5, logical(sens_idx))';
%                     plot_data.label = data_grad.label(sens_idx);
%                     if t == nTrial
%                         saveas(gca, [img_dir, bands{f}, '_wpli_bio.png'], 'png')
%                     end
                    % save
                    v_wpli = reshape(wpli(logical(bio_idx), logical(sens_idx)), [], 1);
                    wpli_bio(:,t) = v_wpli;
                end
                if any(any(isnan(wpli_bio)))
                    error('This matrix contains nans')
                end
                save([save_dir, curr_band, '_bio_wpli.mat'], 'wpli_bio')
            end
%          catch
%              errors{cnte,1} = [d.folder, '/', d.name];
%              cnte = cnte + 1;
%          end
    end
end
end

