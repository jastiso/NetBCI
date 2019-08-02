%% Get trial numbers
% this script is made solely to address an apparent discrepancy inthe
% number of trials between that MEG data, and the behavioral data.

%% Load things

addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830')
addpath(genpath('/Users/stiso/Documents/MATLAB/eeglab_current'))
addpath(genpath('/Users/stiso/Documents/Code/NetBCI/matlab/'))

top_dir = '/Users/stiso/Documents/MATLAB/NetBCI/';
data_dir = '/Users/stiso/Resilio Sync/NETBCI.RAW/DataBase/1_Signals/2_Segmentation/2_MEG/';
raw_dir = '/Users/stiso/Resilio Sync/NETBCI.RAW/DataBase/1_Signals/0_RawData/2_MEG/';
sessions = [{'Session1'}, {'Session2'}, {'Session3'}, {'Session4'}];
condition = [{'test01'}, {'test02'}, {'test03'}, {'test04'}, {'test05'}, {'test06'}]; % not including rest
subjs = [1:20];

res = 2; % number of entries per run (15 or 16 trials)
behavior_all = cell(numel(subjs),1);
errors = [];
% load behavior
load([top_dir, 'Behavior/behavior_updated_trials.mat']);

trials = zeros(2,20); % number of modalities (MEG,behavior), number of subjects
trials_breakdown = cell(20,2);
diff = {};

%% Behavior

cnte = 0;
for i = subjs
    for j = 1:numel(sessions)
        sess = ['Sess', num2str(j)];
        for k = 1:numel(condition)
            eval(['curr = behavior_updated.BCI.Perf.', sess, '.Trials.Result{i}(k,:);']);
            eval(['idx = behavior_updated.BCI.Perf.', sess, '.Trials.TargetCode{i}(k,:);']);
            curr = curr(idx == 2);
            if numel(curr) > 16
                warning('too big')
            end
            % get number of trials in each
            nTrials = numel(curr);
            trials(1,i) = trials(1,i) + nTrials;
            trials_breakdown{i,1} = [trials_breakdown{i,1}, nTrials];
        end
    end
end

%% MEG

for i = 1:numel(sessions)
    session = sessions{i};
    for j = 1:numel(condition)
        cond = condition{j};
        for k = subjs
            subj = k;
            ext = sprintf('%03d',subj);
            disp(ext)
            % load data
            d = dir([data_dir, session, '/*', cond, '*', '/Seg_MEG_Subj', ext, '_Ses', session(end), '_', cond, '.mat']);
            d_raw = dir([raw_dir, session, '/*', cond, '*', '/RawData_MEG_Subj', ext, '_Ses', session(end), '_', cond, '.mat']);
            
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
            
            trials(2,k) = trials(2,k) + numel(data.trial);
            trials_breakdown{k,2} = [trials_breakdown{k,2}, numel(data.trial)];
            % flag if something is different
            if numel(data.trial) ~= trials_breakdown{k,1}(numel(trials_breakdown{k,2}))
                diff{numel(diff)+1} = ['Subj', ext, '_Ses', session(end), '_', cond];
            end
        end
    end
end
