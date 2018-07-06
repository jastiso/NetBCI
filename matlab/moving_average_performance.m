%% Get a moving average of the BCI performance at different resolutions

% @author JStiso jeni.stiso@gmail.com
% Change Log
% July 5 - wrote script
%

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

res = 1; % number of entries per run (15 or 16 trials)
behavior_all = cell(numel(subjs),1);
errors = [];
% load behavior
load([top_dir, 'Behavior/behavior_updated_trials.mat']);

%% get moving average

cnte = 0;
for i = subjs
    for j = 1:numel(sessions)
        sess = ['Sess', num2str(j)];
        for k = 1:numel(condition)
            try
                eval(['curr = behavior_updated.BCI.Perf.', sess, '.Trials.Result{i}(k,:);']);
                eval(['idx = behavior_updated.BCI.Perf.', sess, '.Trials.TargetCode{i}(k,:);']);
                curr = curr(idx == 1);
                % get number of trials in each mvoing average thing
                nTrials = numel(curr);
                nToAvg = zeros(1,res);
                nToAvg(1:(end-1)) = floor(nTrials/res);
                nToAvg(end) = nTrials - sum(nToAvg);
                m_avg = zeros(1,nTrials);
                cnt = 1;
                for m = 1:numel(nToAvg)
                    m_avg(cnt:(cnt+nToAvg(m)-1)) = mean(curr(cnt:(cnt+nToAvg(m)-1)));
                    cnt = cnt + nToAvg(m);
                end
                behavior_all{i} = [behavior_all{i}, m_avg*100];
            catch
                cnte = cnte + 1;
                errors{cnte} = ['Subj', num2str(i), '_' sess];
            end
        end
    end
end
save([top_dir, 'Behavior/behavior_all_movavg.mat'], 'behavior_all');

% plot
b_mat = [];
figure(1); clf
hold on
for i = subjs
    b_mat = [b_mat; unique(behavior_all{i})];
    plot(behavior_all{i}, 'linewidth', 2)
    %hist(behavior_all{i})
end
hold off



