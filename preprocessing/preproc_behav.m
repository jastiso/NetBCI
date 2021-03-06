%% Get behavrioral measures


%% Define Global Variables

addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830')
addpath(genpath('/Users/stiso/Documents/MATLAB/eeglab_current'))

top_dir = '/Users/stiso/Documents/MATLAB/NetBCI/';
data_dir = '/Users/stiso/Resilio Sync/NETBCI.RAW/DataBase/';
sessions = [{'Session1'}, {'Session2'}, {'Session3'}, {'Session4'}];
condition = [{'test01'}, {'test02'}, {'test03'}, {'test04'}, {'test05'}, {'test06'}]; % not including rest
nSubj = 20;


%% Behavioral Analysis

%load data
load([data_dir, '0_Behavior/behavior_20Subjects.mat']);
% Performance - calculated from number of hits and misses, excluding the
% first two trials, which were used to train the classifier
%get vector of performance over time by subject
nRuns = 6;
nSess = 4;
nTp = 15;
perf = zeros(nRuns*nSess*nTp,nSubj); % number of runs, * number of timepoints * number of sessions
perf_real_time = zeros(nRuns*nSess,nSubj);
% for i = 1:nSubj
%     perf(:,i) = interp_perf( behavior_updated, i, nRuns, nTp );
%     
% end
% save([top_dir, 'Behavior/perf.mat'], 'perf')
for i = 1:nSubj
    cnt = 1;
    tmp = behavior_updated.BCI.Perf.Sess1.Runs{i};
    perf_real_time(cnt:cnt+nRuns-1,i) = tmp;
    cnt = cnt+6;
    tmp = behavior_updated.BCI.Perf.Sess2.Runs{i};
    perf_real_time(cnt:cnt+nRuns-1,i) = tmp;
    cnt = cnt+6;
    tmp = behavior_updated.BCI.Perf.Sess3.Runs{i};
    perf_real_time(cnt:cnt+nRuns-1,i) = tmp;
    cnt = cnt+6;
    tmp = behavior_updated.BCI.Perf.Sess4.Runs{i};
    perf_real_time(cnt:cnt+nRuns-1,i) = tmp;
    
end
time = 1:size(perf_real_time,1);
rhos = zeros(1,nSubj);
for i = 1:nSubj
   rhos(i) = corr(perf_real_time(:,i), time', 'type', 'spearman'); 
end
betas = zeros(1,nSubj);
for i = 1:nSubj
   betas(i) = time'\perf_real_time(:,i); 
end
improvement = max(perf_real_time) - min(perf_real_time);
difference = perf_real_time(end,:) - perf_real_time(1,:);
maximum = max(perf_real_time);
final = mean(perf_real_time((end-5):end,:)); % mean of last sessopm
save([top_dir, 'Behavior/stats.mat'], 'betas', 'rhos', 'improvement', 'difference', 'maximum', 'final')

% plot
figure(1)
plot(perf, 'linewidth', 2); hold on
xlabel('Time'); ylabel('Performance')
plot([6,6], [0,100], 'k'); plot([6*2,6*2], [0,100], 'k'); plot([6*3,6*3], [0,100], 'k');
xlim([1, size(perf,1)])
hold off
saveas(gca, [top_dir, 'Behavior/all_subj_perf.png'], 'png')

% get averages
mean_perf = mean(perf,2);
se = std(perf,0,2)/sqrt(nSubj);
figure(2); clf
plot(mean_perf, 'color', rgb('SteelBlue'), 'linewidth', 2); hold on
shade_plot(1:numel(mean_perf), mean_perf', se', rgb('SlateGrey'), 0.5);
%plot(mean_perf-se, '--', 'color', rgb('SlateGrey'), 'linewidth', 1);
%plot(mean_perf+se, '--', 'color', rgb('SlateGrey'), 'linewidth', 1);
plot([nRuns*nTp,nRuns*nTp], [50,80], 'k'); plot([nRuns*nTp*2,nRuns*nTp*2], [50,80], 'k'); plot([nRuns*nTp*3,nRuns*nTp*3], [50,80], 'k');
xlim([0, numel(mean_perf)])
hold off
saveas(gca, [top_dir, 'Behavior/avg_perf.png'], 'png')

