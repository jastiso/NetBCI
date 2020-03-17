%% Get EOG connectivity
% make consensus graph for time varying EOG connectivity

%% Global Variables

addpath(genpath('/Users/stiso/Documents/MATLAB/npy-matlab-master/'))
addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/')
addpath(genpath('/Users/stiso/Documents/MATLAB/BCT/'))

top_dir = '/Users/stiso/Documents/MATLAB/NetBCI/';
data_dir = '/Users/stiso/Documents/Python/NetBCI/NMF/';
save_dir = [top_dir, 'GroupAvg/wpli/analysis/'];
img_dir = [top_dir, 'GroupAvg/wpli/images/'];

subjs = [1:20];
nSubj = numel(subjs);
nNode = 102;
nEdge = (nNode^2-nNode)/2;
freqs = [7,14;15,30;31,45];
bands = [{'alpha'},{'beta'},{'low_gamma'}];
sensors = [{'grad'}];
sessions = [{'Session1'}, {'Session2'}, {'Session3'}, {'Session4'}];
condition = [{'test01'}, {'test02'}, {'test03'}, {'test04'}, {'test05'}, {'test06'}]; 
nBio = 2;

% initialize
consensus_eog = zeros(nNode,nBio, nSubj, numel(bands));


% get labels for topoplot
labels = [];
load(['/Users/stiso/Resilio Sync/NETBCI.RAW/DataBase/1_Signals/0_RawData/2_MEG/Session1/1_test01/RawData_MEG_Subj001_Ses1_test01.mat'])
data = RawData_MEG_Subj001_Ses1_test01;
clear RawData_MEG_Subj001_Ses1_test01

sens_idx = strcmp(data.meg_sensors.chantype,'megplanar');
grad_label = data.label(sens_idx);
cnt = 1;
for n = 1:2:sum(sens_idx)
    labels{cnt,1} = [grad_label{n}, '+', grad_label{n+1}(end-3:end)];
    cnt = cnt + 1;
end

%% Loop through data

sens = sensors{1};

for i = subjs
    subj = sprintf('%03d', i);
    load([top_dir, 'idx_', sens, '.mat'])
    for k = 1:numel(bands)
        f = bands{k};
        fprintf('\nSubj %03d... band %s...', i, f)
        
        wpli_all = [];
        for n = 1:numel(sessions)
            session = sessions{n};
            for m = 1:numel(condition)
                cond = condition{m};
                % get eog wpli data for all trials
                load([top_dir, session, '/', cond, '/', subj, '/FCmatrices/', f, '_bio_wpli.mat'])
                wpli_all = [wpli_all, wpli_bio];
            end
        end
        wpli_all = mean(wpli_all,2);
        wpli_all = reshape(wpli_all, nNode, nBio);
        % add to consensus
        consensus_eog(:,:,i,k) = wpli_all;
    end
end


%% make and plot consensus - consistency

thr = [.75];
sum_edge = zeros(numel(thr), numel(bands));
max_edge = zeros(numel(thr), numel(bands));

for i = 1:numel(bands)
    curr = consensus_eog(:,:,:,i);
    
    for j = 1:numel(thr)
        c = zeros(nNode,nBio, nSubj);

        % threshold
        for k = 1:nSubj
            for n = 1:nBio
                A = curr(:,n,k);
                [edges] = sort(nonzeros(A),'ascend');
                nEdge = round(numel(edges)*thr);
                rm_edges = edges(1:nEdge);
                for m = 1:numel(rm_edges)
                    A(A == rm_edges(m)) = 0;
                end
                c(:,n,k) = A;
            end
        end
        % get consistency
        c_avg = sum(c ~=0, 3); % ~=0 for consistency
        
        % get number of edges remaining and max
        sum_edge(j,i) = sum(sum(c_avg))/numel(c_avg); % same number in each sg
        max_edge(j,i) = max(max(c_avg));
        
        % save consistent mats for strictest thresh
        if j == 3
            save([save_dir, 'consistent_bio_', bands{i}, '.mat'], 'c_avg')
        end
        
        % add some additional thresholding
        %c_avg(c_avg <= 5) = 0;
        
        % save stricter thresholding to report later
        [r,ci] = find(c_avg > 9);
        
        % topoplot
        for n =1:nBio
            plot_data.powspctrm = c_avg(:,n);
            plot_data.label = labels;
            plot_data.dimord = 'chan_freq';
            plot_data.freq = mean(freqs(i,:));
            plot_data.cfg = [];
            cfg = [];
            cfg.layout = [top_dir, 'layouts/neuromag306cmb.lay'];
            figure(1); clf
            ft_topoplotER(cfg,plot_data); colorbar
            saveas(gca, [img_dir, bands{i}, '_eog_topo', num2str(n), '.png'], 'png')
        end
    end
end

save([save_dir, 'consistent_bio.mat'], 'sum_edge', 'max_edge')
