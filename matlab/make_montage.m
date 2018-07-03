% make montages into vectors


addpath(genpath('/Users/stiso/Documents/MATLAB/npy-matlab-master/'))
addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/')
R_dir = '/Users/stiso/Documents/R/NetBCI/data/gc/';

top_dir = '/Users/stiso/Documents/MATLAB/NetBCI/';
data_dir = '/Users/stiso/Resilio Sync/NETBCI.RAW/DataBase/1_Signals/0_RawData/2_MEG/';
mon_dir = '/Users/stiso/Documents/MATLAB/NetBCI/montages/';

regions = [{'Left_frontal'}, {'Left_occipital'}, {'Left_parietal'}, {'Left_temporal'}, ...
    {'Right_frontal'}, {'Right_occipital'}, {'Right_parietal'}, {'Right_temporal'}, {'Vertex'}];

% load raw data for labels
load([data_dir, 'Session1/1_test01/RawData_MEG_Subj001_Ses1_test01.mat'])
labels = RawData_MEG_Subj001_Ses1_test01.label;

% get grad labels
mag_idx = strcmp(RawData_MEG_Subj001_Ses1_test01.meg_sensors.chantype,'megmag');
grad_idx = ~mag_idx;
grad_label = labels(grad_idx);
cmb_label = [];
cnt = 1;
for i = 1:2:numel(grad_label)
    cmb_label = [cmb_label; grad_label{i}, '+', grad_label{i+1}(end-3:end)];
    cnt = cnt + 1;
end

% get vector for every montage
for i = 6%1:numel(regions)
   load([mon_dir, 'montage_', regions{i}, '_MEG.mat'])
   tmp_idx = false(size(cmb_label,1),2);
   for j = i:numel(Montages.ChanNames)
       curr = Montages.ChanNames{j};
       tmp_idx(:,1) = tmp_idx(:,1) | sum(curr((end-3):end) == char(cmb_label(:,(end-3):end)),2) == 4 ;
       tmp_idx(:,2) = tmp_idx(:,2) | sum(curr((end-3):end) == char(cmb_label(:,4:7)),2) == 4 ;
   end
   idx = tmp_idx(:,1) | tmp_idx(:,2);
   save([mon_dir, regions{i}, '_idx.mat'], 'idx')
end