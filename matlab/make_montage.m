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
for i = 1:numel(regions)
   load([mon_dir, 'montage_', regions{i}, '_MEG.mat'])
   tmp_idx = false(size(cmb_label,1),2);
   for j = i:numel(Montages.ChanNames)
       curr = Montages.ChanNames{j};
       % if both of the gradiometers are in the lobe, it goes in the index
       tmp_idx(:,1) = tmp_idx(:,1) | (sum(curr((end-3):end) == char(cmb_label(:,(end-3):end)),2) == 4);
       tmp_idx(:,2) = tmp_idx(:,2) | (sum(curr((end-3):end) == char(cmb_label(:,4:7)),2) == 4);
   end
   % for a couple montages, this gave the same answer with or or and. For
   % some, its off but only by a little
   idx = tmp_idx(:,1) & tmp_idx(:,2);
   save([mon_dir, regions{i}, '_idx.mat'], 'idx')
end

%% make border indicies

fronto_parietal = [{'MEG0212+0213'}, {'MEG0222+0223'},  {'MEG0632+0633'}, {'MEG1122+1123'}, {'MEG1312+1313'}, {'MEG1322+1323'}, {'MEG1442+1443'},  {'MEG1422+1423'}];
parieto_occipital = [{'MEG2042+2043'}, {'MEG2032+2033'}, {'MEG2312+2313'}, {'MEG2342+2343'}];
midline_frontal = [{'MEG0812+0813'}, {'MEG0912+0913'}];

idx = ismember(cmb_labels, fronto_parietal);
save([mon_dir, 'Frontoparietal_idx.mat'], 'idx')

idx = ismember(cmb_labels, parieto_occipital);
save([mon_dir, 'Parietoccipital_idx.mat'], 'idx')

idx = ismember(cmb_labels, midline_frontal);
save([mon_dir, 'Frontal_idx.mat'], 'idx')