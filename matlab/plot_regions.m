
addpath(genpath('/Users/stiso/Documents/MATLAB/npy-matlab-master/'))
addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/')

top_dir = '/Users/stiso/Documents/MATLAB/NetBCI/';
img_dir = [top_dir, 'GroupAvg/wpli/images/'];

cfg.layout = [top_dir, 'layouts/neuromag306cmb.lay'];
regions = [{'Left_frontal'}, {'Left_occipital'}, {'Left_parietal'}, {'Left_temporal'}, ...
    {'Right_frontal'}, {'Right_occipital'}, {'Right_parietal'}, {'Right_temporal'}, {'Vertex'}];
nNode = 102;

%% Lobes

for i = 1:numel(regions)
    load([top_dir, 'montages/', regions{i}, '_idx.mat'])
    plot_data.powspctrm = idx;
    plot_data.label = cmb_labels;
    plot_data.dimord = 'chan_freq';
    plot_data.freq = 1;
    plot_data.cfg = [];
    figure(1); clf
    ft_topoplotER(cfg,plot_data); colorbar
    saveas(gca, [img_dir, regions{i}, '.png'], 'png')
end

%% Motor

% things for making control set
l_grad2 = ['MEG0242';'MEG0232';'MEG0442';'MEG0432';'MEG0712';'MEG0742';'MEG1842';'MEG1822';'MEG1812';'MEG1622';'MEG1612'];
l_grad3 = ['MEG0243';'MEG0233';'MEG0443';'MEG0433';'MEG0713';'MEG0743';'MEG1843';'MEG1823';'MEG1813';'MEG1623';'MEG1613'];

r_grad2 = ['MEG1332';'MEG1342';'MEG1132';'MEG1142';'MEG0722';'MEG0732';'MEG2212';'MEG2222';'MEG2412';'MEG2232';'MEG2422'];
r_grad3 = ['MEG1333';'MEG1343';'MEG1133';'MEG1143';'MEG0723';'MEG0733';'MEG2213';'MEG2223';'MEG2413';'MEG2233';'MEG2423'];

load(['/Users/stiso/Resilio Sync/NETBCI.RAW/DataBase/1_Signals/2_Segmentation/2_MEG/1_Signals/2_Segmentation/2_MEG/Session1/1_test01/Seg_MEG_Subj001_Ses1_test01.mat'])
labels = Seg_MEG_Subj001_Ses1_test01.MI.label;
load(['/Users/stiso/Resilio Sync/NETBCI.RAW/DataBase/1_Signals/0_RawData/2_MEG/Session1/1_test01/RawData_MEG_Subj001_Ses1_test01.mat'])
mag_idx = strcmp(RawData_MEG_Subj001_Ses1_test01.meg_sensors.chantype,'megmag');
grad_label = labels(~mag_idx);

cnt = 1;
cmb_labels = cell(nNode,1);
for i = 1:2:numel(grad_label)
    cmb_labels{cnt} = [grad_label{i}, '+', grad_label{i+1}(end-3:end)];
    cnt = cnt + 1;
end
lm_label = cell(size(l_grad2,1),i);
for i = 1:size(l_grad2,1)
    lm_label{i} = [l_grad2(i,:), '+', l_grad3(i,end-3:end)];
end
rm_label = cell(size(r_grad2,1),i);
for i = 1:size(r_grad2,1)
    rm_label{i} = [r_grad2(i,:), '+', r_grad3(i,end-3:end)];
end

% get idx
lm = zeros(nNode,1);
for i = 1:size(l_grad2,1)
    curr_idx = find(contains(cmb_labels,lm_label{i}));
    lm(curr_idx) = 1;
end
idx = lm;
save('/Users/stiso/Documents/MATLAB/NetBCI/montages/Left_motor_idx.mat', 'idx');
rm = zeros(nNode,1);
for i = 1:size(r_grad2,1)
    curr_idx = find(contains(cmb_labels,rm_label{i}));
    rm(curr_idx) = 1;
end
idx = rm;
save('/Users/stiso/Documents/MATLAB/NetBCI/montages/Right_motor_idx.mat', 'idx');

% plot

plot_data.powspctrm = lm;
plot_data.label = cmb_labels;
plot_data.dimord = 'chan_freq';
plot_data.freq = 1;
plot_data.cfg = [];
figure(1); clf
ft_topoplotER(cfg,plot_data); colorbar
saveas(gca, [img_dir, 'lm.png'], 'png')

plot_data.powspctrm = rm;
plot_data.label = cmb_labels;
plot_data.dimord = 'chan_freq';
plot_data.freq = 1;
plot_data.cfg = [];
figure(2); clf
ft_topoplotER(cfg,plot_data); colorbar
saveas(gca, [img_dir, 'rm.png'], 'png')