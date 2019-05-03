%% Model Validation

% we are making assumptions that the regularization makes this suitable to
% think of as a pathway for information to flow

% To test this, we will simulate t+1 power based on t power and high/low
% SGs. See if these simulations are more similar when expression is high

addpath(genpath('/Users/stiso/Documents/MATLAB/npy-matlab-master/'))
addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/')
addpath('/Users/stiso/Documents/MATLAB/control_example/')
R_dir = '/Users/stiso/Documents/R/NetBCI/data/wpli/';

top_dir = '/Users/stiso/Documents/MATLAB/NetBCI/';
data_dir = '/Users/stiso/Documents/Python/NetBCI/NMF/param/';
save_dir = [top_dir, 'GroupAvg/wpli/analysis/'];
img_dir = [top_dir, 'GroupAvg/wpli/images/'];
% make directories
if ~exist(img_dir, 'dir')
    mkdir(img_dir)
end
% make directories
if ~exist(save_dir, 'dir')
    mkdir(save_dir)
end
% make directories
if ~exist(R_dir, 'dir')
    mkdir(R_dir)
end

% things for making control set
grad2 = ['MEG0242';'MEG0232';'MEG0442';'MEG0432';'MEG0712';'MEG0742';'MEG1842';'MEG1822';'MEG1812';'MEG1622';'MEG1612'];
grad3 = ['MEG0243';'MEG0233';'MEG0443';'MEG0433';'MEG0713';'MEG0743';'MEG1843';'MEG1823';'MEG1813';'MEG1623';'MEG1613'];

Subj = [1:20];
nSubj = numel(Subj);
freqs = [7,14;15,30;31,45;55,70];
bands = [{'alpha'}, {'beta'}, {'low_gamma'}];
sensors = [{'grad'}];

nNode = 102;
nEdges = (nNode^2-nNode)/2;
load([save_dir, 'noise_sg.mat']);

% initialize
high_corr = cell(nSubj, numel(bands));
low_corr = cell(nSubj, numel(bands));

%% Make control set

% get labels
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

load([top_dir, 'montages/Left_motor_idx.mat'])
B = idx;

% relax B
B(~B) = 1e-5;


%% Loop through data

sens = 'grad';

cnt = 1; % for order
for i = Subj
    s_idx = find(i == Subj);
    subj = sprintf('%03d', i);
    
    for k = 1:numel(bands)
        f = bands{k};
        fprintf('\nSubj %03d... band %s...', i, f)

        % get labels
        band_order{cnt} = f;
        subj_order{cnt} = subj;
        cnt = cnt + 1;

        % get subgraph data
        if strcmp(data_dir(end-5:end-1),'param')
            subset = readNPY([data_dir, subj, '/', sens, 'wpli_', f, '_subset.npy']);
            coeff = readNPY([data_dir, subj, '/',sens, 'wpli_', f, '_coeff.npy']);
        else
            subset = readNPY([data_dir, subj, '/', sens, '/wpli_', f, '_subset.npy']);
            coeff = readNPY([data_dir, subj, '/',sens, '/wpli_', f, '_coeff.npy']);
        end
        
        % remove noise SG
        idx = noise_sg{k,i};
        coeff = coeff(~idx,:);
        subset = subset(~idx,:);
        
        b_exp = subset(:,end);
        [tmp,tmp_idx] = sort(b_exp,'descend');
        bSG3 = tmp_idx(3); % third biggest element
        [~, nzSG] = min(nonzeros((b_exp))); % smallest nonzero - this is how we will opporationalize "low"
        
        
        % scale to be stable
        %3rd
        high3_mat = get_sg_matrix(nNode, subset(bSG3,:));
        high3_mat = high3_mat./eigs(high3_mat,1) - eye(nNode).*1.001;
        %lowest
        low_mat = get_sg_matrix(nNode, subset(nzSG,:));
        low_mat = low_mat./eigs(low_mat,1) - eye(nNode).*1.001;
        
        % load states
        load([top_dir, 'states/' subj, '_states.mat'])
        
        % initialize time seris of correlation values
        nTrial = size(states{k},2);
        high_ts = zeros(1, nTrial - 1);
        low_ts = zeros(1, nTrial - 1);
        
        for t = 2:nTrial
            curr = states{k}(:,(t-1));
            % get u
            u = repmat(ones(nNode,1).*1,1, 1000);

            % get t+1 prediction
            pred_high = openLoopControl(high3_mat, diag(B), curr, u);
            pred_low = openLoopControl(low_mat, diag(B), curr, u);
            %get correlation
            high_ts(t-1) = max(corr(pred_high, states{k}(:,(t)), "type", "pearson"));
            low_ts(t-1) = max(corr(pred_low, states{k}(:,(t)), "type",  "pearson"));
        end
        high_corr{i,k} = high_ts;
        low_corr{i,k} = low_ts;
    end
end

[h,p] = ttest(cellfun(@(x) median(x), high_corr(:,2)));
p
[h,p] = ttest(cellfun(@(x) median(x), low_corr(:,2)));
p
%save([R_dir, 'grad/corrs.mat'], 'band_order', 'subj_order')

%% Correlate with expression

high_exp_corr = zeros(nSubj,numel(bands));
low_exp_corr = zeros(nSubj, numel(bands));
high_max = zeros(nSubj,numel(bands));
low_max = zeros(nSubj, numel(bands));
high_min = zeros(nSubj,numel(bands));
low_min = zeros(nSubj, numel(bands));

sens = 'grad';

cnt = 1; % for order
for i = Subj
    s_idx = find(i == Subj);
    subj = sprintf('%03d', i);
    
    for k = 1:numel(bands)
        f = bands{k};
        fprintf('\nSubj %03d... band %s...', i, f)

        % get labels
        band_order{cnt} = f;
        subj_order{cnt} = subj;
        cnt = cnt + 1;

        % get subgraph data
        if strcmp(data_dir(end-5:end-1),'param')
            subset = readNPY([data_dir, subj, '/', sens, 'wpli_', f, '_subset.npy']);
            coeff = readNPY([data_dir, subj, '/',sens, 'wpli_', f, '_coeff.npy']);
        else
            subset = readNPY([data_dir, subj, '/', sens, '/wpli_', f, '_subset.npy']);
            coeff = readNPY([data_dir, subj, '/',sens, '/wpli_', f, '_coeff.npy']);
        end
        
        % remove noise SG
        idx = noise_sg{k,i};
        coeff = coeff(~idx,:);
        subset = subset(~idx,:);

        [tmp,tmp_idx] = sort(b_exp,'descend');
        bSG3 = tmp_idx(3); % third biggest element
        [~, nzSG] = min(nonzeros((b_exp))); % smallest nonzero - this is how we will opporationalize "low"
        
        high_exp = coeff(bSG3,2:end);
        low_exp = coeff(nzSG,2:end);
        
        high_exp_corr(i,k) = corr(high_exp', high_corr{i,k}');
        low_exp_corr(i,k) = corr(low_exp', low_corr{i,k}');
        
        % max and min
        [~,maxI_high] = max(high_exp);
        [~,minI_high] = min(high_exp);
        [~,maxI_low] = max(high_exp);
        [~,minI_low] = min(high_exp);
        
        high_max(i,k) = high_corr{i,k}(maxI_high);
        low_max(i,k) = low_corr{i,k}(maxI_low);
        high_min(i,k) = high_corr{i,k}(minI_high);
        low_min(i,k) = low_corr{i,k}(minI_low);
    end
end

figure(1); clf;
imagesc(high_exp_corr); colorbar
saveas(gca, [img_dir, 'high_te_corr.png'], 'png')

figure(2); clf;
imagesc(low_exp_corr); colorbar
saveas(gca, [img_dir, 'low_te_corr.png'], 'png')

[h,p] = ttest(high_exp_corr(:,2));
p
[h,p] = ttest(low_exp_corr(:,2));
p

diff_high = high_max - high_min;
diff_low = low_max - low_min;

figure(3); clf;
imagesc(diff_high); colorbar
saveas(gca, [img_dir, 'high_max_min_corr.png'], 'png')

figure(4); clf;
imagesc(diff_low); colorbar
saveas(gca, [img_dir, 'low_max_min_corr.png'], 'png')

[h,p] = ttest(high_max(:,2), high_min(:,2));
p
[h,p] = ttest(low_max(:,2), low_min(:,2));
p

%% Repeat for switched expression

high_exp_corr_switch = zeros(nSubj,numel(bands));
low_exp_corr_switch = zeros(nSubj, numel(bands));
high_max = zeros(nSubj,numel(bands));
low_max = zeros(nSubj, numel(bands));
high_min = zeros(nSubj,numel(bands));
low_min = zeros(nSubj, numel(bands));

sens = 'grad';

cnt = 1; % for order
for i = Subj
    s_idx = find(i == Subj);
    subj = sprintf('%03d', i);
    
    for k = 1:numel(bands)
        f = bands{k};
        fprintf('\nSubj %03d... band %s...', i, f)

        % get labels
        band_order{cnt} = f;
        subj_order{cnt} = subj;
        cnt = cnt + 1;

        % get subgraph data
        if strcmp(data_dir(end-5:end-1),'param')
            subset = readNPY([data_dir, subj, '/', sens, 'wpli_', f, '_subset.npy']);
            coeff = readNPY([data_dir, subj, '/',sens, 'wpli_', f, '_coeff.npy']);
        else
            subset = readNPY([data_dir, subj, '/', sens, '/wpli_', f, '_subset.npy']);
            coeff = readNPY([data_dir, subj, '/',sens, '/wpli_', f, '_coeff.npy']);
        end
        
        % remove noise SG
        idx = noise_sg{k,i};
        coeff = coeff(~idx,:);
        subset = subset(~idx,:);

        [tmp,tmp_idx] = sort(b_exp,'descend');
        bSG3 = tmp_idx(3); % third biggest element
        [~, nzSG] = min(nonzeros((b_exp))); % smallest nonzero - this is how we will opporationalize "low"
        
        high_exp = coeff(bSG3,2:end);
        low_exp = coeff(nzSG,2:end);
        
        % switched exp
        high_exp_corr_switch(i,k) = corr(low_exp', high_corr{i,k}');
        low_exp_corr_switch(i,k) = corr(high_exp', low_corr{i,k}');
        
        % max and min
        [~,maxI_high] = max(high_exp);
        [~,minI_high] = min(high_exp);
        [~,maxI_low] = max(high_exp);
        [~,minI_low] = min(high_exp);
        
        high_max(i,k) = high_corr{i,k}(maxI_high);
        low_max(i,k) = low_corr{i,k}(maxI_low);
        high_min(i,k) = high_corr{i,k}(minI_high);
        low_min(i,k) = low_corr{i,k}(minI_low);
    end
end

figure(1); clf;
imagesc(high_exp_corr_switch); colorbar
saveas(gca, [img_dir, 'high_te_corr_switch.png'], 'png')

figure(2); clf;
imagesc(low_exp_corr_switch); colorbar
saveas(gca, [img_dir, 'low_te_corr_switch.png'], 'png')

[h,p] = ttest(high_exp_corr_switch(:,2));
p
[h,p] = ttest(low_exp_corr_switch(:,2));
p

[h,p] = ttest2(high_exp_corr_switch(:,2), high_exp_corr(:,2));
p
[h,p] = ttest2(low_exp_corr_switch(:,2), low_exp_corr(:,2));
p

diff_high = high_max - high_min;
diff_low = low_max - low_min;

figure(3); clf;
imagesc(diff_high); colorbar
saveas(gca, [img_dir, 'high_max_min_corr.png'], 'png')

figure(4); clf;
imagesc(diff_low); colorbar
saveas(gca, [img_dir, 'low_max_min_corr.png'], 'png')

[h,p] = ttest(high_max(:,2), high_min(:,2));
p
[h,p] = ttest(low_max(:,2), low_min(:,2));
p
