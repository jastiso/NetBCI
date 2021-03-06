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
pred_states = cell(nSubj, numel(bands));
pred_temp_exp_null = cell(nSubj, numel(bands));
pred_temp_exp = cell(nSubj, numel(bands));
pred_temp_exp_3 = cell(nSubj, numel(bands));
pred_temp_exp_low = cell(nSubj, numel(bands));
pred_temp_exp_3_only = cell(nSubj, numel(bands));
pred_temp_exp_low_only = cell(nSubj, numel(bands));

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
        
        nSG = numel(subset(:,end));
        b_exp = subset(:,end);
        [tmp,tmp_idx] = sort(b_exp,'descend');
        bSG3 = tmp_idx(3); % third biggest element
        [~, nzSG] = min(nonzeros((b_exp))); % smallest nonzero - this is how we will opporationalize "low"
       
        
        % scale to be stable
        all_sg = zeros(nNode,nNode,nSG);
        for n = 1:nSG
            curr_mat = get_sg_matrix(nNode, subset(n,:));
            all_sg(:,:,n) = curr_mat./eigs(curr_mat,1) - eye(nNode).*1.001;
        end
        
        
        % load states
        load([top_dir, 'states/' subj, '_states.mat'])
        
        % initialize time seris of correlation values
        nTrial = size(states{k},2);
        ts = zeros(nNode, nTrial);
        
        curr_pred = zeros(nNode, nSG, nTrial);
        for n = 1:nSG
            curr_mat = all_sg(:,:,n);
            for t = 2:nTrial
                curr = states{k}(:,(t-1));
                % get u
                u = repmat(ones(nNode,1).*1,1, 1000);
                
                % get t+1 prediction
                pred = openLoopControl(curr_mat, diag(B), curr, u);
                
                %get correlation
                [~,idx] = max(corr(pred, states{k}(:,(t)), "type", "pearson"));
                ts(:,t) = pred(:,idx);
            end
            curr_pred(:,n,:) = ts;
        end
        pred_states{i,k} = curr_pred;
        
        % get predicted temporal coeffs
        pte = zeros(nSG,nTrial);
        pte_3 = zeros(1,nTrial);
        pte_low = zeros(1,nTrial);
        for t = 2:nTrial
            curr = states{k}(:,(t-1));
            A = pred_states{i,k}(:,:,t);
            A_3 = pred_states{i,k}(:,bSG3,t);
            A_low = pred_states{i,k}(:,nzSG,t);
            
            % get x_dot
            x_dot = states{k}(:,(t)) - curr;
            
            %get predicted temp_exp
            pte(:,t) = pinv(A)*x_dot; 
            pte_3(:,t) = pinv(A_3)*x_dot; 
            pte_low(:,t) = pinv(A_low)*x_dot; 
        end
        
        % zscore across graphs
        pte = zscore(pte,0,1);
        coeff1 = zscore(coeff,0,1);
        
        % get correlations
        curr_corr = [];
        curr_corr_null = [];
        for n = 1:nSG
           r=n;
           while r==n
               r = randi(nSG);
           end
           curr_corr = [curr_corr, corr((pte(n,:)'), (coeff(n,:)'))];
           curr_corr_null = [curr_corr_null, corr((pte(n,:)'), (coeff(r,:)'))];
        end
        curr_corr_3 = corr((pte(bSG3,:)'), (coeff(bSG3,:)'));
        curr_corr_low = corr((pte(nzSG,:)'), (coeff(nzSG,:)'));
        
        curr_corr_3_only = corr((pte_3'), (coeff(bSG3,:)'));
        curr_corr_low_only = corr((pte_low'), (coeff(nzSG,:)'));
        
        pred_temp_exp{i,k} = mean(curr_corr);
        pred_temp_exp_null{i,k} = mean(curr_corr_null);
        pred_temp_exp_3{i,k} = curr_corr_3;
        pred_temp_exp_low{i,k} = curr_corr_low;
        
        pred_temp_exp_3_only{i,k} = curr_corr_3_only;
        pred_temp_exp_low_only{i,k} = curr_corr_low_only;
    end
end


%save([R_dir, 'grad/corrs.mat'], 'band_order', 'subj_order')

%% Plot

imagesc(cell2mat(pred_temp_exp))

empirical = [pred_temp_exp{:,2}];
p = mult_comp_perm_t1(empirical')

save([R_dir, 'model_validation.mat'], 'empirical', 'null', 'low', 'high3')


