%% Temporal Exp stats for high, high2, high3, low

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
regions = [{'Left_frontal'}, {'Left_occipital'}, {'Left_parietal'}, {'Left_temporal'}, ...
    {'Right_frontal'}, {'Right_occipital'}, {'Right_parietal'}, {'Right_temporal'}, {'Vertex'}, {'Left_motor'}, {'Right_motor'}];
%regions = [{'Left_motor'}, {'Right_motor'}];

% load behavior
load([top_dir, 'Behavior/stats'])

nNode = 102;
nEdges = (nNode^2-nNode)/2;
load([save_dir, 'noise_sg.mat']);

% initialize
e_high = [];
e_high2 = [];
e_high3 = [];
e_low = [];
e_zero = [];
band_order = {};
subj_order = {};
slope = [];
p_high = [];
p_high2 = [];
p_high3 = [];
p_low = [];
p_zero = [];
m_high = [];
m_high2 = [];
m_high3 = [];
m_low = [];
m_zero = [];

%% Loop through data

% this used to index over mag and grad, but we dont look at mag
sens = sensors{1};
R_dir_s = [R_dir, sens, '/'];

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
        slope(cnt) = betas(i);
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
        [~,tmp_idx] = sort(b_exp,'descend');
        bSG = tmp_idx(1);
        bSG2 = tmp_idx(2); % second biggest element
        bSG3 = tmp_idx(3); % third biggest element
        [~, nzSG] = min(nonzeros((b_exp))); % smallest nonzero - this is how we will opporationalize "low"
        
        % if you actualy want to look at 0s
        if sum(b_exp == 0) <= 1
            [~,nbSG] = min(b_exp);
        else
            nbSG = find(b_exp == 0);
        end
        nZero = numel(nbSG);
        
        % get energy
        e_high = [e_high, sum(coeff(bSG,:).^2)];
        e_high2 = [e_high2, sum(coeff(bSG2,:).^2)];
        e_high3 = [e_high3, sum(coeff(bSG3,:).^2)];
        e_low = [e_low, sum(coeff(nzSG,:).^2)];
        tmp = zeros(1,nZero);
        for j = 1:nZero
            tmp(j) = sum(coeff(nbSG(j),:).^2);
        end
        e_zero = [e_zero, mean(tmp)];
        
        % get index of max and max
        [m,n] = max(coeff(bSG,:));
        p_high = [p_high, n/numel(coeff(bSG,:))];
        m_high = [m_high, m];
        [m,n] = max(coeff(bSG2,:));
        p_high2 = [p_high2, n/numel(coeff(bSG2,:))];
        m_high2 = [m_high2, m];
        [m,n] = max(coeff(bSG3,:));
        p_high3 = [p_high3, n/numel(coeff(bSG3,:))];
        m_high3 = [m_high3, m];
        [m,n] = max(coeff(nzSG,:));
        p_low = [p_low, n/numel(coeff(nzSG,:))];
        m_low = [m_low, m];
        tmp = zeros(2,nZero);
        for j = 1:nZero
            [m,n] = max(coeff(nbSG(j),:));
            tmp(1,j) = n/numel(coeff(nbSG(j),:));
            tmp(2,j) = m;
        end
        p_zero = [p_zero, mean(tmp(1,:))];
        m_zero = [m_zero, mean(tmp(2,:))];
    end
end

save([R_dir, 'grad/temp_exp.mat'], 'band_order', 'subj_order', 'slope', 'e_high', 'e_high2', 'e_high3', 'e_low', 'e_zero',...
    'p_high', 'p_high2', 'p_high3', 'p_low', 'p_zero', 'm_high', 'm_high2', 'm_high3', 'm_low', 'm_zero')
