%% Count triangles
% Count the number of triangles in each Subgraph, and in each individual
% connectivity matrix

addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830')
addpath(genpath('/Users/stiso/Documents/MATLAB/eeglab_current'))

%% First for WPLI data

R_dir = '/Users/stiso/Documents/R/NetBCI/data/wpli/';
top_dir = '/Users/stiso/Documents/MATLAB/NetBCI/';
save_dir = '/Users/stiso/Documents/MATLAB/NetBCI/NMF/';
sessions = [{'Session1'}, {'Session2'}, {'Session3'}, {'Session4'}];
condition = [{'test01'}, {'test02'}, {'test03'}, {'test04'}, {'test05'}, {'test06'}]; % not including rest
subjs = [1:20];
freqs = [3,6;7,14;15,30;31,45];
bands = [{'alpha'},{'beta'},{'low_gamma'}];
eType = 'grad';
nNode = 102;

subj_order = {};
band_order = {};

%possible triangles
C = ones(nNode) - eye(nNode);
pTriangles = sum(diag(C^3));

nTri_wpli = [];
cnt = 0;
for f  = 1:numel(bands)
    freq = bands{f};
    for i = subjs
        A = [];
        subj = sprintf('%03d', i);
        subj_tri = 0;
        subj_tri_bin = 0;
        
        % get labels
        
        cnt = cnt + 1;
        band_order{cnt} = freq;
        subj_order{cnt} = subj;
        
        
        for j = 1:numel(sessions)
            sess = sessions{j};
            fprintf('\nSubj %03d... band %s...', i, freq)
            for k = 1:numel(condition)
                cond = condition{k};
                
                load([top_dir, sess, '/', cond, '/', subj, '/FCmatrices/NMF_', freq, '_', eType, '_wpli.mat']);
                
                % loop through every trial
                for t = 1:size(wpli_grad,2)
                    A = wpli_grad(:,t);
                    % these are fully weighted, and will always have full
                    % number of triangles
                    A = get_sg_matrix(nNode, A);
                    
                    %binarize
                    A(logical(eye(size(A)))) = 0;
                    A_bin = A;
                    A_bin(A_bin~=0) = 1;
                    
                    %weight
                    A = A.^(1/3);
                    
                    % count triangles
                    subj_tri_bin = mean([subj_tri_bin, sum(diag(A_bin^3))./pTriangles]);
                    subj_tri = mean([subj_tri, sum(diag(A^3))]);
                end
            end
        end
        nTri_wpli(cnt) = subj_tri;
        nTri_wpli_bin(cnt) = subj_tri_bin;
    end
end
save([R_dir, 'bin_tri_wpli.mat'], 'subj_order', 'band_order', 'nTri_wpli_bin');
save([R_dir, 'wei_tri_wpli.mat'], 'subj_order', 'band_order', 'nTri_wpli');


%% Now for subgraphs
top_dir = '/Users/stiso/Documents/MATLAB/NetBCI/';
data_dir = '/Users/stiso/Documents/Python/NetBCI/NMF/param/';
save_dir = [top_dir, 'GroupAvg/wpli/analysis/'];
Subj = [1:20];

load([save_dir, 'noise_sg.mat']);
sens = 'grad';
subj_order = {};
band_order = {};

% initialize
high_nTri = [];
high2_nTri = [];
high3_nTri = [];
low_nTri = [];
% binary
high_nTri_bin = [];
high2_nTri_bin = [];
high3_nTri_bin = [];
low_nTri_bin = [];

cnt = 1; % for order
for i = Subj
    subj = sprintf('%03d', i);
    
    for k = 1:numel(bands)
        f = bands{k};
        fprintf('\nSubj %03d... band %s...', i, f)
        
        % get labels
        band_order{cnt} = f;
        subj_order{cnt} = subj;
        
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
        
        % scale to be stable
        % 1st
        high_mat = get_sg_matrix(nNode, subset(bSG,:));
        %2nd
        high2_mat = get_sg_matrix(nNode, subset(bSG2,:));
        %3rd
        high3_mat = get_sg_matrix(nNode, subset(bSG3,:));
        %lowest
        low_mat = get_sg_matrix(nNode, subset(nzSG,:));
        
        % get triangles
        %binarize
        high_mat(logical(eye(size(high_mat)))) = 0;
        high_bin = high_mat;
        high_bin(high_bin~=0) = 1;
        high2_mat(logical(eye(size(high2_mat)))) = 0;
        high2_bin = high2_mat;
        high2_bin(high2_bin~=0) = 1;
        high3_mat(logical(eye(size(high3_mat)))) = 0;
        high3_bin = high3_mat;
        high3_bin(high3_bin~=0) = 1;
        low_mat(logical(eye(size(low_mat)))) = 0;
        low_bin = low_mat;
        low_bin(low_bin~=0) = 1;
        
        % weight
        high_mat = high_mat.^1/3;
        high2_mat = high2_mat.^1/3;
        high3_mat = high3_mat.^1/3;
        low_mat = low_mat.^1/3;
        
        % count triangles
        high_nTri_bin(cnt) = sum(diag(high_bin^3))./pTriangles;
        high2_nTri_bin(cnt) = sum(diag(high2_bin^3))./pTriangles;
        high3_nTri_bin(cnt) = sum(diag(high3_bin^3))./pTriangles;
        low_nTri_bin(cnt) = sum(diag(low_bin^3))./pTriangles;
        
        % count triangles
        high_nTri(cnt) = sum(diag(high_mat^3));
        high2_nTri(cnt) = sum(diag(high2_mat^3));
        high3_nTri(cnt) = sum(diag(high3_mat^3));
        low_nTri(cnt) = sum(diag(low_mat^3));
        
        cnt = cnt + 1;
        
    end
end

%save things
save([R_dir, 'bin_tri_sg.mat'], 'subj_order', 'band_order', 'high_nTri_bin', 'high2_nTri_bin', 'high3_nTri_bin', 'low_nTri_bin');
save([R_dir, 'wei_tri_sg.mat'], 'subj_order', 'band_order', 'high_nTri', 'high2_nTri', 'high3_nTri', 'low_nTri');
