%% Subgraph lobe centrality
% abs(# edges in lobe - # edges between lobes)/# edges within lobe + #
% edges between lobes
% 
% CHange Log
% June 19th, fixed bug thatused coef instead of subset to select max and
% min graphs

addpath(genpath('/Users/stiso/Documents/MATLAB/npy-matlab-master/'))
addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/')
R_dir = '/Users/stiso/Documents/R/NetBCI/data/gc/';

top_dir = '/Users/stiso/Documents/MATLAB/NetBCI/';
data_dir = '/Users/stiso/Documents/Python/NetBCI/NMF/gc/';
save_dir = [top_dir, 'GroupAvg/gc/analysis/'];
img_dir = [top_dir, 'GroupAvg/gc/images/'];
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

Subj = [1:18,20];
nSubj = numel(Subj);
freqs = [7,14;15,30;31,45;55,70];
bands = [{'alpha'}, {'beta'}, {'low_gamma'}];
sensors = [{'grad'}];
regions = [{'Left_frontal'}, {'Left_occipital'}, {'Left_parietal'}, {'Left_temporal'}, ...
    {'Right_frontal'}, {'Right_occipital'}, {'Right_parietal'}, {'Right_temporal'}, {'Vertex'}];

nNode = 102;
nEdges = (nNode^2-nNode)/2;
edge_similarity = zeros(nEdges, nSubj-1, numel(bands), numel(sensors));
% initialize
high_participation = cell(nSubj*numel(bands)*numel(regions), 4);
low_participation = cell(nSubj*numel(bands)*numel(regions), 4);
% load behavior
load([top_dir, 'Behavior/stats'])

%% Loop through data
for j = 1:numel(sensors)
    sens = sensors{j};
    R_dir_s = [R_dir, sens, '/'];
    
    cnth = 0;
    cntl = 0;
    for i = Subj
        s_idx = find(i == Subj);
        subj = sprintf('%03d', i);
        
        for k = 1:numel(bands)
            f = bands{k};
            
            % get subgraph data
            subset = readNPY([data_dir, subj, '/', sens, '/gc_', f, '_subset.npy']);
            coeff = readNPY([data_dir, subj, '/',sens, '/gc_', f, '_coeff.npy']);
            b_exp = subset(:,end);
            [~,bSG] = max(b_exp);
            [~,nbSG] = min(b_exp);
            nSG = size(subset,1);
             
            for n = 1:nSG
                % get expressin into matrix
                curr_SG = subset(n,1:end-1);
                node_exp = get_sg_matrix(nNode, curr_SG);
                if n == bSG
                    for m = 1:numel(regions)
                        cnth = cnth + 1;
                        load([top_dir, 'montages/', regions{m}, '_idx.mat'])
                        curr_exp = node_exp(idx,:);
                        p = (sum(sum(curr_exp(:,idx),2)) - sum(sum(curr_exp(:,~idx),2)))/(sum(sum(curr_exp(:,idx),2)) + sum(sum(curr_exp(:,~idx),2)));
                        high_participation{cnth,1} = subj;
                        high_participation{cnth,2} = f;
                        high_participation{cnth,3} = regions{m};
                        high_participation{cnth,4} = p;
                        high_participation{cnth,5} = betas(i);
                    end
                elseif n == nbSG
                    for m = 1:numel(regions)
                        cntl = cntl + 1;
                        load([top_dir, 'montages/', regions{m}, '_idx.mat'])
                        curr_exp = node_exp(idx,:);
                        p = (sum(sum(curr_exp(:,idx),2)) - sum(sum(curr_exp(:,~idx),2)))/(sum(sum(curr_exp(:,idx),2)) + sum(sum(curr_exp(:,~idx),2)));
                        low_participation{cntl,1} = subj;
                        low_participation{cntl,2} = f;
                        low_participation{cntl,3} = regions{m};
                        low_participation{cntl,4} = p;
                    end
                end
            end
        end
    end
    subj_h = high_participation(:,1); band_h = high_participation(:,2); region_h = high_participation(:,3);
    part_h = high_participation(:,4); slope = high_participation(:,5);
    subj_l = low_participation(:,1); band_l = low_participation(:,2); region_l = low_participation(:,3);
    part_l = low_participation(:,4);
    save([R_dir_s, 'high_participation.mat'], 'subj_h', 'band_h', 'region_h', 'part_h', 'slope')
    save([R_dir_s, 'low_participation.mat'], 'subj_l', 'band_l', 'region_l', 'part_l', 'slope')

end