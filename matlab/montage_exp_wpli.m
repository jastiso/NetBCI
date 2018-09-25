%% Regional expression

addpath(genpath('/Users/stiso/Documents/MATLAB/npy-matlab-master/'))
addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/')
R_dir = '/Users/stiso/Documents/R/NetBCI/data/wpli/';

top_dir = '/Users/stiso/Documents/MATLAB/NetBCI/';
data_dir = '/Users/stiso/Documents/Python/NetBCI/NMF/';
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

Subj = [1:20];
nSubj = numel(Subj);
freqs = [7,14;15,30;31,45;55,70];
bands = [{'alpha'}, {'beta'}, {'low_gamma'}, {'gamma'}];
sensors = [{'grad'}];
regions = [{'Left_frontal'}, {'Left_occipital'}, {'Left_parietal'}, {'Left_temporal'}, ...
    {'Right_frontal'}, {'Right_occipital'}, {'Right_parietal'}, {'Right_temporal'}, {'Vertex'}];

nNode = 102;
nEdges = (nNode^2-nNode)/2;
edge_similarity = zeros(nEdges, nSubj-1, numel(bands), numel(sensors));
% initialize
mod_exp_high = cell(nSubj*numel(bands)*numel(regions), 4);
mod_exp_low = cell(nSubj*numel(bands)*numel(regions), 4);
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
            subset = readNPY([data_dir, subj, '/', sens, '/wpli_', f, '_subset.npy']);
            coeff = readNPY([data_dir, subj, '/',sens, '/wpli_', f, '_coeff.npy']);
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
                        curr_exp = sum(sum(node_exp(idx,idx)))/numel(nonzeros(idx))^2;
                        mod_exp_high{cnth,1} = subj;
                        mod_exp_high{cnth,2} = f;
                        mod_exp_high{cnth,3} = regions{m};
                        mod_exp_high{cnth,4} = curr_exp;
                        mod_exp_high{cnth,5} = betas(i);

                    end
                elseif n == nbSG
                    for m = 1:numel(regions)
                        cntl = cntl + 1;
                        load([top_dir, 'montages/', regions{m}, '_idx.mat'])
                        curr_exp = sum(sum(node_exp(idx,idx)))/numel(nonzeros(idx))^2;
                        mod_exp_low{cntl,1} = subj;
                        mod_exp_low{cntl,2} = f;
                        mod_exp_low{cntl,3} = regions{m};
                        mod_exp_low{cntl,4} = curr_exp;
                    end
                end
            end

        end
    end
    subj_h = mod_exp_high(:,1); band_h = mod_exp_high(:,2); region_h = mod_exp_high(:,3);
    exp_h = mod_exp_high(:,4); slope = mod_exp_high(:,5);
    subj_l = mod_exp_low(:,1); band_l = mod_exp_low(:,2); region_l = mod_exp_low(:,3);
    exp_l = mod_exp_low(:,4);
    save([R_dir_s, 'mod_exp_high.mat'], 'subj_h', 'band_h', 'region_h', 'exp_h', 'slope')
    save([R_dir_s, 'mod_exp_low.mat'], 'subj_l', 'band_l', 'region_l', 'exp_l', 'slope')

end
    