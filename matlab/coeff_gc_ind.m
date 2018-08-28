%% run temporal expression statistics for gc
% independent null

% Change Log
% July 10 - updated this to run statistics for all but the highest...since
% the highest tends to look exaclty like the behavior


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

Subj = [1:20];
nSubj = numel(Subj);
freqs = [7,14;15,30;31,45;55,70];
bands = [{'alpha'}, {'beta'}];
sensors = [{'grad'}];

nNode = 102;
nEdges = (nNode^2-nNode)/2;
edge_similarity = zeros(nEdges, nSubj-1, numel(bands), numel(sensors));
% load behavior
load([top_dir, 'Behavior/perf'])

%% Loop through data
for j = 1:numel(sensors)
    sens = sensors{j};
    R_dir_s = [R_dir, sens, '/'];
    
    E_high = zeros(nSubj,numel(bands));
    E_low = zeros(nSubj,numel(bands));
    E_perf = zeros(nSubj, numel(bands));
    E_other = zeros(nSubj, numel(bands));
    E_exp_corr = [];
    H_high = zeros(nSubj-1,numel(bands));
    H_low = zeros(nSubj-1,numel(bands));
    H_perf = zeros(nSubj-1, numel(bands));
    H_other = zeros(nSubj - 1, numel(bands));
    H_exp_corr = [];
    M_high = zeros(nSubj-1,numel(bands));
    M_low = zeros(nSubj-1,numel(bands));
    M_perf = zeros(nSubj-1, numel(bands));
    M_other = zeros(nSubj - 1, numel(bands));
    M_exp_corr = [];
    S_high = zeros(nSubj-1,numel(bands));
    S_low = zeros(nSubj-1,numel(bands));
    S_perf = zeros(nSubj-1, numel(bands));
    S_other = zeros(nSubj - 1, numel(bands));
    S_exp_corr = [];
    band_order = [];
    subj_order = [];
    cnt = 0;
    for i = Subj
        s_idx = find(i == Subj);
        subj = sprintf('%03d', i);
        
        load([top_dir, 'idx_', sens, '.mat'])
        for k = 1:numel(bands)
            f = bands{k};
            
            % get subgraph data
            subset = readNPY([data_dir, subj, '/', sens, '/ind_gc_', f, '_subset.npy']);
            coeff = readNPY([data_dir, subj, '/',sens, '/ind_gc_', f, '_coeff.npy']);
            b_exp = coeff(:,end);
            [~,bSG] = max(b_exp);
            [~,nbSG] = min(b_exp);
            nSG = size(subset,1);
            
            % energy
            [ high, low, e_perf, other, exp_corr, ~, ~ ] = ...
                get_energy(subset, coeff, perf, i, k);
            E_high(s_idx,k) = high; E_low(s_idx,k) = low; E_other(s_idx,k) = other;
            E_perf(s_idx,k) = e_perf;
            E_exp_corr = [E_exp_corr; exp_corr];
            clear high low e_perf other exp_corr
            
%             % entropy
%             [ high, low, h_perf, other, exp_corr, ~, ~ ] = ...
%                 get_entropy(subset, coeff, perf, i, k);
%             H_high(s_idx,k) = high; H_low(s_idx,k) = low; H_other(s_idx,k) = other;
%             H_perf(s_idx,k) = h_perf; H_exp_corr = [H_exp_corr; exp_corr];
%             clear high low h_perf other exp_corr
            
            % mean
            [ high, low, m_perf, other, exp_corr, ~, ~ ] = ...
                get_mean(subset, coeff, perf, i, k);
            M_high(s_idx,k) = high; M_low(s_idx,k) = low; M_other(s_idx,k) = other;
            M_perf(s_idx,k) = m_perf; M_exp_corr = [M_exp_corr; exp_corr];
            
            clear high low h_perf other exp_corr
            
            %skew
            [ high, low, s_perf, other, exp_corr, b_order, s_order ] = ...
                get_skew(subset, coeff, perf, i, k);
            S_high(s_idx,k) = high; S_low(s_idx,k) = low; S_other(s_idx,k) = other;
            S_perf(s_idx,k) = s_perf; S_exp_corr = [S_exp_corr; exp_corr];
            band_order = [band_order; b_order]; subj_order = [subj_order; s_order];
            clear high low e_perf other exp_corr

        end
    end
    
    
    % plot
    % energy
    for i = 1:numel(bands)
        figure(1); clf
        boxplot([E_high(:,i), E_low(:,i), E_perf(:,i), E_other(:,i)], [{'High'}, {'Low'}, {'Behavior'}, {'Other'}])
        pause(0.01)
    end
    idx = E_exp_corr(:,2)~=0;
    E_exp_corr = E_exp_corr(idx,:);
    ebo = band_order(idx);
    eso = subj_order(idx);
    %figure(2)
    %corrplot(E_exp_corr, 'Var', {'E', 'BE'}, 'type', 'Spearman', 'testR', 'on')
    %pause(0.01)
    
    
    % save data
    save([R_dir_s, 'ind_E.mat'], 'E_high', 'E_low', 'E_other');
    save([R_dir_s, 'ind_E_corr.mat'], 'E_exp_corr', 'ebo', 'eso');
    
    
%     % entropy
%     for n = 1:numel(bands)
%         figure(1); clf
%         boxplot([H_high(:,n), H_low(:,n), H_perf(:,n), H_other(:,n)], [{'High'}, {'Low'}, {'Behavior'}, {'Other'}])
%         pause(0.01)
%     end
%     figure(2)
%     H_exp_corr = H_exp_corr(idx,:);
%     hbo = band_order(idx);
%     hso = subj_order(idx);
%     corrplot(H_exp_corr, 'Var', {'H', 'BE'}, 'type', 'Spearman', 'testR', 'on')
%     pause(0.01)
%     
%     % save data
%     save([R_dir_s, 'ind_H.mat'], 'H_high', 'H_low', 'H_other');
%     save([R_dir_s, 'ind_H_corr.mat'], 'H_exp_corr', 'hbo', 'hso');
    
    % mean
    for n = 1:numel(bands)
        figure(1); clf
        boxplot([M_high(:,n), M_low(:,n), M_other(:,n)], [{'High'}, {'Low'},  {'Other'}])
        pause(0.01)
    end
    figure(2)
    M_exp_corr = M_exp_corr(idx,:);
    mbo = band_order(idx);
    mso = subj_order(idx);
    %corrplot(M_exp_corr, 'Var', {'H', 'BE'}, 'type', 'Spearman', 'testR', 'on')
    %pause(0.01)
    
    % save data
    save([R_dir_s, 'ind_M.mat'], 'M_high', 'M_low', 'M_other');
    save([R_dir_s, 'ind_M_corr.mat'], 'M_exp_corr', 'mbo', 'mso');
    
    
    
    % skew
    for n = 1:numel(bands)
        figure(1); clf
        boxplot([S_high(:,n), S_low(:,n), S_other(:,n)], [{'High'}, {'Low'}, {'Other'}])
        pause(0.01)
    end
    S_exp_corr = S_exp_corr(idx,:);
    sbo = band_order(idx);
    sso = subj_order(idx);
    %figure(2)
    %corrplot(S_exp_corr, 'Var', {'S', 'BE'}, 'type', 'Spearman', 'testR', 'on')
    %pause(0.01)
    % save data
    save([R_dir_s, 'ind_S.mat'], 'S_high', 'S_low', 'S_other');
    save([R_dir_s, 'ind_S_corr.mat'], 'S_exp_corr', 'sbo', 'sso');
    
end
% %% Repeat for Baseline
% Subj = [2,4:16];
% nSubj = numel(Subj);
% % Loop through data
% E_high = zeros(nSubj,numel(bands));
% E_low = zeros(nSubj,numel(bands));
% E_perf = zeros(nSubj, numel(bands));
% E_other = zeros(nSubj, numel(bands));
% E_exp_corr = [];
% H_high = zeros(nSubj-1,numel(bands));
% H_low = zeros(nSubj-1,numel(bands));
% H_perf = zeros(nSubj-1, numel(bands));
% H_other = zeros(nSubj - 1, numel(bands));
% H_exp_corr = [];
% S_high = zeros(nSubj-1,numel(bands));
% S_low = zeros(nSubj-1,numel(bands));
% S_perf = zeros(nSubj-1, numel(bands));
% S_other = zeros(nSubj - 1, numel(bands));
% S_exp_corr = [];
% high_corr = [];
% band_order = [];
% subj_order = [];
% cnt = 0;
% for i = Subj
%     subj = sprintf('%03d', i);
%     for j = 1:numel(sensors)
%         sens = sensors{j};
%         load([top_dir, 'idx_', sens, '.mat'])
%         for k = 1:numel(bands)
%             f = bands{k};
%
%             % get subgraph data
%             subset = readNPY([data_dir, subj, '/', sens, '/ac_', f, '_subset_bl.npy']);
%             coeff = readNPY([data_dir, subj, '/',sens, '/ac_', f, '_coeff_bl.npy']);
%             b_exp = coeff(:,end);
%             [~,bSG] = max(b_exp);
%             [~,nbSG] = min(b_exp);
%             nSG = size(subset,1);
%
%             % energy
%             [ high, low, e_perf, other, exp_corr, b_order, s_order ] = ...
%                 get_energy(coeff, subset, perf, i, k);
%             E_high(i,k) = high; E_low(i,k) = low; E_other(i,k) = other;
%             E_perf(i,k) = e_perf;
%             band_order = [band_order; b_order];
%             subj_order = [subj_order; s_order];
%             high_corr = [high_corr, low, coeff(nbSG, end)];
%             clear high low e_perf other
%
%             % entropy
%             [ high, low, h_perf, other, exp_corr, b_order, s_order ] = ...
%                 get_energy(coeff, subset, perf, i, k);
%             H_high(i,k) = high; H_low(i,k) = low; H_other(i,k) = other;
%             H_perf(i,k) = h_perf;
%             clear high low h_perf other
%
%             %skew
%             [ high, low, s_perf, other, exp_corr, b_order, s_order ] = ...
%                 get_energy(coeff, subset, perf, i, k);
%             S_high(i,k) = high; S_low(i,k) = low; S_other(i,k) = other;
%             S_perf(i,k) = s_perf;
%             clear high low e_perf other
%         end
%     end
% end
%
% % plot
% % energy
% for i = 1:numel(bands)
%     figure(1); clf
%     boxplot([E_high(:,i), E_low(:,i), E_perf(:,i), E_other(:,i)], [{'High'}, {'Low'}, {'Behavior'}, {'Other'}])
%     saveas(gca, [save_dir, bands{i}, '_energy_bl.png'], 'png')
% end
% %corrplot(E_exp_corr, 'Var', {'E', 'BE'}, 'type', 'Spearman', 'testR', 'on')
% %saveas(gca, [save_dir, 'E_BE_corr_all_bl.png'], 'png')
% %corrplot(high_corr, 'Var', {'E', 'BE'}, 'type', 'Spearman', 'testR', 'on')
% %saveas(gca, [save_dir, 'E_BE_corr_low_bl.png'], 'png')
% % save data
% save([R_dir, 'E_bl.mat'], 'E_high', 'E_low', 'E_other');
% save([R_dir, 'E_corr_bl.mat'], 'E_exp_corr', 'band_order', 'subj_order');
%
%
% % entropy
% for i = 1:numel(bands)
%     figure(1); clf
%     boxplot([H_high(:,i), H_low(:,i), H_perf(:,i), H_other(:,i)], [{'High'}, {'Low'}, {'Behavior'}, {'Other'}])
%     saveas(gca, [save_dir, bands{i}, '_entropy_bl.png'], 'png')
% end
% % save data
% save([R_dir, 'H_bl.mat'], 'H_high', 'H_low', 'H_other');
% save([R_dir, 'H_corr_bl.mat'], 'H_exp_corr', 'band_order', 'subj_order');
%
%
%
% % skew
% for i = 1:numel(bands)
%     figure(1); clf
%     boxplot([S_high(:,i), S_low(:,i), S_perf(:,i), S_other(:,i)], [{'High'}, {'Low'}, {'Behavior'}, {'Other'}])
%     saveas(gca, [save_dir, bands{i}, '_skew_bl.png'], 'png')
% end
% % save data
% save([R_dir, 'S_bl.mat'], 'S_high', 'S_low', 'S_other');
% save([R_dir, 'S_corr_bl.mat'], 'S_exp_corr', 'band_order', 'subj_order');

