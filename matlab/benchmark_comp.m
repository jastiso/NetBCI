%% Network Comparison

% @author JStiso jeni.stiso@gmail.com


addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830')
addpath(genpath('/Users/stiso/Documents/MATLAB/BCT'))

top_dir = '/Users/stiso/Documents/MATLAB/NetBCI/';
data_dir = [top_dir, 'DataBase/1_Signals/2_Segmentation/2_MEG/'];
raw_dir = [top_dir, 'DataBase/1_Signals/0_RawData/2_MEG/'];
sessions = [{'Session1'}, {'Session2'}, {'Session3'}, {'Session4'}];
condition = [{'test01'}, {'test02'}, {'test03'}, {'test04'}, {'test05'}, {'test06'}, {'RS1'}, {'RS2'}];
nSubj = 15;
nShuffle = 100;
thresholds = 0%:.1:.9;

emp_path_length = zeros(numel(sessions), numel(condition), nSubj, numel(thresholds));
emp_cluster_coef = zeros(numel(sessions), numel(condition), nSubj, numel(thresholds));
rand_path_length = zeros(numel(sessions), numel(condition), nSubj, numel(thresholds), nShuffle);
rand_cluster_coef = zeros(numel(sessions), numel(condition), nSubj, numel(thresholds), nShuffle);
reg_path_length = zeros(numel(sessions), numel(condition), nSubj, numel(thresholds), 1);
reg_cluster_coef = zeros(numel(sessions), numel(condition), nSubj, numel(thresholds), 1);
emp_skewness = zeros(numel(sessions), numel(condition), nSubj);
rand_skewness = zeros(numel(sessions), numel(condition), nSubj,  nShuffle);
reg_skewness = zeros(numel(sessions), numel(condition), nSubj);
emp_G = zeros(numel(sessions), numel(condition), nSubj);
rand_G = zeros(numel(sessions), numel(condition), nSubj,  nShuffle);
reg_G = zeros(numel(sessions), numel(condition), nSubj);
%% Start Loop

for i = 1:numel(sessions)
    session = sessions{i};
    for j = 1:numel(condition)
        cond = condition{j};
        for k = 1:nSubj
            ext = sprintf('%03d',k);
            img_dir = [top_dir, session, '/', cond, '/', ext, '/diagnostics/'];
            save_dir = [top_dir, session, '/', cond, '/', ext, '/FCmatrices/'];
            if exist([save_dir, 'beta_imag_coh.mat'], 'file') ~= 0
                load([save_dir, 'beta_imag_coh.mat'])
                A = coherence.cohspctrm;
                A(A<0) = 0;
                [sort_weight, idx] = sort(abs(reshape(A, [],1)),'descend');
                n_zero = sum(sum(A==0));
                
                
                % collect global efficiency, and path length for all graphs at
                % all thresholds
                for t = 1:numel(thresholds)
                    fprintf('\n%d', t)
                    thresh = thresholds(t);
                    connectivity= A;
                    % make comparison networks (regular and random)
                    cut_off = round(numel(idx)*thresh);
                    connectivity(idx(1:cut_off)) = 0;
                    figure(1); clf
                    imagesc(connectivity); colorbar
                    pause(.01)
                    
                    N = size(connectivity, 1);
                    w = nonzeros(connectivity); % weights
                    p = numel(w);
                    nswaps = 2^4;
                    
                    % rand
                    rand_graph_wei = zeros(N,N,nShuffle);
                    for m = 1:nShuffle
                        rand_graph = randmio_dir(connectivity,nswaps);
                        rand_graph_wei(:,:,m) = rand_graph;
                        %rand_graph_wei(:,:,j) = fcn_match_strength(rand_graph,w,stren,temp,dec,energyfcn,tolerance);
                    end
                    reg_graph_wei = zeros(N,N,1);
                    for m = 1:1
                        reg_graph = makeringlatticeCIJ(N,p);
                        reg_graph(reg_graph ~= 0) = w(randperm(p));
                        reg_graph_wei(:,:,m) = reg_graph;
                        %rand_graph_wei(:,:,j) = fcn_match_strength(rand_graph,w,stren,temp,dec,energyfcn,tolerance);
                    end
                    
                    if thresh == 0
                        % compare degree distributions
                        [N1, X] = hist(sum(connectivity), 20);
                        [N2, X1] = hist(reshape(sum(reg_graph_wei), N, []), 20);
                        [N3, X2] = hist(reshape(sum(rand_graph_wei), N, []), 20);
                        figure('position', [1, 1, 1000, 300]); clf;
                        subplot(1,3,1);
                        bar(X, N1, 'facecolor', 'r', 'edgecolor', 'r','facealpha', .5); hold on
                        title('Empirical Degree Distribution'); xlabel('Degree'); ylabel('Number of Nodes')
                        subplot(1,3,3)
                        bar(X2, N3, 'facecolor', 'k', 'edgecolor', 'k','facealpha', .5);
                        title('Random Degree Distribution'); xlabel('Degree'); ylabel('Number of Nodes')
                        tmp_lim = get(gca, 'xlim');
                        subplot(1,3,2)
                        bar(X1, N2, 'facecolor', 'b', 'edgecolor', 'b','facealpha', .5);
                        title('Regular Degree Distribution'); xlabel('Degree'); ylabel('Number of Nodes')
                        xlim(tmp_lim)
                        saveas(gca, [img_dir, session, cond, '_', num2str(k),  '_dd.png'], 'png');
                        emp_skewness(i,j,k) = skewness(sum(connectivity));
                        rand_skewness(i,j,k,:) = skewness(reshape(sum(rand_graph_wei), N, []));
                        reg_skewness(i,j,k,:) = skewness(reshape(sum(reg_graph_wei), N, []));
                        
                        D = diag(sum(connectivity));
                        G = expm(D^-.5*connectivity*D^-.5);
                        D = diag(sum(reg_graph_wei));
                        G_reg = expm(D^-.5*reg_graph_wei*D^-.5);
                        G_rand = zeros(N,N,nShuffle);
                        for m = 1:nShuffle
                            D = diag(sum(rand_graph_wei(:,:,m)));
                            G_rand(:,:,m) = expm(D^-.5*rand_graph_wei(:,:,m)*D^-.5);
                        end
                        [N1, X] = hist(sum(G), 20);
                        [N2, X1] = hist(sum(G_reg), 20);
                        [N3, X2] = hist(reshape(sum(G_rand), N, []), 20);
                        figure('position', [1, 1, 1000, 300]); clf;
                        subplot(1,3,1);
                        bar(X, N1, 'facecolor', 'r', 'edgecolor', 'r','facealpha', .5); hold on
                        title('Empirical G Distribution'); xlabel('G'); ylabel('Number of Nodes')
                        subplot(1,3,3)
                        bar(X2, N3, 'facecolor', 'k', 'edgecolor', 'k','facealpha', .5);
                        title('Random G Distribution'); xlabel('G'); ylabel('Number of Nodes')
                        tmp_lim = get(gca, 'xlim');
                        subplot(1,3,2)
                        bar(X1, N2, 'facecolor', 'b', 'edgecolor', 'b','facealpha', .5);
                        title('Regular G Distribution'); xlabel('G'); ylabel('Number of Nodes')
                        xlim(tmp_lim)
                        saveas(gca, [img_dir, session, cond, '_', num2str(k),  '_G.png'], 'png');
                        emp_G(i,j,k) = mean(sum(G));
                        rand_G(i,j,k,:) = mean(reshape(sum(rand_graph_wei), N, []));
                        reg_G(i,j,k,:) = mean(reshape(sum(reg_graph_wei), N, []));
                        close all
                    end
                    %reg_graph_wei = fcn_match_strength(reg_graph,w,stren,temp,dec,energyfcn,tolerance);
                    
                    % plot three matrices
                    % get coordinates
                    tmp_reg = makeringlatticeCIJ(50,numel(nonzeros(connectivity(1:20, 1:20))));
                    %tmp_reg(tmp_reg ~= 0) = nonzeros(bin_con(1:0, 1:20));
                    xCenter = 12;
                    yCenter = 10;
                    theta = linspace(0,2*pi,50);
                    radius = 5;
                    x = radius * cos(theta) + xCenter;
                    y = radius * sin(theta) + yCenter;
                    coord = [x;y]';
                    figure('position', [1, 1, 1500, 300])
                    subplot(1,3,1)
                    gplot(connectivity(1:20, 1:20), coord, '-or'); title('Empirical Network')
                    subplot(1,3,2)
                    gplot(tmp_reg, coord, '-ob'); title('Regular Graph')
                    subplot(1,3,3)
                    gplot(rand_graph_wei(1:20, 1:20, 1), coord, '-ok'); title('Random Graph')
                    %saveas(gca, [img_dir, res, '_graph_ex.png'], 'png');
                    close
                    
                    % get global efficiency and path length
                    D = diag(sum(connectivity));
                    G = expm(D^-.5*connectivity*D^-.5);
                    emp_path_length(i,j,k,t) = mean(sum(G));
                    emp_cluster_coef(i,j,k,t) = mean(clustering_coef_wd(connectivity));
                    % for random and reg
                    for m = 1:nShuffle
                        D = diag(sum(rand_graph_wei(:,:,m)));
                        G = expm(D^-.5*rand_graph_wei(:,:,m)*D^-.5);
                        rand_path_length(i,j,k,t,m) = mean(sum(G));
                        rand_cluster_coef(i,j,k,t,m) = mean(clustering_coef_wd(rand_graph_wei(:,:,m)));
                        
                    end
                    D = diag(sum(reg_graph_wei(:,:,1)));
                    G = expm(D^-.5*reg_graph_wei(:,:,1)*D^-.5);
                    reg_path_length(i,j,k,t) = mean(sum(G));
                    reg_cluster_coef(i,j,k,t,1) = mean(clustering_coef_wd(reg_graph_wei(:,:,1)));
                end
                % plot
%                 avg_rand_cc = squeeze(mean(rand_cluster_coef(i,j,k,:,:),5));
%                 sd_rand_cc = std(squeeze(rand_cluster_coef(i,j,k,:,:)),0,2);
%                 avg_reg_cc = squeeze(mean(reg_cluster_coef(i,j,k,:,:),5));
%                 sd_reg_cc = std(squeeze(reg_cluster_coef(i,j,k,:,:)),0,2);
%                 %cc
%                 figure(1), clf;
%                 plot(thresholds, squeeze(emp_cluster_coef(i,j,k,:)), 'k', 'linewidth', 2); hold on
%                 plot(thresholds, avg_rand_cc, 'r', 'linewidth', 2);
%                 plot(thresholds, avg_reg_cc, 'b', 'linewidth', 2);
%                 %shade_plot(thresholds', avg_rand_cc, sd_rand_cc, rgb('Salmon'))
%                 title([session, cond, ' ', num2str(k), ' CC']); xlabel('Threshold'); ylabel('Clustering Coefficient'); hold off
%                 saveas(gca,[img_dir, session, cond, '_', num2str(k), 'benchmark_cc.png'], 'png');
%                 
%                 avg_rand_pl = squeeze(mean(rand_path_length(i,j,k,:,:),5));
%                 sd_rand_pl = std(squeeze(rand_path_length(i,j,k,:,:)),0,2);
%                 avg_reg_pl = squeeze(mean(reg_path_length(i,j,k,:,:),5));
%                 sd_reg_pl = std(squeeze(reg_path_length(i,j,k,:,:)),0,2);
%                 %cc
%                 figure(2)
%                 plot(thresholds, log10(squeeze(emp_path_length(i,j,k,:))), 'k', 'linewidth', 2); hold on
%                 plot(thresholds, log10(avg_rand_pl), 'r', 'linewidth', 2);
%                 plot(thresholds, log10(avg_reg_pl), 'b', 'linewidth', 2);
%                 %shade_plot(thresholds', avg_rand_pl, sd_rand_pl, rgb('Salmon'))
%                 %shade_plot(thresholds', avg_reg_pl, sd_reg_pl, rgb('DeepSkyBlue'))
%                 title([session, cond, ' ', num2str(k), ' PL']); xlabel('Threshold'); ylabel('Path Length'); hold off
%                 saveas(gca,[img_dir, session, cond, '_', num2str(k), 'benchmark_pl.png'], 'png');
            end
        end
    end
end
save([top_dir, 'emp_cc'], 'emp_cluster_coef');
save([top_dir, 'emp_pl'], 'emp_path_length');
save([top_dir, 'rand_cc'], 'rand_cluster_coef');
save([top_dir, 'rand_pl'], 'rand_path_legth');
save([top_dir, 'reg_cc'], 'reg_cluster_coef');
save([top_dir, 'reg_pl'], 'reg_path_legth');
save([top_dir, 'emp_deg_skew'], 'emp_skewness');
save([top_dir, 'rand_deg_skew'], 'rand_skewness');
save([top_dir, 'reg_deg_skew'], 'reg_skewness');
save([top_dir, 'emp_G_skew'], 'emp_G');
save([top_dir, 'rand_G_skew'], 'rand_G');
save([top_dir, 'reg_G_skew'], 'reg_G');
emp_skewness(emp_skewness == 0) = [];
rand_skewness(rand_skewness == 0) = [];
reg_skewness(reg_skewness == 0) = [];
emp_G(emp_G == 0) = [];
rand_G(rand_G == 0) = [];
reg_G(reg_G == 0) = [];

[h,p,ci, stat] = ttest2(emp_skewness, rand_skewness)
[h,p,ci,stat] = ttest2(emp_skewness, reg_skewness)
[h,p,ci,stat] = ttest2(emp_G, rand_G)
[h,p,ci,stat] = ttest2(emp_G, reg_G)