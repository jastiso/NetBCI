%% Get correlation between expression and behavior
% mostly a sonity check

addpath(genpath('/Users/stiso/Documents/MATLAB/npy-matlab-master/'))
addpath(genpath('/Users/stiso/Documents/MATLAB/mult_comp_perm_corr/'))
addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/')
R_dir_o = '/Users/stiso/Documents/R/NetBCI/data/gc/';

top_dir = '/Users/stiso/Documents/MATLAB/NetBCI/';
data_dir = '/Users/stiso/Documents/Python/NetBCI/NMF/gc/';
save_dir = [top_dir, 'GroupAvg/analysis/'];
img_dir = [top_dir, 'GroupAvg/images/'];
% make directories
if ~exist(img_dir, 'dir')
    mkdir(img_dir)
end
% make directories
if ~exist(save_dir, 'dir')
    mkdir(save_dir)
end

Subj = [1:20];
nSubj = numel(Subj);
freqs = [7,14;15,30;31,45;55,70];
bands = [{'alpha'},{'beta'},{'low_gamma'},{'gamma'}];

sensors = [{'grad'}];
% load behavior
load([top_dir, 'Behavior/behavior_all'])
load([top_dir, 'Behavior/stats'])
%initialize
max_corrs = zeros(nSubj, numel(bands), numel(sensors));
min_corrs = zeros(nSubj, numel(bands), numel(sensors));
max_exp = zeros(nSubj, numel(bands), numel(sensors));
min_exp = zeros(nSubj, numel(bands), numel(sensors));
sum_exp = zeros(nSubj, numel(bands), numel(sensors));
sd_exp = zeros(nSubj, numel(bands), numel(sensors));


%% Loop through data

for j = 1:numel(sensors)
    sens = sensors{j};
    R_dir = [R_dir_o, sens, '/'];
    
    for i = Subj
        idx = find(i == Subj);
        subj = sprintf('%03d', i);
        
        for k = 1:numel(bands)
            f = bands{k};
            
            % get subgraph data
            subset = readNPY([data_dir, subj,  '/', sens, '/gc_', f, '_subset.npy']);
            coeff = readNPY([data_dir, subj,  '/', sens, '/gc_', f, '_coeff.npy']);
            
            %get subgraph with highest expression
            [max_e,bsg] = max(subset(:,end));
            [min_e,nbsg] = min(subset(:,end));
            
            % get corr
            max_corrs(idx,k,j) = corr(behavior_all{i}', coeff(bsg,:)');
            min_corrs(idx,k,j) = corr(behavior_all{i}', coeff(nbsg,:)');
            max_exp(idx,k,j) = max_e;
            min_exp(idx,k,j) = min_e;
            sum_exp(idx,k,j) = mean(subset(:,end));
            sd_exp(idx,k,j) = std(subset(:,end));
            figure(1); clf
            hist(subset(:,end))
            pause(0.001)
            saveas(gca, [save_dir, subj, '_', f, '_exp_hist.png'], 'png')
        end
    end
    
    
    % make picture
    figure(1)
    imagesc(max_corrs);
    caxis([min(min(min_corrs)), max(max(max_corrs))])
    colorbar
    xlabel('Frequency')
    ylabel('Subject')
    saveas(gca, [save_dir, 'behavior_exp_corr.png'], 'png')
    save([R_dir, 'max_exp_b_corr.mat'], 'max_corrs')
    save([R_dir, 'min_exp_b_corr.mat'], 'min_corrs')
    
    figure(2)
    imagesc(min_corrs);
    caxis([min(min(min_corrs)), max(max(max_corrs))])
    colorbar
    xlabel('Frequency')
    ylabel('Subject')
    saveas(gca, [save_dir, 'min_behavior_exp_corr.png'], 'png')
    
    % plot realtionship between slope of learning and expression
    for i = 1:numel(bands)
        corrplot([betas(Subj)', improvement(Subj)', max_exp(:,i), sum_exp(:, i), sd_exp(:,i)], 'Var', {'Beta', 'Im', 'Max', 'mean', 'sd'}, 'type', 'pearson','testr', 'on')
        pause(0.001)
    end
    
    figure(1); clf
    corrplot([repmat(betas(Subj)',numel(bands),1), repmat(improvement(Subj)',numel(bands),1),  ...
        reshape(max_exp,[],1),reshape(sum_exp(:, :), [],1), reshape(sd_exp(:, :), [],1)], 'Var',...
        {'Beta',  'Im', 'Max', 'mean', 'sd'}, 'type', 'spearman','testr', 'on')
    
    % save to R_dir
    save([R_dir, 'exp_beh_cor.mat'], 'betas',  'max_exp', 'sum_exp', 'sd_exp');
    
end



% select individual corrplot
figure(1); clf
s1 = plot(betas(Subj)', max_exp(:,1), 'k.', 'color', [.7,.7,.7]);
set(s1, 'MarkerSize', 20, 'LineWidth', 2);
%%% regression line
hold on
l = lsline ;
set(l,'color', 'k', 'LineWidth', 2)
hold off
[p,r] = mult_comp_perm_corr(betas(Subj)', max_exp(:,1))
saveas(gca, [save_dir, 'max_beta_alpha.png'], 'png')

figure(1); clf
s1 = plot(betas(Subj)', sd_exp(:,1), 'k.', 'color', [.7,.7,.7]);
set(s1, 'MarkerSize', 20, 'LineWidth', 2);
%%% regression line
hold on
l = lsline ;
set(l,'color', 'k', 'LineWidth', 2)
hold off
[p,r] = mult_comp_perm_corr(betas(Subj)', sd_exp(:,1))
saveas(gca, [save_dir, 'sd_beta_alpha.png'], 'png')

figure(1); clf
s1 = plot(betas(Subj)', max_exp(:,2), 'k.', 'color', [.7,.7,.7]);
set(s1, 'MarkerSize', 20, 'LineWidth', 2);
%%% regression line
hold on
l = lsline ;
set(l,'color', 'k', 'LineWidth', 2)
hold off
[p,r] = mult_comp_perm_corr(betas(Subj)', max_exp(:,2))
saveas(gca, [save_dir, 'max_beta_beta.png'], 'png')

figure(1); clf
s1 = plot(betas(Subj)', sd_exp(:,2), 'k.', 'color', [.7,.7,.7]);
set(s1, 'MarkerSize', 20, 'LineWidth', 2);
%%% regression line
hold on
l = lsline ;
set(l,'color', 'k', 'LineWidth', 2)
hold off
[p,r] = mult_comp_perm_corr(betas(Subj)', sd_exp(:,2))
saveas(gca, [save_dir, 'sd_beta_beta.png'], 'png')





% %% Now for BL
%
% % Loop through data
%
% for i = 2:5
%     subj = sprintf('%03d', i);
%     for j = 1:numel(sensors)
%         sens = sensors{j};
%         for k = 1:numel(bands)
%             f = bands{k};
%
%             % get subgraph data
%             subset = readNPY([data_dir, subj, '/', f, '_subset_bl.npy']);
%             coeff = readNPY([data_dir, subj, '/', f, '_coeff_bl.npy']);
%
%             %get subgraph with highest expression
%             [max_e,bsg] = max(coeff(:,end));
%             [min_e,nbsg] = min(coeff(:,end));
%
%             % get corr
%             max_corrs(i,k,j) = corr(perf(:,i), subset(bsg,:)');
%             min_corrs(i,k,j) = corr(perf(:,i), subset(nbsg,:)');
%             max_exp(i,k,j) = max_e;
%             min_exp(i,k,j) = min_e;
%             sum_exp(i,k,j) = mean(coeff(:,end));
%             sd_exp(i,k,j) = std(coeff(:,end));
%             hist(coeff(:,end))
%             pause(0.001)
%         end
%     end
% end
%
% % make picture
% figure(1)
% imagesc(max_corrs);
% caxis([min(min(min_corrs)), max(max(max_corrs))])
% colorbar
% xlabel('Frequency')
% ylabel('Subject')
% saveas(gca, [save_dir, 'behavior_exp_corr_bl.png'], 'png')
% save([R_dir, 'max_exp_b_corr_bl.mat'], 'max_corrs')
% save([R_dir, 'min_exp_b_corr_bl.mat'], 'min_corrs')
%
% figure(2)
% imagesc(min_corrs);
% caxis([min(min(min_corrs)), max(max(max_corrs))])
% colorbar
% xlabel('Frequency')
% ylabel('Subject')
% saveas(gca, [save_dir, 'min_behavior_exp_corr_bl.png'], 'png')
%
% % plot realtionship between slope of learning and expression
% for i = 1:numel(bands)
%     figure(i); clf
%     corrplot([betas(2:7)', improvement(2:7)', max_exp(2:end,i), min_exp(2:end,i), sum_exp(2:end, i), sd_exp(2:end,i)], 'Var', {'Beta', 'Im', 'Max', 'Min', 'mean', 'sd'}, 'type', 'spearman','testr', 'on')
%     pause(0.001)
% end
%
% figure(1); clf
% corrplot([repmat(betas(2:7)',4,1), repmat(improvement(2:7)',4,1),  ...
%     reshape(max_exp(2:end,:),[],1), reshape(min_exp(2:end,:),[],1), reshape(sum_exp(2:end, :), [],1), reshape(sd_exp(2:end, :), [],1)], 'Var',...
%     {'Beta',  'Im', 'Max', 'Min', 'mean', 'sd'}, 'type', 'spearman','testr', 'on')
%
