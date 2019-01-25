%% Get lobes of most consistent edges

% Global Variables

addpath(genpath('/Users/stiso/Documents/MATLAB/npy-matlab-master/'))
addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/')
addpath(genpath('/Users/stiso/Documents/MATLAB/BCT/'))

top_dir = '/Users/stiso/Documents/MATLAB/NetBCI/';
data_dir = '/Users/stiso/Documents/Python/NetBCI/NMF/param/';
save_dir = [top_dir, 'GroupAvg/wpli/analysis/'];
img_dir = [top_dir, 'GroupAvg/wpli/images/'];

subjs = [1:20];
nSubj = numel(subjs);
nNode = 102;
nEdge = (nNode^2-nNode)/2;
freqs = [7,14;15,30;31,45;55,70];
bands = [{'alpha'},{'beta'},{'low_gamma'}];
sensors = [{'grad'}];

load([save_dir, 'most_consistent_pr.mat'])
high_pr = c_high_max;
high2_pr = c_high2_max;
high3_pr = c_high3_max;
low_pr = c_low_max;

load([save_dir, 'most_consistent.mat'])
high = c_high_max(:,3);
high2 = c_high2_max(:,3);
high3 = c_high3_max(:,3);
low = c_low_max(:,3);

all_cond = [sum_high; sum_high2; sum_high3; sum_low];
all_cond_pr = [sum_high_pr; sum_high2_pr; sum_high3_pr; sum_low_pr];
cond_diff = all_cond - all_cond_pr;
imagesc(cond_diff); colorbar

data = [sum_high, sum_high2, sum_high3, sum_low, sum_high_pr, sum_high2_pr, sum_high3_pr, sum_low_pr];
cond = [repmat({'high'},1,3), repmat({'high2'},1,3), repmat({'high3'},1,3), repmat({'low'},1,3), repmat({'high'},1,3), repmat({'high2'},1,3), repmat({'high3'},1,3), repmat({'low'},1,3)];
band = repmat([{'alpha'}, {'beta'}, {'low_gamma'}], 1,8);
model = [repmat({'emp'}, 1,12), repmat({'upr'}, 1,12)];

clear X
band = categorical(band'); model = categorical(model'); cond = categorical(cond');
X = table(data',band,model,cond);
fit = fitlm(X, 'Var1~band+model+cond')
anova(fit)
