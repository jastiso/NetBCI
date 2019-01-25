%% Compare Consistency

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

load([save_dir, 'consistent_sums_pr.mat'])
sum_high_pr = sum_high;
sum_high2_pr = sum_high2;
sum_high3_pr = sum_high3;
sum_low_pr = sum_low;

load([save_dir, 'consistent_sums.mat'])
sum_high = sum_high(3,:);
sum_high2 = sum_high2(3,:);
sum_high3 = sum_high3(3,:);
sum_low = sum_low(3,:);

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
