%% Preprocessing

% preprocess MEG data, and get into connectivity matrix
% @author JStiso jeni.stio@gmail.com Nov 2017

% CHange Log
% August 27, 2018 - updated to do wpli, which does not include a baseline
% correction



%% Define Global Variables

addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830')
addpath(genpath('/Users/stiso/Documents/MATLAB/eeglab_current'))
addpath(genpath('/Users/stiso/Documents/Code/NetBCI/matlab/'))

top_dir = '/Users/stiso/Documents/MATLAB/NetBCI/';
data_dir = '/Users/stiso/Resilio Sync/NETBCI.RAW/DataBase/1_Signals/2_Segmentation/2_MEG/';
raw_dir = '/Users/stiso/Resilio Sync/NETBCI.RAW/DataBase/1_Signals/0_RawData/2_MEG/';
sessions = [{'Session1'}, {'Session2'}, {'Session3'}, {'Session4'}];
condition = [{'test01'}, {'test02'}, {'test03'}, {'test04'}, {'test05'}, {'test06'}]; % not including rest
subjs = [1:20];
nNode = 102;

%freqs = [7,14;15,30;31,45;55,70];
%bands = [{'alpha'},{'beta'},{'low_gamma'}, {'gamma'}];
freqs = [8,10];
bands = [{'mu'}];

st = 3;
en = 6; % in seconds, the feedback period: 3-6s

bl_st = 0.2; % 0-1 is baseline
bl_en = 0.8;
num_neg_grad = [];
num_neg_mag = [];
cnte = 1;
% GC params
dt = 100;
lag = 10;
t0 = dt+lag;

%% Start Loop

errors = cell(1,4);

parfor i = subjs
    subj = i;
   [ errors{i} ] = wrapper_wpli( sessions, condition, subj, data_dir, raw_dir, top_dir, bands, freqs, st, en, nNode );
end




