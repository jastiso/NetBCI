%% Preprocessing Bio Channels

% preprocess MEG data and Bio data, gets connectivity betweeen bio and MEG
% @author JStiso jeni.stio@gmail.com Feb 2020





%% Define Global Variables

addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830')
addpath(genpath('/Users/stiso/Documents/MATLAB/eeglab_current'))
addpath(genpath('/Users/stiso/Documents/Code/netBCI/matlab/'))

top_dir = '/Users/stiso/Documents/MATLAB/NetBCI/';
data_dir = '/Users/stiso/Resilio Sync/NETBCI.RAW/DataBase/1_Signals/2_Segmentation/2_MEG/';
raw_dir = '/Users/stiso/Resilio Sync/NETBCI.RAW/DataBase/1_Signals/0_RawData/2_MEG/';
bio_dir = '/Users/stiso/Resilio Sync/NETBCI.RAW/DataBase/1_Signals/2_Segmentation/0_BIO/';

sessions = [{'Session1'}, {'Session2'}, {'Session3'}, {'Session4'}];
condition = [{'test01'}, {'test02'}, {'test03'}, {'test04'}, {'test05'}, {'test06'}]; % not including rest
subjs = [17:20];
nNode = 102;

freqs = [7,14;15,30;31,45];
bands = [{'alpha'},{'beta'},{'low_gamma'}];

st = 3;
en = 6; % in seconds, the feedback period: 3-6s

% bio channel order: ECG, Left EMG, right EMG, Vertical EOG and Horizontal EOG


%% Start Loop

errors = cell(1,4);

parfor i = subjs
    subj = i;
   [ errors{i} ] = wrapper_bio( sessions, condition, subj, bio_dir, data_dir, raw_dir, top_dir, bands, freqs, st, en, nNode );
end