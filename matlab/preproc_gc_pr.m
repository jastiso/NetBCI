%% Preprocessing

% preprocess MEG data, and get into connectivity matrix
% @author JStiso jeni.stio@gmail.com Nov 2017

% CHange log
% June 20, 2018 - changed to baseline correct from intertrial interval and
% added wraper function; switched from time shifted to phase randomized
% null model
% July 1- - changed baseline to be 0-1s


%% Define Global Variables

addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830')
addpath(genpath('/Users/stiso/Documents/MATLAB/eeglab_current'))

top_dir = '/Users/stiso/Documents/MATLAB/NetBCI/';
data_dir = '/Users/stiso/Resilio Sync/NETBCI.RAW/DataBase/1_Signals/2_Segmentation/2_MEG/';
raw_dir = '/Users/stiso/Resilio Sync/NETBCI.RAW/DataBase/1_Signals/0_RawData/2_MEG/';
sessions = [{'Session1'}, {'Session2'}, {'Session3'}, {'Session4'}];
condition = [{'test01'}, {'test02'}, {'test03'}, {'test04'}, {'test05'}, {'test06'}]; % not including rest
subjs = [13];
nNode = 102;
nChunks = 1; %number of sets to divide trials into
freqs = [7,14;15,30;31,45;55,70];
bands = [{'alpha'},{'beta'},{'low_gamma'}, {'gamma'}];
st = 3;
en = 6; % in seconds, the feedback period: 3-6s
bl_st = 0.2; % half o fthe intertrial interval
bl_en = 0.8;
num_neg_grad = [];
num_neg_mag = [];
cnte = 1;
% GC params
dt = 100;
lag = 10;
t0 = dt+lag;

%% Start Loop

errors_pr = cell(1,4);
parfor i = 1:numel(sessions)
    session = sessions{i};
   [ errors_pr{i} ] = wrapper_pr_gc( session, condition, subjs, data_dir, raw_dir, top_dir, bands, freqs, bl_st, bl_en, dt, lag, t0, st, en, nNode );
end

