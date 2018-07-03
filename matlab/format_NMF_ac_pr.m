%% Get data into format ideal for NMF
% namely, edge x time matrix

addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830')
addpath(genpath('/Users/stiso/Documents/MATLAB/eeglab_current'))

top_dir = '/Users/stiso/Documents/MATLAB/NetBCI/';
save_dir = '/Users/stiso/Documents/MATLAB/NetBCI/NMF/';
sessions = [{'Session1'}, {'Session2'}, {'Session3'}, {'Session4'}];
condition = [{'test01'}, {'test02'}, {'test03'}, {'test04'}, {'test05'}, {'test06'}]; % not including rest
subjs = [2:16];
freqs = [3,6;7,14;15,30;31,45;55,70];
bands = [{'all'},{'theta'},{'alpha'},{'beta'},{'low_gamma'}, {'gamma'}];
elecs = [{'grad'}, {'mag'}];
%sf = 102;

% load behavior
load([top_dir, 'Behavior/perf'])

%% Load and concatenate

for e = 1:numel(elecs)
    eType = elecs{e};
    for f  = 1:numel(bands)
        freq = bands{f};
        for i = subjs
            A = [];
            subj = sprintf('%03d', i);
            for j = 1:numel(sessions)
                sess = sessions{j};
                for k = 1:numel(condition)
                    cond = condition{k};
                    load([top_dir, sess, '/', cond, '/', subj, '/FCmatrices/pr/NMF_', freq, '_', eType, '_ac.mat']);
                    if strcmp(freq, 'all')
                        eval(['A = [A, ac_', num2str(eType), '];']);
                    else
                        eval(['A = [A; ac_', num2str(eType), '];']);
                    end
                end
            end
            %scale perf
            eval(['sf = mean(mean(ac_', num2str(eType), '(ac_', num2str(eType), '~=0)));']);
            behavior = perf(:,i)';
            behavior = behavior/mean(behavior);
            behavior = behavior.*sf;
            if strcmp(freq, 'all')
                A = [A; behavior];
            else
            A = [A, behavior']';
            end
            save([save_dir, 'pr_ac_', freq, '_', eType, '_', subj], 'A')
        end
    end
end



