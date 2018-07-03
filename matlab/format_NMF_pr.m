%% Get data into format ideal for NMF for phase rand
% namely, edge x time matrix

addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830')
addpath(genpath('/Users/stiso/Documents/MATLAB/eeglab_current'))

top_dir = '/Users/stiso/Documents/MATLAB/NetBCI/';
save_dir = '/Users/stiso/Documents/MATLAB/NetBCI/NMF/pr/';
sessions = [{'Session1'}, {'Session2'}, {'Session3'}, {'Session4'}];
condition = [{'test01'}, {'test02'}, {'test03'}, {'test04'}, {'test05'}, {'test06'}]; % not including rest
subjs = [2:16];
freqs = [3,6;7,14;15,30;31,45;55,70];
bands = [{'theta'},{'alpha'},{'beta'},{'low_gamma'},{'gamma'}];
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
                    load([top_dir, sess, '/', cond, '/', subj, '/pr/FCmatrices/NMF_', freq, '_', eType, '.mat']);
                    eval(['A = [A, conn_tv_', num2str(eType), '];']);
                end
            end
            %scale perf
            eval(['sf = mean(mean(conn_tv_', num2str(eType), '(conn_tv_', num2str(eType), '~=0)));']);
            behavior = perf(:,i)';
            behavior = behavior/mean(behavior);
            behavior = behavior.*sf;
            A = [A; behavior];
            save([save_dir, freq, '_', eType, '_', subj], 'A')
        end
    end
end



%% make combined freq

for e = 1:numel(elecs)
    eType = elecs{e};
    for i = subjs
        B = [];
            subj = sprintf('%03d', i);
        for f  = 1:numel(bands)
        freq = bands{f};
        load([save_dir, freq, '_', eType, '_', subj])
        B = [B;A(1:end-1,:)];
        
        end
        B = [B; A(end,:)];
        % needs to be called A for python pipeline
        A = B;
        save([save_dir, 'all_', eType, '_', subj], 'A')
    end
    
end

%% Repeat for baseline

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
                    load([top_dir, sess, '/', cond, '/', subj, '/pr/FCmatrices/NMF_', freq, '_', eType, '_bl.mat']);
                    eval(['A = [A, conn_tv_', num2str(eType), '_bl];']);
                end
            end
            %scale perf
            eval(['sf = mean(mean(conn_tv_', num2str(eType), '_bl(conn_tv_', num2str(eType), '_bl~=0)));']);
            behavior = perf(:,i)';
            behavior = behavior/mean(behavior);
            behavior = behavior.*mean(mean(A(A ~=0)));
            behavior = behavior.*sf;
            A = [A; behavior];
            save([save_dir, 'bl_', freq, '_', eType, '_', subj], 'A')
        end
    end
end



%% make combined freq baseline

for e = 1:numel(elecs)
    eType = elecs{e};
    for i = subjs
        B = [];
            subj = sprintf('%03d', i);
        for f  = 1:numel(bands)
        freq = bands{f};
        load([save_dir, 'bl_', freq, '_', eType, '_', subj])
        B = [B;A(1:end-1,:)];
        
        end
        B = [B; A(end,:)];
        % needs to be called A for python pipeline
        A = B;
        save([save_dir, 'bl_all_', eType, '_', subj], 'A')
    end
    
end
