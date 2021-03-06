%% Get data into format ideal for NMF
% namely, edge x time matrix

addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830')
addpath(genpath('/Users/stiso/Documents/MATLAB/eeglab_current'))

top_dir = '/Users/stiso/Documents/MATLAB/NetBCI/';
save_dir = '/Users/stiso/Documents/MATLAB/NetBCI/NMF/';
sessions = [{'Session1'}, {'Session2'}, {'Session3'}, {'Session4'}];
condition = [{'test01'}, {'test02'}, {'test03'}, {'test04'}, {'test05'}, {'test06'}]; % not including rest
subjs = [1:20];
freqs = [3,6;7,14;15,30;31,45];
bands = [{'alpha'},{'beta'},{'low_gamma'}];
elecs = [{'grad'}];

errors = [];
behavior_all = {};
% load behavior
load([top_dir, 'Behavior/behavior_20Subjects.mat']);

%% Load and concatenate

cnte = 0;
for e = 1:numel(elecs)
    eType = elecs{e};
    for f  = 1:numel(bands)
        freq = bands{f};
        for i = subjs
            try
                % get trial numbers
                nTrial = zeros(4,numel(condition));
                A = [];
                subj = sprintf('%03d', i);
                for j = 1:numel(sessions)
                    sess = sessions{j};
                    for k = 1:numel(condition)
                        cond = condition{k};
                        % next time you run this, change from pc to wpli. pc
                        % was a typo in the preproc script
                        load([top_dir, sess, '/', cond, '/', subj, '/FCmatrices/NMF_', freq, '_', eType, '_wpli.mat']);
                        eval(['A = [A, wpli_', num2str(eType), '];']);
                        eval(['nTrial(j,k) = size(wpli_', num2str(eType), ',2);'])
                    end
                end
                if any(any(isnan(A)))
                    display([num2str(i), ' ', num2str(f)])
                    warning('This matrix has NaNs')
                end
                if any(any(~isreal(A)))
                    display([num2str(i), ' ', num2str(f)])
                    warning('This matrix has complex numbers')
                end
                
                %scale perf
                % the mean edge weight, graphs are fully weighted, so non-zeros gets rid of lower diagonal
                sf = mean(mean((nonzeros(triu(A,1)))));
                
                behavior = interp_perf(behavior_updated, i, 6, nTrial);
                behavior_all{i} = behavior;
                behavior = behavior./mean(behavior);
                behavior = behavior.*sf;
                A = [A; behavior];
                save([save_dir, 'wpli_', freq, '_', eType, '_', subj], 'A')
            catch
                cnte = cnte + 1;
                errors{cnte} = [top_dir, sess, '/', cond, '/', subj, '/FCmatrices/NMF_', freq, '_', eType, '_wpli.mat'];
            end
        end
    end
end

save([top_dir, 'Behavior/behavior_all.mat'], 'behavior_all');

%% Repeat for uniform phase randomization

cnte = 0;
errors_pr = [];
nTrial = 0;

for e = 1:numel(elecs)
    eType = elecs{e};
    for f  = 1:numel(bands)
        freq = bands{f};
        for i = subjs
            try
                A = [];
                subj = sprintf('%03d', i);
                for j = 1:numel(sessions)
                    sess = sessions{j};
                    for k = 1:numel(condition)
                        cond = condition{k};
                        load([top_dir, sess, '/', cond, '/', subj, '/FCmatrices/NMF_', freq, '_', eType, '_wpli_pr.mat']);
                        
                        eval(['A = [A, wpli_', num2str(eType), '];']);
                        eval(['nTrial(j,k) = size(wpli_', num2str(eType), ',2);'])
                    end
                    
                end
                if any(any(isnan(A)))
                    display([num2str(i), ' ', num2str(f)])
                    warning('This matrix has NaNs')
                end
                if any(any(~isreal(A)))
                    display([num2str(i), ' ', num2str(f)])
                    warning('This matrix has complex numbers')
                end
                %scale perf
                sf = mean(mean((nonzeros(triu(A,1)))));  
                
                behavior = interp_perf(behavior_updated, i, 6, nTrial);
                behavior_all{i} = behavior;
                behavior = behavior./mean(behavior);
                behavior = behavior.*sf;
                A = [A; behavior];
                save([save_dir, 'pr_wpli_', freq, '_', eType, '_', subj], 'A')
            catch
                cnte = cnte + 1;
                errors_pr{cnte} = [top_dir, sess, '/', cond, '/', subj, '/FCmatrices/NMF_', freq, '_', eType, 'pr_wpli.mat'];
            end
        end
    end
end



%% Repeat for independent phase randomiztion

cnte = 0;
errors_ind = [];
nTrial = 0;

for e = 1:numel(elecs)
    eType = elecs{e};
    for f  = 1:numel(bands)
        freq = bands{f};
        for i = subjs
            try
                A = [];
                subj = sprintf('%03d', i);
                for j = 1:numel(sessions)
                    sess = sessions{j};
                    for k = 1:numel(condition)
                        cond = condition{k};
                        load([top_dir, sess, '/', cond, '/', subj, '/FCmatrices/NMF_', freq, '_', eType, '_ind.mat']);
                        
                        eval(['A = [A, gc_', num2str(eType), '];']);
                        eval(['nTrial(j,k) = size(gc_', num2str(eType), ',2);'])
                        
                    end
                    
                end
                if any(any(isnan(A)))
                    display([num2str(i), ' ', num2str(f)])
                    warning('This matrix has NaNs')
                end
                if any(any(~isreal(A)))
                    display([num2str(i), ' ', num2str(f)])
                    warning('This matrix has complex numbers')
                end
                %scale perf
                sf = mean(mean((nonzeros(triu(A,1)))));  
                
                
                behavior = interp_perf(behavior_updated, i, 6, nTrial);
                behavior_all{i} = behavior;
                behavior = behavior./mean(behavior);
                behavior = behavior.*sf;
                A = [A; behavior];
                save([save_dir, 'ind_gc_', freq, '_', eType, '_', subj], 'A')
            catch
                cnte = cnte + 1;
                errors_ind{cnte} = [top_dir, sess, '/', cond, '/', subj, '/FCmatrices/NMF_', freq, '_', eType, 'ind_gc.mat'];
            end
        end
    end
end