%% Get data into format ideal for NMF
% namely, edge x time matrix

addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830')
addpath(genpath('/Users/stiso/Documents/MATLAB/eeglab_current'))

top_dir = '/Users/stiso/Documents/MATLAB/NetBCI/';
save_dir = '/Users/stiso/Documents/MATLAB/NetBCI/NMF/';
sessions = [{'Session1'}, {'Session2'}, {'Session3'}, {'Session4'}];
condition = [{'test01'}, {'test02'}, {'test03'}, {'test04'}, {'test05'}, {'test06'}]; % not including rest
subjs = [1:20];
freqs = [3,6;7,14;15,30;31,45;55,70];
bands = [{'alpha'},{'beta'},{'low_gamma'}, {'gamma'}];
elecs = [{'grad'}];
%sf = 102;
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
                        % next time you run this, change from pc to gc. pc
                        % was a typo in the preproc script
                        load([top_dir, sess, '/', cond, '/', subj, '/FCmatrices/NMF_', freq, '_', eType, '_gc.mat']);
                        eval(['A = [A, gc_', num2str(eType), '];']);
                        eval(['nTrial(j,k) = size(gc_', num2str(eType), ',2);'])
                        if any(any(isnan(A)))
                            display([num2str(i), ' ', num2str(f)])
                            warning('This matrix has NaNs')
                        end
                    end
                end
                %scale perf
                eval(['sf = mean(mean(gc_', num2str(eType), '(gc_', num2str(eType), '~=0)));']);
                
                behavior = interp_perf(behavior_updated, i, 6, nTrial);
                behavior_all{i} = behavior;
                behavior = behavior.*sf;
                A = [A; behavior];
                save([save_dir, 'gc_', freq, '_', eType, '_', subj], 'A')
            catch
                cnte = cnte + 1;
                errors{cnte} = [top_dir, sess, '/', cond, '/', subj, '/FCmatrices/NMF_', freq, '_', eType, '_gc.mat'];
            end
        end
    end
end

save([top_dir, 'Behavior/behavior_all.mat'], 'behavior_all');

%% Repeat for time shift

cnte = 0;
errors_ts = [];
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
                        load([top_dir, sess, '/', cond, '/', subj, '/FCmatrices/ts/NMF_', freq, '_', eType, '_gc.mat']);
                        
                        eval(['A = [A, gc_', num2str(eType), '];']);
                        eval(['nTrial(j,k) = size(gc_', num2str(eType), ',2);'])
                         
                    end
                end
                %scale perf
                eval(['sf = mean(mean(gc_', num2str(eType), '(gc_', num2str(eType), '~=0)));']);
                
                behavior = interp_perf(behavior_updated, i, 6, nTrial);
                behavior = behavior.*sf;
                A = [A; behavior];
                save([save_dir, 'ts_gc_', freq, '_', eType, '_', subj], 'A')
                if any(any(isnan(A)))
                    display([num2str(i), ' ', num2str(f)])
                    warning('This matrix has NaNs')
                end
            catch
                cnte = cnte + 1;
                errors_ts{cnte} = [top_dir, sess, '/', cond, '/', subj, '/FCmatrices/NMF_', freq, '_', eType, 'ts_gc.mat'];
            end
        end
    end
end



