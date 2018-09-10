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
                A = [];
                subj = sprintf('%03d', i);
                for j = 1:numel(sessions)
                    sess = sessions{j};
                    for k = 1:numel(condition)
                        cond = condition{k};
                        load([top_dir, sess, '/', cond, '/', subj, '/FCmatrices/NMF_', freq, '_', eType, '_wpli_lowres.mat']);
                        eval(['A = [A, wpli_', num2str(eType), '];']);
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
                
                behavior = [behavior_updated.BCI.Perf.Sess1.Runs{i}, behavior_updated.BCI.Perf.Sess2.Runs{i}, behavior_updated.BCI.Perf.Sess3.Runs{i}, behavior_updated.BCI.Perf.Sess4.Runs{i}];
                behavior = behavior./mean(behavior);
                behavior = behavior.*sf;
                A = [A; behavior];
                save([save_dir, 'wpli_lowres_', freq, '_', eType, '_', subj], 'A')
            catch
                cnte = cnte + 1;
                errors{cnte} = [top_dir, sess, '/', cond, '/', subj, '/FCmatrices/NMF_', freq, '_', eType, '_wpli_lowres.mat'];
            end
        end
    end
end
