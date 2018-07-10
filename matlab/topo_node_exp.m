%% Topoplot Subgraph expression
% make a topoplot of each nodes total participation in each subgraph
% this moves things into the nodes space, and might washout effects of edge
% space


%% Global Variables

addpath(genpath('/Users/stiso/Documents/MATLAB/npy-matlab-master/'))
addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/')

top_dir = '/Users/stiso/Documents/MATLAB/NetBCI/';
data_dir = '/Users/stiso/Documents/Python/NetBCI/NMF/gc/';

subjs = [1,2,6,7,9:20];
nSubj = numel(subjs);
nNode = 102;
nEdge = (nNode^2-nNode)/2;
freqs = [7,14;15,30;31,45;55,70];
bands = [{'gamma'}];
sensors = [{'grad'}];
% threshold for edge visualization
thr = 1;
consensus = zeros(nNode,nNode, nSubj-1, numel(bands));

%% Loop through data

for i = subjs
    subj = sprintf('%03d', i);
    for j = 1:numel(sensors)
        sens = sensors{j};
        load([top_dir, 'idx_', sens, '.mat'])
        for k = 1:numel(bands)
            f = bands{k};
            img_dir = [top_dir, 'NMF/', subj, '/', sens, '/images/'];
            % make directories
            if ~exist(img_dir, 'dir')
                mkdir(img_dir)
            end
            save_dir = [top_dir, 'NMF/', subj, '/', sens, '/analysis/'];
            % make directories
            if ~exist(save_dir, 'dir')
                mkdir(save_dir)
            end
            
            % get subgraph data
            coeff = readNPY([data_dir, subj,  '/', sens, '/gc_', f, '_coeff.npy']);
            subset = readNPY([data_dir, subj,  '/', sens, '/gc_', f, '_subset.npy']);
            err = readNPY([data_dir, subj,  '/', sens, '/gc_', f, '_err.npy']);
            % get mag index
            labels = [];
            try
                load(['/Users/stiso/Resilio Sync/NETBCI.RAW/DataBase/1_Signals/0_RawData/2_MEG/Session1/1_test01/RawData_MEG_Subj', subj, '_Ses1_test01.mat'])
                eval(['data = RawData_MEG_Subj', subj, '_Ses1_test01;']);
                eval(['clear RawData_MEG_Subj', subj, '_Ses1_test01'])
            catch
                load(['/Users/stiso/Resilio Sync/NETBCI.RAW/DataBase/1_Signals/0_RawData/2_MEG/Session1/1_test01/RawData_MEG_Subj016_Ses1_test01.mat'])
                eval(['data = RawData_MEG_Subj016_Ses1_test01;']);
                eval(['clear RawData_MEG_Subj016_Ses1_test01'])
            end
            if strcmp(sens,'mag')
                sens_idx = strcmp(data.meg_sensors.chantype,'megmag');
                labels = data.label(sens_idx);
            else
                sens_idx = strcmp(data.meg_sensors.chantype,'megplanar');
               grad_label = data.label(sens_idx);
                cnt = 1;
                for n = 1:2:sum(sens_idx)
                    labels{cnt,1} = [grad_label{n}, '+', grad_label{n+1}(end-3:end)];
                    cnt = cnt + 1;
                end
            end
            
            
            % get subgraph with highest expresssion in behavior
            [~,bsg] = max(subset(:,end));
            
            % sizes
            nTime = size(coeff,2);
            nSG = size(subset,1);
            nSens = numel(labels);
            
            for n = 1:nSG
                % get expressin into matrix
                curr_SG = subset(n,1:end-1);
                node_exp = get_sg_matrix(nSens, curr_SG);
                sg_deg = sum(node_exp)';
                if n == bsg
                    consensus(:,:,i,k) = node_exp;
                end
                % save edge file
                % first threshold
                v_graph = reshape(node_exp((node_exp ~= 0)), 1, []);
                v_graph = sort(v_graph);
                n_thrsh = floor(numel(v_graph)*thr);
                v_thrsh = v_graph(end-n_thrsh+1);
                node_exp(node_exp<v_thrsh) = 0;
                dlmwrite([save_dir, f, '_graph_', num2str(n), '.edge'], node_exp, '\t')
                
                % get into fieldtrip format for topoplot
                cfg = [];
                cfg.style = 'straight';
                if strcmp('mag',sens)
                    cfg.layout = [top_dir, 'layouts/neuromag306', sens, '.lay'];
                    fid = fopen([top_dir, 'layouts/neuromag306', sens, '.txt'], 'r');
                else
                    cfg.layout = [top_dir, 'layouts/neuromag306cmb.lay'];
                    fid = fopen([top_dir, 'layouts/neuromag306cmb.txt'], 'r');
                end

                
                plot_data.powspctrm = sg_deg;
                plot_data.label = labels;
                plot_data.dimord = 'chan_freq';
                plot_data.freq = 6;
                plot_data.cfg = [];
                figure(1); clf
                ft_topoplotER(cfg,plot_data); colorbar
                saveas(gca, [img_dir, f, '_', num2str(n), '.png'], 'png')
                
                % sidebar - make node file
                bv_label = textscan(fid, '%d %10f %10f %d %d %s', 'Delimiter', '\n');
                fclose(fid);
                write_bv_node([save_dir, f, '_graph_', num2str(n), '.node'], bv_label{:,2},bv_label{:,3},double(bv_label{:,4}), 1, sg_deg);
            end
            
        end
    end
end


%% plot consensus - consistency (but without distance)
thrs = linspace(0.01, 1, 100);
nE_thr = zeros(numel(thrs),numel(bands));
consistent_edges = zeros(numel(thrs), numel(bands));
for k = 1:numel(thrs)
    thrsh_mats = consensus;
    for i = 1:nSubj
        for j = 1:numel(bands)
            thresh = thrs(k);
            curr = thrsh_mats(:,:,i,j);
            thresh = round(thresh*sum(sum(curr ~= 0)));
            tmp = curr(curr ~= 0);
            tmp = sort(tmp, 'descend');
            curr(curr < tmp(thresh)) = 0;
            nE_thr(k,j) = sum(sum(curr ~= 0));
            thrsh_mats(:,:,i,j) = curr;
        end
    end
    for m = 1:numel(bands)
        G = fcn_group_avg4(thrsh_mats(:,:,:,m),nSubj-1);
        consistent_edges(k,m) = ((sum(sum(G))/2)./nE_thr(k,m))*100;
        if k == numel(thrs)
            avg_consensus = mean(thrsh_mats(:,:,:,m),3);
            avg_consensus(~logical(G)) = 0;
            plot_data.powspctrm = sum(avg_consensus)';
            plot_data.label = labels;
            plot_data.dimord = 'chan_freq';
            plot_data.freq = 6;
            plot_data.cfg = [];
            figure(1); clf
            ft_topoplotER(cfg,plot_data); colorbar
            saveas(gca, ['/Users/stiso/Documents/MATLAB/NetBCI/GroupAvg/gc/images/', bands{m}, '_consensus_topo.png'], 'png')
        end
    end
end
% plot
plot(thrs,consistent_edges, 'linewidth', 2)
legend(bands)
saveas(gca, [img_dir, 'consensus_thresh.png'], 'png')

