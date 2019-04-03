%% Bin consistency by lobe

% Global Variables

addpath(genpath('/Users/stiso/Documents/MATLAB/npy-matlab-master/'))
addpath('/Users/stiso/Documents/MATLAB/fieldtrip-20170830/')
addpath(genpath('/Users/stiso/Documents/MATLAB/BCT/'))

top_dir = '/Users/stiso/Documents/MATLAB/NetBCI/';
data_dir = '/Users/stiso/Documents/Python/NetBCI/NMF/param/';
save_dir = [top_dir, 'GroupAvg/wpli/analysis/'];
img_dir = [top_dir, 'GroupAvg/wpli/images/'];

%regions = [{'Vertex'},{'Left_frontal'}, {'Frontal'} ,{'Right_frontal'},{'Frontoparietal'},   ...
%    {'Left_parietal'},{'Right_parietal'},{'Parietoccipital'}  {'Left_occipital'}, {'Right_occipital'}, ...
%    {'Left_temporal'},{'Right_temporal'}, {'Left_motor'}, {'Right_motor'}];
regions = [{'Left_frontal'}, {'Left_motor'},{'Left_parietal'} , {'Left_temporal'},{'Left_occipital'},  ...
    {'Right_frontal'}, {'Right_motor'},{'Right_parietal'},{'Right_temporal'}, {'Right_occipital'}];
%regions = [{'Left_frontal'}, {'Right_frontal'},{'Left_motor'}, {'Right_motor'},{'Left_parietal'} ,  ...
%    {'Right_parietal'},{'Left_temporal'},{'Right_temporal'},{'Left_occipital'},  {'Right_occipital'}];

subjs = [1:20];
nSubj = numel(subjs);
nNode = 102;
nEdge = (nNode^2-nNode)/2;
freqs = [7,14;15,30;31,45;55,70];
bands = [{'alpha'},{'beta'},{'low_gamma'}];
sensors = [{'grad'}];

load([save_dir, 'most_consistent_pr.mat'])
high_pr = c_high_max;
high2_pr = c_high2_max;
high3_pr = c_high3_max;
low_pr = c_low_max;


% intialize
high_edges = zeros(numel(regions), numel(regions), numel(bands));
high2_edges = zeros(numel(regions), numel(regions), numel(bands));
high3_edges = zeros(numel(regions), numel(regions), numel(bands));
low_edges = zeros(numel(regions), numel(regions), numel(bands));

high_edges_pr = zeros(numel(regions), numel(regions), numel(bands));
high2_edges_pr = zeros(numel(regions), numel(regions), numel(bands));
high3_edges_pr = zeros(numel(regions), numel(regions), numel(bands));
low_edges_pr = zeros(numel(regions), numel(regions), numel(bands));

labels = cell(numel(regions), numel(regions));

%% Fill in data

for k = 1:numel(bands)
    load([save_dir, 'consistent_mats_', bands{k}, '.mat'])
    high = c_high_avg;
    high2 = c_high2_avg;
    high3 = c_high3_avg;
    low = c_low_avg;
    for i = 1:numel(regions)
        load([top_dir, 'montages/', regions{i}, '_idx.mat'])
        idx1 = logical(idx);
        for j = 1:numel(regions)
            load([top_dir, 'montages/', regions{j}, '_idx.mat'])
            labels{i,j} = [regions{i}, '_', regions{j}];
            idx = logical(idx);
            high_edges(i,j,k) = mean(mean(high(idx1,idx)));
            high2_edges(i,j,k) = mean(mean(high2(idx1,idx)));
            high3_edges(i,j,k) = mean(mean(high3(idx1,idx)));
            low_edges(i,j,k) = mean(mean(low(idx1,idx)));
        end
    end
end

%make colormaps
nColors = 55;
x = [1,round(nColors/2),nColors];
%gamma_cmap1 = [linspace(1,172/255,55)',linspace(1,65/255,55)',linspace(1,0/255,55)'];
red = [1,232/255,172/255]; green = [1,200/255,65/255]; blue = [1,134/255,0/255];
gamma_cmap = [interp1(x,red,1:nColors,'cubic')',interp1(x,green,1:nColors,'cubic')',interp1(x,blue,1:nColors,'cubic')'];

red = [1,192/255,0/255]; green = [1,217/255,89/255]; blue = [1,224/255,115/255];
alpha_cmap = [interp1(x,red,1:nColors,'cubic')',interp1(x,green,1:nColors,'cubic')',interp1(x,blue,1:nColors,'cubic')'];

red = [1,217/255,151/255]; green = [1,143/255,0/255]; blue = [1,147/255,0/255];
beta_cmap = [interp1(x,red,1:nColors,'cubic')',interp1(x,green,1:nColors,'cubic')',interp1(x,blue,1:nColors,'cubic')'];

for i = 1:numel(bands)
    if i == 1
        %cmap = cbrewer("seq", "Blues", 55);
        cmap = alpha_cmap;
    elseif i == 2
        %cmap = cbrewer("seq", "Reds", 55);
        cmap = beta_cmap;
    elseif i == 3
        cmap = gamma_cmap;
    end
            
    imagesc(high_edges(:,:,i));colorbar
    title([bands{i}, '_high']);
    %caxis([1,6]);
    colormap(cmap);
    saveas(gca, [img_dir, bands{i}, '_high_lobes.png'], 'png')
    
    imagesc(high2_edges(:,:,i));colorbar
    title([bands{i}, '_high2']);
    %caxis([1,6]);
    colormap(cmap);
    saveas(gca, [img_dir, bands{i}, '_high2_lobes.png'], 'png')
    
    imagesc(high3_edges(:,:,i));colorbar
    title([bands{i}, '_high3']);
    %caxis([1,6]);
    colormap(cmap);
    saveas(gca, [img_dir, bands{i}, '_high3_lobes.png'], 'png')
    
    imagesc(low_edges(:,:,i));colorbar
    title([bands{i}, '_low']);
    %caxis([1,6]);
    colormap(cmap);
    saveas(gca, [img_dir, bands{i}, '_low_lobes.png'], 'png')
    
end

%% Get strongest edges?

% intiailize, cell arrays
high_ranked = cell((numel(regions)^2 - numel(regions))/2 + numel(regions), numel(bands)*2);
high2_ranked = cell((numel(regions)^2 - numel(regions))/2 + numel(regions), numel(bands)*2);
high3_ranked = cell((numel(regions)^2 - numel(regions))/2 + numel(regions), numel(bands)*2);
low_ranked = cell((numel(regions)^2 - numel(regions))/2 + numel(regions), numel(bands)*2);

for i = 1:numel(bands)
    curr_h = high_edges(:,:,i);
    curr_2 = high2_edges(:,:,i);
    curr_3 = high3_edges(:,:,i);
    curr_l = low_edges(:,:,i);
    
    [edge_h,idx_h] = sort(curr_h(triu(true(10))),'descend');
    [edge_2,idx_2] = sort(curr_2(triu(true(10))),'descend');
    [edge_3,idx_3] = sort(curr_3(triu(true(10))),'descend');
    [edge_l,idx_l] = sort(curr_l(triu(true(10))),'descend');
   
    v_labels = labels(triu(true(10)));
    
    high_ranked(:,i) = v_labels(idx_h);
    high2_ranked(:,i) = v_labels(idx_2);
    high3_ranked(:,i) = v_labels(idx_3);
    low_ranked(:,i) = v_labels(idx_l);
    high_ranked(:,i+numel(bands)) = num2cell(edge_h);
    high2_ranked(:,i+numel(bands)) = num2cell(edge_2);
    high3_ranked(:,i+numel(bands)) = num2cell(edge_3);
    low_ranked(:,i+numel(bands)) = num2cell(edge_l);

end

save([save_dir, 'lobe_ranks'], 'high_ranked','high2_ranked','high3_ranked','low_ranked')

%% Fill in data - PR

for k = 1:numel(bands)
    load([save_dir, 'consistent_mats_', bands{k}, '_pr.mat'])
    high = c_high_avg;
    high2 = c_high2_avg;
    high3 = c_high3_avg;
    low = c_low_avg;
    for i = 1:numel(regions)
        load([top_dir, 'montages/', regions{i}, '_idx.mat'])
        idx1 = logical(idx);
        for j = 1:numel(regions)
            load([top_dir, 'montages/', regions{j}, '_idx.mat'])
            labels{i,j} = [regions{i}, '_', regions{j}];
            idx = logical(idx);
            high_edges(i,j,k) = mean(mean(high(idx1,idx)));
            high2_edges(i,j,k) = mean(mean(high2(idx1,idx)));
            high3_edges(i,j,k) = mean(mean(high3(idx1,idx)));
            low_edges(i,j,k) = mean(mean(low(idx1,idx)));
        end
    end
end

%make colormaps
nColors = 55;
x = [1,round(nColors/2),nColors];
%gamma_cmap1 = [linspace(1,172/255,55)',linspace(1,65/255,55)',linspace(1,0/255,55)'];
red = [1,232/255,172/255]; green = [1,200/255,65/255]; blue = [1,134/255,0/255];
gamma_cmap = [interp1(x,red,1:nColors,'cubic')',interp1(x,green,1:nColors,'cubic')',interp1(x,blue,1:nColors,'cubic')'];

red = [1,192/255,0/255]; green = [1,217/255,89/255]; blue = [1,224/255,115/255];
alpha_cmap = [interp1(x,red,1:nColors,'cubic')',interp1(x,green,1:nColors,'cubic')',interp1(x,blue,1:nColors,'cubic')'];

red = [1,217/255,151/255]; green = [1,143/255,0/255]; blue = [1,147/255,0/255];
beta_cmap = [interp1(x,red,1:nColors,'cubic')',interp1(x,green,1:nColors,'cubic')',interp1(x,blue,1:nColors,'cubic')'];

for i = 1:numel(bands)
    if i == 1
        %cmap = cbrewer("seq", "Blues", 55);
        cmap = alpha_cmap;
    elseif i == 2
        %cmap = cbrewer("seq", "Reds", 55);
        cmap = beta_cmap;
    elseif i == 3
        cmap = gamma_cmap;
    end
            
    imagesc(high_edges(:,:,i));colorbar
    title([bands{i}, '_high']);
    %caxis([1,6]);
    colormap(cmap);
    saveas(gca, [img_dir, bands{i}, '_high_lobes_pr.png'], 'png')
    
    imagesc(high2_edges(:,:,i));colorbar
    title([bands{i}, '_high2']);
    %caxis([1,6]);
    colormap(cmap);
    saveas(gca, [img_dir, bands{i}, '_high2_lobes_pr.png'], 'png')
    
    imagesc(high3_edges(:,:,i));colorbar
    title([bands{i}, '_high3']);
    %caxis([1,6]);
    colormap(cmap);
    saveas(gca, [img_dir, bands{i}, '_high3_lobes_pr.png'], 'png')
    
    imagesc(low_edges(:,:,i));colorbar
    title([bands{i}, '_low']);
    %caxis([1,6]);
    colormap(cmap);
    saveas(gca, [img_dir, bands{i}, '_low_lobes_pr.png'], 'png')
    
end

%% Get strongest edges?

% intiailize, cell arrays
high_ranked = cell((numel(regions)^2 - numel(regions))/2 + numel(regions), numel(bands)*2);
high2_ranked = cell((numel(regions)^2 - numel(regions))/2 + numel(regions), numel(bands)*2);
high3_ranked = cell((numel(regions)^2 - numel(regions))/2 + numel(regions), numel(bands)*2);
low_ranked = cell((numel(regions)^2 - numel(regions))/2 + numel(regions), numel(bands)*2);

for i = 1:numel(bands)
    curr_h = high_edges(:,:,i);
    curr_2 = high2_edges(:,:,i);
    curr_3 = high3_edges(:,:,i);
    curr_l = low_edges(:,:,i);
    
    [edge_h,idx_h] = sort(curr_h(triu(true(10))),'descend');
    [edge_2,idx_2] = sort(curr_2(triu(true(10))),'descend');
    [edge_3,idx_3] = sort(curr_3(triu(true(10))),'descend');
    [edge_l,idx_l] = sort(curr_l(triu(true(10))),'descend');
   
    v_labels = labels(triu(true(10)));
    
    high_ranked(:,i) = v_labels(idx_h);
    high2_ranked(:,i) = v_labels(idx_2);
    high3_ranked(:,i) = v_labels(idx_3);
    low_ranked(:,i) = v_labels(idx_l);
    high_ranked(:,i+numel(bands)) = num2cell(edge_h);
    high2_ranked(:,i+numel(bands)) = num2cell(edge_2);
    high3_ranked(:,i+numel(bands)) = num2cell(edge_3);
    low_ranked(:,i+numel(bands)) = num2cell(edge_l);

end

save([save_dir, 'lobe_ranks_pr'], 'high_ranked','high2_ranked','high3_ranked','low_ranked')