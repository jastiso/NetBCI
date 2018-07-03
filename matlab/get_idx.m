%% NMF index
% make index to move things between electrode coords and nmf coords

nNode = 204;
nEdge = (nNode^2-nNode)/2;
idx = zeros(nNode); 
cnt = 0;
for col = 1:nNode
    for row = 1:col
        if row ~= col
            cnt = cnt + 1;
            idx(row,col) = cnt;
        end
    end
end
imagesc(idx)
% check if this setup is correct
tmp = ones(nNode);
vect = reshape(idx(logical(triu(tmp,1))),1,[]);
save('/Users/stiso/Documents/MATLAB/NetBCI/idx_grad.mat', 'idx')



nNode = 102;
nEdge = (nNode^2-nNode)/2;
idx = zeros(nNode); 
cnt = 0;
for col = 1:nNode
    for row = 1:col
        if row ~= col
            cnt = cnt + 1;
            idx(row,col) = cnt;
        end
    end
end
imagesc(idx)
% check if this setup is correct
tmp = ones(nNode);
vect = reshape(idx(logical(triu(tmp,1))),1,[]); % how I reshape this in the data
save('/Users/stiso/Documents/MATLAB/NetBCI/idx_mag.mat', 'idx')

