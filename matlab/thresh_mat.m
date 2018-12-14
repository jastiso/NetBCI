function [A_thresh] = thresh_mat(A,thr)
%remove the thr% weakest connections in A
% this assumes A is symmetric

[edges] = sort(nonzeros(tril(A,-1)),'ascend');
nEdge = round(numel(edges)*thr);
rm_edges = edges(1:nEdge);
A_thresh = tril(A,-1);
for i = 1:numel(rm_edges)
    A_thresh(A_thresh == rm_edges(i)) = 0;
end
A_thresh = A_thresh+A_thresh';

end

