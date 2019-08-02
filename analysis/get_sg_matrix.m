function [ node_exp ] = get_sg_matrix( nSens, curr_SG )
% turn in to matrix
% curr_CG is upper triangle - diagonal
node_exp = zeros(nSens);
cnt = 0;
for col = 1:nSens
    for row = 1:col
        if row ~= col
            cnt = cnt + 1;
            node_exp(row,col) = curr_SG(cnt);
        end
    end
end
node_exp = node_exp+node_exp';
end

