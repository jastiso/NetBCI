n = 5152;
x = linspace(0.01,1,24);
A = zeros(n,numel(x));
for i = 1:n
    m = (rand(1) - 0.5)./0.5;
    A(i,:) = x.*m;
end

save('/Users/stiso/Documents/MATLAB/NetBCI/NMF/reg_A.mat', 'A');

A = rand(n,numel(x));

save('/Users/stiso/Documents/MATLAB/NetBCI/NMF/rand_A.mat', 'A');