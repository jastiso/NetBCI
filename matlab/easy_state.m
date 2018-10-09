function [vect,val, small_vect] = easy_state(A,B,C,D)
% This calculates the easiest to reach state for a given system (see help
% sys)
% essentially it's the largest eigen vector of the grammian

% get system
sys = ss(A,B,C,D);
% get grammian           
grammian = gram(sys, 'c');
            
% get largest eigenvector
[vect, val] = eigs(grammian,1);

% get smallest eigen vector
[V, vals] = eig(grammian);
vals = diag(vals);
[~,idx] = min(vals);
small_vect = V(:,idx);

end

