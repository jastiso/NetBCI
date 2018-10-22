function [vect,val, small_vect, vect2] = easy_state(A,B,C,D)
% This calculates the easiest to reach state for a given system (see help
% sys)
% essentially it's the largest eigen vector of the grammian

% get system
sys = ss(A,B,C,D);
% get grammian           
grammian = gram(sys, 'c');
            
% get largest eigenvector
[Vect, val] = eigs(grammian,2);
vect = Vect(:,1);
vect2 = Vect(:,2);

% get smallest eigen vector
[V, vals] = eig(grammian);
vals = diag(vals);
[~,idx] = min(vals);
small_vect = V(:,idx);

end

