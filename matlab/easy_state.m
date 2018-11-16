function [vect,val, small_vect, vect2, vect3] = easy_state(A,B,C,D)
% This calculates the easiest to reach state for a given system (see help
% sys)
% essentially it's the largest eigenvector of the grammian

% get system
sys = ss(A,B,C,D);
% get grammian           
grammian = gram(sys, 'c');
            
% get largest eigenvector
[Vect, val] = eigs(grammian,4);
vect = Vect(:,1);
vect2 = Vect(:,2);
vect3 = ones(size(vect));%Vect(:,3);

% get smallest eigen vector
[V, vals] = eig(grammian);
vals = diag(vals);
[~,idx] = min(vals);
small_vect = V(:,idx);

end

