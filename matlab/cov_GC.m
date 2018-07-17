function [GC, pairs] = covGC_time(X, dt, lag, t0)
% [GC, pairs] = covGC_time(X, dt, lag, t0)
%
% Computes single-trials covariance-based Granger Causality for gaussian variables
% 
% X   = data arranged as channels x samples
% dt  = duration of the time window for covariance correlation in samples
% lag = number of samples for the lag within each trial
% t0  = zero time in samples
%
% GC  = Granger Causality arranged as (number of pairs) x (3 directionalities (pair(:,1)->pair(:,2), pair(:,2)->pair(:,1), instantaneous))
% pairs = indices of sources arranged as number of pairs x 2
%
% -------------------- Total Granger interdependence ----------------------
% Total Granger interdependence:
% TGI = GC(x,y)
% TGI = sum(GC,2):
% TGI = GC(x->y) + GC(y->x) + GC(x.y)
% TGI = GC(x->y) + GC(y->x) + GC(x.y) = Hycy + Hxcx - Hxxcyy
% This quantity can be defined as the Increment of Total
% Interdependence and it can be calculated from the different of two
% mutual informations as follows
%
% ----- Relations between Mutual Informarion and conditional entropies ----
% % I(X_i+1,X_i|Y_i+1,Y_i) = H(X_i+1) + H(Y_i+1) - H(X_i+1,Y_i+1)
% Ixxyy   = log(det_xi1) + log(det_yi1) - log(det_xyi1);
% % I(X_i|Y_i) = H(X_i) + H(Y_i) - H(X_i, Y_i)
% Ixy     = log(det_xi) + log(det_yi) - log(det_yxi);
% ITI(np) = Ixxyy - Ixy;
%
% Reference
% Brovelli A, Chicharro D, Badier JM, Wang H, Jirsa V (2015) JNeurosci
%
% Copyright of Andrea Brovelli (Jan 2015)

% Change Log
% JStiso: added warning for returning complex numbers

% Data parameters. Size = sources x time points
[nSo, nTi] = size(X);

% Select a single window according to index t0
ind_t = [t0-dt+1 : t0]';

% Create indeces for all lags
ind_t = repmat(ind_t, [1 lag+1]) - repmat( [0:lag], [dt 1]);

% Pairs between sources
[pairs(:,1), pairs(:,2)] = find( tril(ones(nSo),-1) == 1 );
nPairs = size(pairs,1);

% Init
GC    = zeros(nPairs,3);
count = 1; tic

% Normalisation coefficient for gaussian entropy
C = log(2*pi*exp(1));

% Loop over number of pairs
for np = 1:nPairs
    
    % Extract data for a given pair of sources
    x = squeeze(X(pairs(np,1),ind_t));
    y = squeeze(X(pairs(np,2),ind_t));
    % Reshape to trials x dt x lags
    x = reshape(x, [dt lag+1]);
    y = reshape(y, [dt lag+1]);
    
    % ---------------------------------------------------------------------
    % Conditional Entropies
    % ---------------------------------------------------------------------
    % Hycy: H(Y_i+1|Y_i) = H(Y_i+1) - H(Y_i)
    det_yi1  = det(cov(y));
    det_yi   = det(cov(y(:,2:end)));
    Hycy     = log(det_yi1) - log(det_yi); 
    % Hycx: H(Y_i+1|X_i,Y_i) = H(Y_i+1,X_i,Y_i) - H(X_i,Y_i)
    det_yxi1 = det(cov([ y x(:,2:end) ]));
    det_yxi  = det(cov([ y(:,2:end) x(:,2:end) ]));
    Hycx     = log(det_yxi1) - log(det_yxi);
    % Hxcx: H(X_i+1|X_i) = H(X_i+1) - H(X_i)
    det_xi1  = det(cov(x));
    det_xi   = det(cov(x(:,2:end)));
    Hxcx     = log(det_xi1) - log(det_xi);
    % Hxcy: H(X_i+1|X_i,Y_i) = H(X_i+1,X_i,Y_i) - H(X_i,Y_i)
    det_xyi1 = det(cov([ x y(:,2:end) ]));
    Hxcy     = log(det_xyi1) - log(det_yxi);      
    % Hxxcyy: H(X_i+1,Y_i+1|X_i,Y_i) = H(X_i+1,Y_i+1,X_i,Y_i) - H(X_i,Y_i)
    det_xyi1 = det(cov([ x y ]));
    Hxxcyy   = log(det_xyi1) - log(det_yxi);
    
    % ---------------------------------------------------------------------
    % Causality measures
    % ---------------------------------------------------------------------   
    % GC(pairs(:,1)->pairs(:,2))
    GC(np,1) = Hycy - Hycx;
    % GC(pairs(:,2)->pairs(:,1))
    GC(np,2) = Hxcx - Hxcy;
    % GC(x.y)
    GC(np,3) = Hycx + Hxcy - Hxxcyy;
    
    % Counter
    if count/100000 == 1
        tempo = toc;
        display([' Elapsed time ' num2str(tempo) 's   % nPairs = ' num2str(round(np/nPairs*1000)/10) ])
        count = 0;        tic
    end
    count = count + 1;
    
    if any(any(~isreal(GC)))
        warning('Your GC has imaginary numbers. This might be due precision errors in calculating the determinant.')
    end
end

