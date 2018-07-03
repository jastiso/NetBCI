function [ perf ] = interp_perf( behavior_updated, i, nRuns, nTp )
perf = [];

% sess 1
tmp = behavior_updated.BCI.Perf.Sess1.Runs{i};
interp_tmp = [];
for i = 1:nRuns
    interp_tmp = [interp_tmp,repmat(tmp(i),1,nTp(1,i))];
end
perf = [perf,interp_tmp];


% sess 2
tmp = behavior_updated.BCI.Perf.Sess2.Runs{i};
interp_tmp = [];
for i = 1:nRuns
    interp_tmp = [interp_tmp,repmat(tmp(i),1,nTp(2,i))];
end
perf = [perf,interp_tmp];

% sess 3
tmp = behavior_updated.BCI.Perf.Sess3.Runs{i};
interp_tmp = [];
for i = 1:nRuns
    interp_tmp = [interp_tmp,repmat(tmp(i),1,nTp(3,i))];
end
perf = [perf,interp_tmp];

% ses 4
tmp = behavior_updated.BCI.Perf.Sess4.Runs{i};
interp_tmp = [];
for i = 1:nRuns
    interp_tmp = [interp_tmp,repmat(tmp(i),1,nTp(4,i))];
end
perf = [perf,interp_tmp];
end

