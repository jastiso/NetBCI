function [ data ] = baseline_correct( data, bl )
% get data as z score of bl power

%   takes fieldtrip data as input
nBl = size(bl.powspctrm,1);
nTrial = size(data.powspctrm,1);
if nBl == nTrial
    bl_mu = mean(bl.powspctrm,4);
    bl_sd = std(bl.powspctrm,0,4);
    data.powspctrm = (data.powspctrm - bl_mu)./bl_sd;
elseif nBl < nTrial
    bl_mu = mean(bl.powspctrm,3);
    bl_sd = std(bl.powspctrm,0,3);
    data.powspctrm = (data.powspctrm(1:nBl,:,:) - bl_mu)./bl_sd;
elseif nBl > nTrial
    bl_mu = mean(bl.powspctrm(1:nTrial,:,:),3);
    bl_sd = std(bl.powspctrm(1:nTrial,:,:),0,3);
    data.powspctrm = (data.powspctrm - bl_mu)./bl_sd;
end
    
end

