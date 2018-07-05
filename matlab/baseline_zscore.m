function [ data ] = baseline_zscore( data, bl )
% get data as z score of bl power
% data and bl are the output of fieldtrip ft_freqalaysis function - they
% need to have a 4D .powspctrm field where the fourth dimension is time


%   takes fieldtrip data as input
nBl = size(bl.powspctrm,1);
nTrial = size(data.powspctrm,1);
    
bl_mu = mean(bl.powspctrm,4);
bl_sd = std(bl.powspctrm,0,4);
data.powspctrm = (data.powspctrm - bl_mu)./bl_sd;

    
end
