function [ wpli ] = get_window_wpli( data, srate, st, en, freqs, t)
%% Get wpli

% wpli for single trials
% ovlp should be a decimal, 0-1
% st and en are in seconds
% winLength should be ~5*your longest wavelength, in ms

%@author JStiso jeni.stiso@gmail.com

% Change Log
% August 27th, 2018 - created script
% August 29th - updated to calculate wpli over a sliding window
% August 30th - updated to use wavelets


%% get cross spectral density
% The product of the Fourier transform of one signal with the conjugate of the Fourier Transform of the signal gives the Cross-Spectral density (CSD).

% this will be estimated with wavelet coherence, doesn't assume
% stationarity, and gets you lots of time points
% defining constants

cfg = [];
cfg.method     = 'wavelet';                
cfg.width      = 6; 
cfg.output     = 'powandcsd';	
cfg.foi        = freqs;	                
cfg.toi        = st:1/srate:en;
cfg.pad        = 'nextpow2';
cfg.trials     = t;
wave = ft_freqanalysis(cfg, data);

% % get time windows, using 5 windows, as discussed below
% % changed to 3 because we have so little data
% csd_ovlp = .5;
% csd_nWin = 3;
% csd_winLength = floor(winLength*convert/(csd_nWin - ovlp));
% 
nSens = size(data.trial{1},1);
nEdge = (nSens^2 - nSens)/2;
% csd = zeros(nWin, nEdge, numel(freqs));
cnt = 0;
label = zeros(nEdge,2);
for i = 1:nSens
    for j = (i+1):nSens
        cnt = cnt +1;
        %for k = 1:nWin
            % channels should be columns, opposite of what is normally done
            % with brain data. this also matters more for pwelch than cpsd
            label(cnt,:) = [i,j];
            %[curr] = cpsd(squeeze(x(k,i,:))',squeeze(x(k,j,:))',csd_winLength,csd_ovlp,freqs,srate);
            %csd(k,cnt,:) = curr;
        %end
     end
 end



%% Get wPLI
% this is code for getting the debiased wPLI - but with 750 timepoints I
% dont think we need to debias it
% with the exception of the abs at the end, this is the same as the
% fieldtrip function
% The debiased WPLI-square estimator is computed by (i) computing the
% imaginary components of the cross-spectra (Fig. 7),
% X = imag(wave.crsspctrm);
% % (ii) computing the average imaginary component of the cross-spectra
% % (Fig. 8, top), and
% num = zeros(size(X,1),size(X,2));
% for i = 1:size(X,1)
%     for j = (i+1):size(X,1)
%         num = num + squeeze(X(:,:,i)).*squeeze(X(:,:,i));
%     end
% end
% % (iii) normalizing by the computed average over the magnitudes of the
% % imaginary component of the cross-spectra (Fig. 8, bottom).
% denom = zeros(size(num));
% for i = 1:size(X,1)
%     for j = (i+1):size(X,1)
%         denom = denom + abs(squeeze(X(i,:,:)).*squeeze(X(j,:,:)));
%     end
% end
% wpli_v = num./denom;

% wpli - not debiased
% before abs, this gives the same output as fieldtrip
X = imag(wave.crsspctrm);
wpli_v = abs(nanmean(X,3))./nanmean(abs(X),3);

%average over freqs
wpli_v = abs(mean(wpli_v,2));

wpli = zeros(nSens);
for i = 1:nSens
    for j = (i+1):nSens
        row_idx = find(label(:,1)==i);
        col_idx = find(label(:,2)==j);
        idx = intersect(row_idx,col_idx);
        wpli(i,j) = wpli_v(idx);
    end
end

wpli = wpli+wpli';

%% Testing some wavelets

% cfg = [];
% cfg.method     = 'wavelet';                
% cfg.width      = 6; 
% cfg.output     = 'powandcsd';	
% cfg.foi        = 7:1:14;	                
% cfg.toi        = 3:1/srate:6;
% wave = ft_freqanalysis(cfg, data_grad);

% Check that output is the same as fieldtrip - it is
% cfg = [];
% cfg.method = 'wpli';
% wave.dimord = 'rpt_chancmb_freq';
% wave.crsspctrm = permute(wave.crsspctrm, [3,1,2]);
% wave.powspctrm = permute(wave.powspctrm, [3,1,2]);
% ft = ft_connectivityanalysis(cfg,wave);

end

