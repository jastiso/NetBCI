function [ wpli ] = get_window_wpli( data, srate, st, en, freqs )
%% Get wpli

% wpli for single trials

%@author JStiso jeni.stiso@gmail.com

% Change Log
% August 27th, 2018 - created script

%% Housekeeping 
% defining constants

% get time windows, using 5 windows, as discussed below
ovlp = .5;
nWin = 5;
samples = data.time{1}((st*srate):(en*srate));
winLength = floor(numel(samples)/(nWin - ovlp));

%% Fourier spectra with welch's method
% "Regarding wPLI and ImC, the developed function follows Fieldtrip?s 
% implementation for wPLI. In it, Fourier spectra is first computed for 
% each signal through a Welch?s method with ?5 windows (The length of 
% each signal is divided by 4.5 and the result is transformed into an 
% integer number of windows) with 50% 
% overlapping between them."  - ?Garc�a-Prieto et al 2017
%
% It seems like welch's methos just means splitting soething up into
% windows

% get indices of each window
% trl = st:(winLength*ovlp):(en-winLength);
% trl = [trl; trl + (winLength-1)];
% trl = [trl; zeros(1,size(trl,2))]';
% cfg = [];
% cfg.trl = trl;
% data = ft_redefinetrial(cfg,data);
% 
% cfg = [];
% cfg.method = 'mtmfft';
% cfg.taper = 'hanning';
% cfg.keeptrials = 'yes';
% cfg.foilim = freqs;
% pow = ft_freqanalysis(cfg,data);

%x = data.trial{1}'; % channels should be columns, opposite of what is normally done
%nfft = 2^nextpow2(size(x,1));
%[psd,f] = pwelch(x,winLength,ovlp,nfft,srate);

%% combine gradiometers
% not sure what step this should happen in....

%psd_cmb = combine_gradiometers(psd,0); % 0 flags this as PSD

%% get cross spectral density
% The product of the Fourier transform of one signal with the conjugate of the Fourier Transform of the signal gives the Cross-Spectral density (CSD). 

x = data.trial{1}(:,(st*srate):(en*srate))'; % channels should be columns, opposite of what is normally done
nSens = size(x,2);
%nfft = 2^nextpow2(size(x,1));
csd = cell(nSens,nSens);
for i = 1:nSens
    for j = i:nSens
        [curr] = cpsd(x(:,i),x(:,j),winLength,ovlp,freqs,srate);
        csd{i,j} = curr;
    end
end


%% Get wPLI

wpli = zeros(nSens,nSens);
for i = 1:nSens
    for j = (i+1):nSens
        curr = csd{i,j};
        wpli(i,j) = abs(nanmean(imag(curr)))/nanmean(abs(imag(curr)));
    end
end
wpli = wpli + wpli';


end

