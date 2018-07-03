function [ wpli ] = get_window_wpli( data, srate, freqs )
% Get wpli in sliding windows

nElec = size(numel(data.labels));
nFreq = size(freqs,1);

win = 751%3*srate;
ovlp = 0;%floor(win);
curr = chunk_mag.trial{1};
wpli = zeros(nElec, nElec, nFreq);
for i = 1:nElec
    for j = 1:nElec
        [Pxy, f] = cpsd(curr(i,:), curr(j,:), win, ovlp, freqs(1,1):freqs(end,end), srate);
        imC = imag(Pxy);
        for k = 1:nFreq
            f_idx = f >= freqs(k,1) & f < freqs(k,2);
            wpli(i,j,k) = abs(mean(imC(f_idx)))/mean(abs(imC(f_idx)));
        end
    end
end

for i = 1:nFreq
    imagesc(squeeze(wpli(:,:,i)))
    colorbar
    pause(1)
end
end

