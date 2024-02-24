function [freq,modulus] = WaveletTransform(time,signal,minFreq, ...
    maxFreq,options)
% Help truncated in this screenshot
% Citation
% --------
% K.J. Moore, M. Kurt, M. Eriten, D.M. McFarland, L.A. Bergman, 
% A.F. Vakakis, "Wavelet-Bounded Empirical Mode Decomposition for
% Measured Time Series Analysis," Mechanical Systems and Signal 
% Processing, 99:14–29, 2018.
% https://dx.doi.org/10.1016/j.ymssp.2017.06.005

arguments
    time (:,1) double
    signal (:,1) double
    minFreq double
    maxFreq double
    options.numFreq double = 100;
    options.motherWaveletFreq double = 2;
end

% Transform Parameters
dt = time(2)-time(1);
lengthSignal = length(signal);
df = (maxFreq-minFreq)/options.numFreq; % Frequency Resolution
freq = minFreq:df:maxFreq; % Prescribe the frequencies of interest
waveletScale = options.motherWaveletFreq./freq; % Compute scales

% FFT Parameters
nfourier = 2^nextpow2(lengthSignal); % Zero-filling
npt = nfourier/2;
Fourierfreq = 1/dt*(0:npt-1)/nfourier; % Frequency vector

% Compute FFT of xnew
signalFFT = fft(signal,nfourier);
signalFFT(npt+1:end) = [];

% Vectorized Computation of Wavelet Transform of signal
core2 = bsxfun(@times,conj(bsxfun(@times,(2^0.5)*exp(-0.5* ...
    (2*pi*(bsxfun(@times,Fourierfreq',waveletScale)- ...
    options.motherWaveletFreq)).^2),sqrt(waveletScale))), ...
    signalFFT);

% Assert Admissibility Condition
if minFreq == 0
    core2(:,1) = 0;
end
result = ifft(core2,nfourier);
result1 = result(1:lengthSignal,:);
modulus  = abs(result1);
end