



%%% Attenuation %%%

alpha = 20e-3;
distance = 10e5;   % 100 km
eta = exp(-alpha*distance);



% Parameters
M = 32; % Fixed M value for Hadamard/Fourier Machine

% Array of mean photon numbers
n_R = [0.01,0.1, 1]

%% Channel parameters
SNRdB = 0:1:20; % SNR in dB array
SNR = 10 .^ (SNRdB / 10); % Convert SNR from dB to linear scale



