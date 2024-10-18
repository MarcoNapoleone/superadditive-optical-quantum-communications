%   ex04: discrete-time additive white Gaussian noise channel,
%       thus at one sample per symbol. Target: estimate the signal-to-noise
%       ratio and check with theory.

clear
clc

%% Global parameters

Nsymb = 1024;           % number of symbols

%% Tx parameters
symbrate = 10;          % symbol rate [Gbaud].
modfor = '16qam';       % modulation format

%% Channel parameters
SNRdB = 0:1:10;         % signal to noise ratio [dB]

%% Init global variables
inigstate(Nsymb);       % initialize global variable

%% Discrete-time AWGN channel

rng(1);
pat = pattern2(Nsymb,'rand',struct('format',modfor));

ak = pat2samp(pat,modfor);
% The average power of a modulated pattern is 1 by default. Here 
% mean(abs(ak).^2) is slightly different to 1 because the sequence is 
% random and of finite length. 
for k=1:length(SNRdB)
    % The noise variance must be referred to the signal power of an
    % infinitely-long signal, which is exactly 1 by default, hence the line
    % below for sigma2:
    sigma2 = 10^(-SNRdB(k)/10); % complex noise variance
    
    nk = sqrt(sigma2/2)*(randn(Nsymb,1)+1i*randn(Nsymb,1)); % noise
    
    rk = ak + nk; % AWGN
    
    SNRdBhat(k) = samp2snr(pat,rk,modfor);
end

%% Plot
plot(SNRdB,SNRdBhat,'o',SNRdB,SNRdB,'k--')
grid on
xlabel('SNR [dB]')
ylabel('SNRhat [dB]')
legend('Estimation','Theory')

% To improve the accuracy increase the number of symbols.