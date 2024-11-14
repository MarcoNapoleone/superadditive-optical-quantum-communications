%   ex12: Monte Carlo estimation of the bit error-rate (BER) of a
%         discrete-time Additive white Gaussian noise (AWGN) channel
%         with a quadrature-phase shift keying (QPSK) signal.
%         Like any Monte Carlo estimation, the accuracy can be improved by
%         increasing the number of observations, i.e., symbols in this
%         simulation.

clear
clc

%% Global parameters

Nsymb = 1024;           % number of symbols   ()

%% Tx parameters
symbrate = 10;          % symbol rate [Gbaud]. (bitrate)
modfor = 'qpsk';        % modulation format    (OnOff keying)

%% Channel parameters
SNRdB = 5:1:12;         % signal to noise ratio [dB]        (Rumnore apm ASE noise) > 3 dB, Photodiod reciever rumore shoot, AWGN

%% Rx parameters
rx.type = 'bin';        % get binary patterns

%% Monte Carlo BER parameters

mc.stop = [0.1 68];     % stop at relative error 0.1 with confidence 68%
mc.maxsamp = 1e8;       % maximum number of samples for Monte Carlo

%% Init
inigstate(Nsymb);       % initialize global variable

avgber = zeros(1,length(SNRdB));    % average bit-error rate
stdber = zeros(1,length(SNRdB));    % std of avgber (accuracy of estimation)

%% Discrete-time AWGN channel

rng(1);
[pat, patbin] = pattern2(Nsymb,'rand',struct('format',modfor));

ak = pat2samp(pat,modfor);

for k=1:length(SNRdB)
    sigma2 = 10^(-SNRdB(k)/10); % complex noise variance
    cond = true; 
    while cond % when cond=false the desired accuracy for avgber has been reached
        
        nk = sqrt(sigma2/2)*(randn(Nsymb,1)+1i*randn(Nsymb,1));        
        rk = ak + nk; % AWGN
        
        patbinhat = samp2pat(rk,modfor,rx); % Rx bits
        % now count bit-errors, i.e., estimate the BER
        [cond,avgber(k),nruns,stdber(k)] = meanestimate(patbin ~= patbinhat,mc);
        % display info
        fprintf('SNR = %4.1f [dB]. Bits observed = %6d, <BER> = %.2e, std(BER)/<BER> = %.3f.\n',...
            SNRdB(k),nruns,avgber(k),stdber(k)/avgber(k));
    end
    fprintf('\n')
end

%% Plot

SNRlin = 10.^(SNRdB/10);
bertheory = 0.5*erfc(sqrt(SNRlin/2)); % only for QPSK
semilogy(SNRdB,avgber,'o',SNRdB,bertheory,'k--')
grid on
xlabel('SNR [dB]')
ylabel('BER')
legend('estimation','theory')
