%   ex13: Monte Carlo estimation of the bit error-rate of a
%         modulated signal in Additive white Gaussian noise.

clear
clc

%% Global parameters

Nsymb = 1024;            % number of symbols
Nt = 32;                 % number of discrete points per symbol

%% Tx parameters
symbrate = 10;           % symbol rate [Gbaud].
tx.rolloff = 0.01;       % pulse roll-off
modfor = 'qpsk';         % modulation format

%% Channel parameters
SNRdB = 5:1:10;          % signal to noise ratio [dB] over symbrate bandwidth

%% Rx parameters
rx.modformat = modfor;   % modulation format
rx.sync.type = 'da';     % time-recovery method
rx.oftype = 'gauss';     % optical filter type
rx.obw = Inf;            % optical filter bandwidth normalized to symbrate
rx.eftype = 'rootrc';    % optical filter type
rx.ebw = 0.5;            % electrical filter bandwidth normalized to symbrate
rx.epar = tx.rolloff;    % electrical filter extra parameters
rx.type = 'bin';	     % binary pattern

%% Monte Carlo BER parameters

mc.stop = [0.1 68];      % stop at relative error 0.1 with confidence 68%
mc.maxsamp = 1e8;        % maximum number of samples for Monte Carlo

%% Init
Nsamp = Nsymb*Nt;        % overall number of samples
fs = symbrate*Nt;        % sampling rate [GHz]
inigstate(Nsamp,fs);     % initialize global variables: Nsamp and fs.

avgber = zeros(1,length(SNRdB));    % average bit-error rate
stdber = zeros(1,length(SNRdB));    % std of avgber (accuracy of estimation)

%% Discrete-time AWGN channel

rng(1);
[pat, patbin] = pattern2(Nsymb,'rand',struct('format',modfor));

sig = digitalmod(patbin,modfor,symbrate,'rootrc',tx);

for k=1:length(SNRdB)
    sigma2 = 10^(-SNRdB(k)/10)*fs/symbrate;
    cond = true; 
    while cond % when cond=false the desired accuracy for avgber has been reached
        
        noise = sqrt(sigma2/2)*(randn(Nsamp,1)+1i*randn(Nsamp,1));
        rsig = sig + noise; % AWGN
        
        akhat = rxdsp(rsig,symbrate,patbin,rx);
        patbinhat = samp2pat(akhat,modfor,rx); % Rx bits
        % now count bit-errors
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
