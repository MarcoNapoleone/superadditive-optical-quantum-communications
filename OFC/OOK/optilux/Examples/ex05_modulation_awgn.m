%   ex05: detection of a linearly modulated digital signal with a receiver.
%       Here we have an oversampling factor of Nt samples, necessary to
%       correctly describe a channel (although here the channel is ideal
%       for the sake of simplicity). We thus move from 1 sample/symbol in
%       the digital domain to Nt samples/symbol in the channel (analog
%       domain), and finally return back to 1 sample/symbol by sampling
%       within rxdsp.

clear
clc

%% Global parameters

Nsymb = 1024;           % number of symbols
Nt = 32;                % number of discrete points per symbol

%% Tx parameters
symbrate = 10;          % symbol rate [Gbaud]
tx.rolloff = 0.01;      % pulse roll-off
modfor = 'qpsk';        % modulation format

%% Channel parameters
SNRdB = 0:1:10;          % signal to noise ratio [dB] over symbrate bandwidth

%% Rx parameters
rx.modformat = modfor;   % modulation format
rx.sync.type = 'da';     % time-recovery method
rx.eftype = 'rootrc';    % electrical filter type
rx.ebw = 0.5;            % electrical filter bandwidth normalized to symbrate
rx.epar = tx.rolloff;    % filter parameters (roll-off in this case)

%==========================================================================

%% Init global variables
Nsamp = Nsymb*Nt;        % overall number of samples
fs = symbrate*Nt;        % sampling rate [GHz]
inigstate(Nsamp,fs);     % initialize global variables: Nsamp and fs.

%% Tx side

rng(1);
pat = pattern2(Nsymb,'rand',struct('format',modfor));

sig = digitalmod(pat,modfor,symbrate,'rootrc',tx);

%% AWGN channel + Rx
for k=1:length(SNRdB)
    sigma2 = 10^(-SNRdB(k)/10)*fs/symbrate;
    % Why "*fs/symbrate",i.e., Nt? Because we have Nt samples per symbol, 
    % hence a bigger bandwidth of a factor Nt and thus more AWGN on the 
    % simulation bandwidth than in the signal bandwidth.
    noise = sqrt(sigma2/2)*(randn(Nsamp,1)+1i*randn(Nsamp,1));
    
    rsig = sig + noise; % AWGN
    
    % rxdsp implements here a baseband receiver
    akhat = rxdsp(rsig,symbrate,pat,rx);      % rx symbols after hard-decision
    SNRdBhat(k) = samp2snr(pat,akhat,modfor); % estimated SNR
end

%% Plot
plot(SNRdB,SNRdBhat,'o',SNRdB,SNRdB,'k--')
grid on
xlabel('SNR [dB]')
ylabel('SNRhat [dB]')
legend('Estimation','Theory')
