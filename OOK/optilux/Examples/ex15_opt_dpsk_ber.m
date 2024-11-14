%   ex15: Monte Carlo estimation of the bit-error rate of a
%         differential-phase shift keying (DPSK) in Additive
%         white Gaussian noise.

clear
clc

%% Global parameters

Nsymb = 1024;           % number of symbols
Nt = 32;                % number of discrete points per symbol

%% Tx parameters
symbrate = 10;          % symbol rate [Gbaud].
tx.rolloff = 0.01;      % pulse roll-off
tx.emph = 'asin';       % digital-premphasis type
modfor = 'dpsk';        % modulation format
PdBm = 0;               % power [dBm]
lam = 1550;             % carrier wavelength [nm]

%% Channel parameters
SNRdB = 3:1:8;          % signal to noise ratio [dB] over symbrate bandwidth

%% Rx parameters
rx.modformat = modfor;      % modulation format
rx.sync.type = 'da';        % time-recovery method
rx.oftype = 'movavg';       % optical filter type
rx.obw = 2;                 % optical filter bandwidth normalized to symbrate
rx.eftype = 'gauss';        % optical filter type
rx.ebw = Inf;               % electrical filter bandwidth normalized to symbrate
rx.type = 'bin';            % binary pattern

%% Monte Carlo BER parameters

mc.stop = [0.1 95];      % stop at relative error 0.1 with confidence 68%
mc.maxsamp = 1e8;        % maximum number of samples for Monte Carlo

%% Init
Nsamp = Nsymb*Nt;        % overall number of samples
fs = symbrate*Nt;        % sampling rate [GHz]
inigstate(Nsamp,fs);     % initialize global variables: Nsamp and fs.

avgber = zeros(1,length(SNRdB));    % average bit-error rate
stdber = zeros(1,length(SNRdB));    % std of avgber (accuracy of estimation)

%% Tx side

Plin = 10.^(PdBm/10);   % [mW]
E = lasersource(Plin,lam,struct('pol','single'));  % y-pol does not exist

rng(1);
[patx, patbinx] = pattern2(Nsymb,'rand',struct('format',modfor));
% differential encoder for the (D)PSK
patbinex = diffencoder(patbinx,modfor);

[sigx, normx] = digitalmod(patbinex,modfor,symbrate,'costails',tx);

E   = mzmodulator(E, sigx,struct('norm',normx));

Er = E;
for k=1:length(SNRdB)
    sigma2 = 10^(-SNRdB(k)/10)*fs/symbrate;
    cond = true;
    while cond % when cond=false the desired accuracy for avgber has been reached
        
        noise = sqrt(sigma2/2)*(randn(Nsamp,1)+1i*randn(Nsamp,1));
        Er.field = E.field + noise; % AWGN
        
        rsigx = rxfrontend(Er,lam,symbrate,rx);      % front-end
        
        % please note: patbinx in rxdsp!
        akxhat = rxdsp(rsigx,symbrate,patbinx,rx);
        patbinxhat = samp2pat(akxhat,modfor,rx); % Rx bits
        
        % now count bit-errors (reject first bit: differential decoder status)
        [cond,avgber(k),nruns,stdber(k)] = ...
            meanestimate(patbinx(2:end,:) ~= patbinxhat(2:end,:),mc);
        % display info
        fprintf('SNR = %4.1f [dB]. Bits observed = %6d, <BER> = %.2e, std(BER)/<BER> = %.3f.\n',...
            SNRdB(k),nruns,avgber(k),stdber(k)/avgber(k));
    end
    fprintf('\n')
end

%% Plot

SNRlin = 10.^(SNRdB/10);
bertheory = 0.5*exp(-SNRlin); % only for DPSK
semilogy(SNRdB,avgber,'o',SNRdB,bertheory,'k--')
grid on
xlabel('SNR [dB]')
ylabel('BER')
legend('estimation','theory')

% The accuracy at high SNRs can be increased by increasing the number of
% symbols