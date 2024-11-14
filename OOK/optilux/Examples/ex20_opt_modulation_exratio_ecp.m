%   ex20: Estimate the eye-closure penalty of on-off keying transmission in
%         presence of Additive white Gaussian noise.
%
%         Note: The Mach-Zehnder modulator here has an extinction ratio, 
%           hence we set the average power by an optical amplifier.

clear
% clc

%% Global parameters

Nsymb = 1024;           % number of symbols
Nt = 32;                % number of discrete points per symbol

%% Tx parameters
symbrate = 10;          % symbol rate [Gbaud].
tx.rolloff = 0.2;       % pulse roll-off
tx.emph = 'asin';       % digital-premphasis type
tx.exratio = 13;        % Mach-Zehnder extinction ratio [dB]
modfor = 'ook';         % modulation format
PdBm = 0;               % power [dBm]
lam = 1550;             % carrier wavelength [nm]

%% Channel parameters
SNRdB = -10:1:30;         % signal to noise ratio [dB] over symbrate bandwidth

%% Rx parameters
rx.modformat = modfor;      % modulation format
rx.sync.type = 'da';        % time-recovery method
rx.oftype = 'gauss';        % optical filter type
rx.obw = 1.8;               % optical filter bandwidth normalized to symbrate

%% Init
Nsamp = Nsymb*Nt;        % overall number of samples
fs = symbrate*Nt;        % sampling rate [GHz]
inigstate(Nsamp,fs);     % initialize global variables: Nsamp and fs.

eyeop = zeros(1,length(SNRdB));    % average bit-error rate
ecp   = zeros(1,length(SNRdB));    % std of avgber (accuracy of estimation)

%% Tx side

Plin = 10.^(PdBm/10);   % [mW]
E = lasersource(Plin,lam);

rng(1);
[patx, patbinx] = pattern2(Nsymb,'rand',struct('format',modfor));

sigx = digitalmod(patbinx,modfor,symbrate,'costails',tx);

E   = mzmodulator(E, sigx);
E   = ampliflat(E,struct('outpower',PdBm)); % set time-averaged power

Pest = mean(sum(abs(E.field).^2)); % average power [mW]
fprintf('Time-average power = %.4f [dBm]\n',10*log10(Pest))

% Now evaluate eye opening in back-to-back (b2b). It is necessary to get
% the eye-closure penalty later.
rsig = rxfrontend(E,lam,symbrate,rx);
eyeopb2b = evaleye(patx,rsig,symbrate,modfor);

% Now the eye after AWGN channel
Er = E; % put aside
for k=1:length(SNRdB)
    sigma2 = 10^(-SNRdB(k)/10)*Plin*fs/symbrate;
    
    noise = sqrt(sigma2/2)*(randn(Nsamp,1)+1i*randn(Nsamp,1));
    Er.field = E.field + noise; % AWGN
    
    rsig = rxfrontend(Er,lam,symbrate,rx);      % front-end
    
    eyeop(k) = evaleye(patx,rsig,symbrate,modfor); % eye-opening [dB]
    ecp(k) = eyeopb2b - eyeop(k); % eye-closure penalty [dB]
    
    %     % display info
    %     fprintf('SNR = %4.1f [dB]. ECP = %.2f [dB].\n',...
    %         SNRdB(k),ecp(k));
    %     fprintf('\n')
end

%% Plot

SNRlin = 10.^(SNRdB/10);
plot(SNRdB,ecp)
grid on
xlabel('SNR [dB]')
ylabel('ECP [dB]')
title('with AWGN')
