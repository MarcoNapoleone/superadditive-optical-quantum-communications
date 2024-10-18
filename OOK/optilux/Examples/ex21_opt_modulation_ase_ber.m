%   ex21: Monte Carlo estimation of the bit-error rate of a
%   polarization-division multiplexed signal propagating in an optical link
%   with only amplified-spontaneous emission (ASE) noise.
%
%   Note: you will read on the screen the warning "Too many samples..."
%       for the first spans, because with a small ASE the error
%       probability is close to zero. For the same reason, it is normal to
%       observe a bad accuracy in the final result because, for the sake of
%       speed, the estimation is stopped after observing just 20000
%       samples.

clear

global CONSTANTS
clc

%% Global parameters

Nsymb = 1024;           % number of symbols
Nt = 32;                % number of discrete points per symbol

%% Tx parameters
symbrate = 10;          % symbol rate [Gbaud].
tx.rolloff = 0.01;      % pulse roll-off
tx.emph = 'asin';       % digital-premphasis type
modfor = 'qpsk';        % modulation format
PdBm = 0;               % power [dBm]
lam = 1550;             % carrier wavelength [nm]

%% Channel parameters

Nspan = 50;                % number of spans

% Transmission fiber (a pure lossy element)
ft.length     = 150E3;      % length [m]
ft.lambda     = 1550;       % wavelength [nm] of fiber parameters
ft.alphadB    = 0.2;        % attenuation [dB/km]
ft.disp       = 0;          % dispersion [ps/nm/km] @ ft.lambda
ft.slope      = 0;          % slope [ps/nm^2/km] @ ft.lambda
ft.pmdpar     = 0;          % PMD parameter [ps/sqrt(km)]
ft.aeff       = 80;         % effective area [um^2]
ft.n2         = 0;          % nonlinear index [m^2/W]

% Optical amplifier
amp.gain = ft.length*ft.alphadB*1e-3;   % gain [dB]
amp.f = 6;                              % noise figure [dB]

%% Rx parameters
rx.modformat = modfor;      % modulation format
rx.sync.type = 'da';        % time-recovery method
rx.oftype = 'gauss';        % optical filter type
rx.obw = Inf;               % optical filter bandwidth normalized to symbrate
rx.eftype = 'rootrc';       % optical filter type
rx.ebw = 0.5;               % electrical filter bandwidth normalized to symbrate
rx.epar = tx.rolloff;
rx.type = 'bin';            % binary pattern

%% Monte Carlo BER parameters

% We voluntarily stop earlier to speed up the simulation. As a result, a
% WARNING of low accuracy will appear on the screen after the first spans,
% for which the bit-error-rate is expected to be very low. After the last
% spans the accumulated noise is so large that the bit-error rate is
% estimated accurately within the maximum limit of mc.maxsamp.

mc.stop = [0.1 68];     % stop at relative error 0.1 with confidence 68%
mc.maxsamp = 2e4;       % maximum number of samples for Monte Carlo

%% Init
Nsamp = Nsymb*Nt;        % overall number of samples
fs = symbrate*Nt;        % sampling rate [GHz]
inigstate(Nsamp,fs);     % initialize global variables: Nsamp and fs.

avgber = zeros(1,Nspan);    % average bit-error rate
stdber = zeros(1,Nspan);    % std of avgber (accuracy of estimation)

%% Tx side

Plin = 10.^(PdBm/10);   % [mW]
E = lasersource(Plin,lam);  % y-pol does not exist

[Ex,Ey] = pbs(E); % split in two orthogonal polarizations

rng(1);
[patx, patbinx] = pattern2(Nsymb,'rand',struct('format',modfor));
[paty, patbiny] = pattern2(Nsymb,'rand',struct('format',modfor));

[elecx, normx] = digitalmod(patx,modfor,symbrate,'rootrc',tx);
[elecy, normy] = digitalmod(paty,modfor,symbrate,'rootrc',tx);

Ex   = iqmodulator(Ex, elecx,struct('norm',normx));
Ey   = iqmodulator(Ey, elecy,struct('norm',normy));

Esp = pbc(Ex,Ey); % combine creating a PDM signal

patbin = [patbinx patbiny]; % concatenate the patterns
for ns=1:Nspan
    cond = true;
    while cond % when cond=false the desired accuracy for avgber has been reached
        
        E = Esp; % get field outgoing the previous span
        E = fiber(E,ft); % iterate over the copy E, such that Esp is untouched
        E = ampliflat(E,amp);
        
        rsig = rxfrontend(E,lam,symbrate,rx);      % front-end
        
        % please note: patbinx in rxdsp!
        akhat = rxdsp(rsig,symbrate,[patx paty],rx);
        patbinhat = samp2pat(akhat,modfor,rx); % Rx bits
        SNRdBhat(ns) = samp2snr(patbin,akhat,modfor);
        
        % now count bit-errors
        [cond,avgber(ns),nruns,stdber(ns)] = meanestimate(patbin ~= patbinhat,mc);
        % display info
        fprintf('span = %d/%d. <BER> = %.2e, std(BER)/<BER> = %.3f.\n',...
            ns,Nspan,avgber(ns),stdber(ns)/avgber(ns));
    end
    Esp = E;
    fprintf('\n')
end

%% Comparison with theory

% Theoretical SNR and BER

CLIGHT  = CONSTANTS.CLIGHT;  % speed of light [m/s]
HPLANCK = CONSTANTS.HPLANCK; % Planck's constant [J*s]

Fac58 = 10*log10(HPLANCK*CLIGHT/E.lambda*1e9*12.5e9)+30; % The famous -58 dBm @ Bo=12.5 GHz
N0B = Fac58 + amp.gain + amp.f + 10*log10(1:Nspan) + 10*log10(symbrate/12.5);
SNRdB = PdBm - N0B; % theoretical SNR [dB] @ bandwidth=symbrate.

SNRlin = 10.^(SNRdB/10);
bertheory = 0.5*erfc(sqrt(SNRlin/2)); % only for QPSK

distance = (1:Nspan)*ft.length*1e-3; % propagated distance [km]
semilogx(distance,ber2q(avgber),'o',distance,ber2q(bertheory),'k--')
grid on
xlabel('distance [km]')
ylabel('Q-factor [dB]')
legend('estimation','theory')
title(sprintf('estimation stopped after %d samples',mc.maxsamp))
