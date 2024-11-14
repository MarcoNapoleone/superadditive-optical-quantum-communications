%   ex06: generation and detection of a single-polarization electric field, 
%        by default on the x-polarization. Here we thus focus on a optical
%        signal and its peculiarities (the laser, Mach Zehnder modulation, 
%        etc).

clear
clc

%% Global parameters

Nsymb = 1024;           % number of symbols
Nt = 32;                % number of discrete points per symbol

%% Tx parameters
symbrate = 10;          % symbol rate [Gbaud]
tx.rolloff = 0.01;      % pulse roll-off
tx.emph = 'asin';       % digital-premphasis type
modfor = 'bpsk';        % modulation format
PdBm = 0;               % power [dBm]       
lam = 1550;             % carrier wavelength [nm]  

%% Rx parameters
rx.modformat = modfor;   % modulation format
rx.sync.type = 'da';     % time-recovery method
rx.oftype = 'gauss';     % optical filter type
rx.obw = Inf;            % optical filter bandwidth normalized to symbrate
rx.eftype = 'rootrc';    % optical filter type
rx.ebw = 0.5;            % electrical filter bandwidth normalized to symbrate
rx.epar = tx.rolloff;

%==========================================================================

%% Init global variables
Nsamp = Nsymb*Nt;        % overall number of samples
fs = symbrate*Nt;        % sampling rate [GHz]
inigstate(Nsamp,fs);     % initialize global variables: Nsamp and fs.

%% Tx side

% We first generate a constant wave laser. We just want
% one-polarization (x-polarization) for the sake of simplicity, hence we 
% add the option 'single' to lasersource (default is dual-polarization).
Plin = 10.^(PdBm/10);   % [mW]
E = lasersource(Plin,lam,struct('pol','single'));  % y-pol does not exist

rng(1);
patx = pattern2(Nsymb,'rand',struct('format',modfor));

[sigx,normx] = digitalmod(patx,modfor,symbrate,'rootrc',tx);

% modulate the BPSK signal with a Mach-Zehnder modulator
E   = mzmodulator(E, sigx,struct('norm',normx));

%% Rx side

% We first select the channel by returning back in base-band through the
% front-end, then we perform some basic digital-signal processing to 
% extract the symbols.

rsigx = rxfrontend(E,lam,symbrate,rx);   % front-end
akxhat = rxdsp(rsigx,symbrate,patx,rx);  % basic dsp
patxhat = samp2pat(akxhat,modfor);

%% Measurements
polar(angle(akxhat),abs(akxhat),'o')

errors = sum(patx ~= patxhat);
fprintf('PdBm: %.2f [dBm]. Pmes: %.2f [dBm]\n',...
    PdBm,10*log10(mean(abs(E.field).^2)))
fprintf('\nErrors: %d\n\n',errors)

