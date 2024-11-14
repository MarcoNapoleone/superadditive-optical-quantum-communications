%   ex19: Optical propagation within a dispersion-managed (DM) optical link
%   in nonlinear regime.
%
%    -        _____                          \ x Nspan
%   /        /     \                 |\       |  
%   |        | ft  |                 | \      |  
%   |        \     /                 |  \     |
%  -|--------------------------------|   \ --------
%   |                      /    \    |   /    |      
%   |                      | fc |    |  /     |  
%   |                      \____/    | /      |  
%    \                               |/      /   
%     -
%
%   There are many fiber parameters untouched here but set by fiber.m to 
%   default values, such as:
%
%   ft.coupling = 'pol' (alternatively can be 'none')
%   ft.dphi1max = 20; (max FWM phase shift per step)
%   ft.stepupd = 'cle' (alternatively can be the nonlinear phase criterion
%                       'nlp')
%   ft.nplates = 1 (1 because here the fiber is solved by the Manakov
%           equation without PMD. Otherwise, the default is ft.nplates=100)

clear
clc

%% Global parameters

Nsymb = 1024;           % number of symbols
Nt = 32;                % number of discrete points per symbol

%% Tx parameters
symbrate = 10;          % symbol rate [Gbaud].
tx.rolloff = 0.2;       % pulse roll-off
tx.emph = 'asin';       % digital-premphasis type
modfor = 'ook';         % modulation format
PdBm = 3;               % power [dBm]
lam = 1550;             % carrier wavelength [nm]

%% Channel parameters

RDPS  = 0;                  % residual dispersion per span [ps/nm]
Nspan = 5;                  % number of spans

% Transmission fiber
ft.length     = 100E3;      % length [m]
ft.lambda     = 1550;       % wavelength [nm] of fiber parameters
ft.alphadB    = 0.2;        % attenuation [dB/km]
ft.disp       = 17;         % dispersion [ps/nm/km] @ ft.lambda
ft.slope      = 0.057;      % slope [ps/nm^2/km] @ ft.lambda
ft.pmdpar     = 0;          % PMD parameter [ps/sqrt(km)]
ft.ismanakov  = true;       % Solve Manakov equation
ft.aeff       = 80;         % effective area [um^2]
ft.n2         = 2.5e-20;    % nonlinear index [m^2/W]
% ft.dphimax    = 20;       % maximum nonlinear phase per step [rad]
ft.dzmax      = 2E4;        % maximum SSFM step size [m]
ft.trace      = true;       % show information on screen

% compensating fiber
fc = ft;                    % same parameters but:
fc.length     = 1e3;        % fixed length [m]
fc.alphadB    = 0.6;        % attenuation [dB/km]
fc.disp       = (RDPS*1e3 - ft.disp*ft.length)/fc.length; % [ps/nm/km] to get RDPS
fc.slope      = 0.057;      % slope [ps/nm^2/km] @ fc.lambda
fc.pmdpar     = 0;          % PMD parameter [ps/sqrt(km)]
fc.aeff       = 20;         % effective area [um^2]
fc.dzmax      = 2E4;        % maximum step size [m]

% Optical amplifier
amp.gain = (ft.length*ft.alphadB + fc.length*fc.alphadB)*1e-3; % gain [dB]

%% Rx parameters
rx.modformat = modfor;      % modulation format
rx.sync.type = 'da';        % time-recovery method
rx.oftype = 'gauss';        % optical filter type
rx.obw = Inf;               % optical filter bandwidth normalized to symbrate
rx.eftype = 'rootrc';       % optical filter type
rx.ebw = 0.5;               % electrical filter bandwidth normalized to symbrate
rx.epar = tx.rolloff;
rx.type = 'bin';            % binary pattern

%% Init
Nsamp = Nsymb*Nt;           % overall number of samples
fs = symbrate*Nt;           % sampling rate [GHz]
inigstate(Nsamp,fs);        % initialize global variables: Nsamp and fs.

%% Tx side

Plin = 10.^(PdBm/10);   % [mW]
E = lasersource(Plin,lam);  

rng(1);
[patx, patbinx] = pattern2(Nsymb,'rand',struct('format',modfor));

[sigx, normx] = digitalmod(patbinx,modfor,symbrate,'costails',tx);

E   = iqmodulator(E, sigx,struct('norm',normx));

%% Channel

plotfield(E,'p---',struct('pol','x'))

for ns=1:Nspan
    [E,parft] = fiber(E,ft);
    [E,parfc] = fiber(E,fc);
    E = ampliflat(E,amp);
end
plotfield(E,'p---','r',struct('pol','x'))

legend('Tx','Rx')