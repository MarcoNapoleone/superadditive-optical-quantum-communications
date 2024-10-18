%   ex08: generation of a polarization division multiplexing (PDM) signal.
%   We still use rc pulses for the sake of simplicity.

clear
clc

%% Global parameters

Nsymb = 1024;           % number of symbols
Nt = 32;                % number of discrete points per symbol

%% Tx parameters
symbrate = 10;          % symbol rate [Gbaud]
tx.rolloff = 0.01;      % pulse roll-off
tx.emph = 'asin';       % digital-premphasis type
modfor = 'qpsk';        % modulation format
PdBm = 0;               % power [dBm]       
lam = 1550;             % carrier wavelength [nm]  

%% Init global variables
Nsamp = Nsymb*Nt;        % overall number of samples
fs = symbrate*Nt;        % sampling rate [GHz]
inigstate(Nsamp,fs);     % initialize global variables: Nsamp and fs.

%% Tx side

Plin = 10.^(PdBm/10);   % [mW]
E = lasersource(Plin,lam);  % electric field

[Ex,Ey] = pbs(E); % split in two orthogonal polarizations

rng(1);
% now separately modulate the two polarizations
patx = pattern2(Nsymb,'rand',struct('format',modfor));
paty = pattern2(Nsymb,'rand',struct('format',modfor));

[elecx, normx] = digitalmod(patx,modfor,symbrate,'rc',tx);
[elecy, normy] = digitalmod(paty,modfor,symbrate,'rc',tx);

Ex   = iqmodulator(Ex, elecx,struct('norm',normx));
Ey   = iqmodulator(Ey, elecy,struct('norm',normy));

E = pbc(Ex,Ey); % combine the two polarizations creating a PDM signal

%% Plot
polar(angle(E.field(1:Nt:end,1)),abs(E.field(1:Nt:end,1)),'o')
grid on, hold on
polar(angle(E.field(1:Nt:end,2)),abs(E.field(1:Nt:end,2)),'*')
legend('x','y')
