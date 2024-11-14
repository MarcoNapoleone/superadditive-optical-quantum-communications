%   ex09: WDM transmission, method 1/2
%   
%   Note: the "0" frequency is by convention in the barycenter of the
%       spectrum. Since we call the multiplexer just once, the "0" is in 
%       the barycenter of the WDM.

clear
clc

%% Global parameters
global GSTATE

Nsymb = 1024;           % number of symbols
Nt = 32;                % number of discrete points per symbol

%% Tx parameters
Nch = 5;                % number of channels
symbrate = 10;          % symbol rate [Gbaud]
tx.rolloff = 0.01;      % pulse roll-off
tx.emph = 'asin';       % digital-premphasis type
modfor = 'qpsk';        % modulation format
PdBm = 0:5:5*(Nch-1);   % power [dBm], tilted for the sake of visualization
lam = 1550;             % carrier wavelength [nm]  
deltaf = 37.5;          % channel spacing [GHz]

%% Conversions
spac = deltaf*lam^2/299792458;  % channel spacing [nm]
Plin = 10.^(PdBm/10);   % [mW]

%% Init global variables
Nsamp = Nsymb*Nt;        % overall number of samples
fs = symbrate*Nt;        % sampling rate [GHz]
inigstate(Nsamp,fs);     % initialize global variables: Nsamp and fs.

%% Tx side

E = lasersource(Plin,lam,spac,Nch);  % electric field

[Ex,Ey] = pbs(E); % split in two orthogonal polarizations

for k=1:Nch
    rng(1e5+k-ceil(Nch/2)); % ceil(..): keep seed for increasing Nch
    patx = pattern2(Nsymb,'rand',struct('format',modfor));
    paty = pattern2(Nsymb,'rand',struct('format',modfor));
    
    [elecx, normx] = digitalmod(patx,modfor,symbrate,'rc',tx);
    [elecy, normy] = digitalmod(paty,modfor,symbrate,'rc',tx);
    % modulate carrier k of the matrix Ex.field
    Ex   = iqmodulator(Ex, elecx,struct('nch',k,'norm',normx));
    Ey   = iqmodulator(Ey, elecy,struct('nch',k,'norm',normy));
end

E = pbc(Ex,Ey); % combine creating a PDM signal

Em = multiplexer(E); % E.field: 2 polarization * Nch carriers: 2*Nch columns
                     % Em.field: 2 columns (polarizations)

%% Plot

plot(fftshift(GSTATE.FN),20*log10(abs(fftshift(fft(Em.field)))))
grid on
xlabel('frequency [GHz]')
ylabel('PSD [dB]')

% the tilted spectrum is the result of the tilted powers.