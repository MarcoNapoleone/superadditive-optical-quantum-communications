%   ex10: WDM transmission, method 2/2
%   
%   Note: the "0" frequency is by convention in the barycenter of the
%       spectrum. Since we call the  multiplexer (Nch-1) times, each time 
%       multiplexing only two signals,the "0" is Nch-dependent. 
%       For instance (assuming frequency spacing equal to 1):
%
%   WDM at freq: call0 (start): 0   1   2   3   4
%                call1        :  0.5    2   3   4
%                call2        :     1.25    3   4
%                call3        :         2.125   4
%                call4        :            3.0625
%
%   After call4 the spectrum is centered at 3.0625 which is set equal to 0.
%   This is not a problem since the position of the zero frequency is just
%   a convention. In absolute terms, the zero frequency corresponds to the 
%   wavelength Em.lambda.
%
%   A closer look to the main figure of this example reveals that the last
%   channel in the frequency domain (first in wavelength) is spaced 
%   deltaf*3.0625 to the zeroth frequency.
%
%   By using the call
%
%       Em = multiplexer(Em,E,struct('lambda',1550));
%
%   the zero frequency is forced to corresponds to the wavelength 1550 nm.

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

lamv = getlambda(lam,spac,Nch); % carrier wavelengths [nm]

for k=1:Nch
    rng(1e5+k-ceil(Nch/2)); % ceil(..): keep seed for increasing Nch    
    E = lasersource(Plin(k),lamv(k));
    
    [Ex,Ey] = pbs(E); % split in two orthogonal polarizations
    patx = pattern2(Nsymb,'rand',struct('format',modfor));
    paty = pattern2(Nsymb,'rand',struct('format',modfor));
    
    [elecx, normx] = digitalmod(patx,modfor,symbrate,'rc',tx);
    [elecy, normy] = digitalmod(paty,modfor,symbrate,'rc',tx);

    Ex  = iqmodulator(Ex, elecx,struct('norm',normx));
    Ey  = iqmodulator(Ey, elecy,struct('norm',normy));
    
    E = pbc(Ex,Ey); % combine creating a PDM signal
    
    if k == 1
        Em = E; % nothing to multiplex together
    else
        Em = multiplexer(Em,E);
%         Em = multiplexer(Em,E,struct('lambda',lam)); % try to force 0th freq
    end
end

%% Plot

plot(fftshift(GSTATE.FN),20*log10(abs(fftshift(fft(Em.field)))))
grid on
xlabel('frequency [GHz]')
ylabel('PSD [dB]')

% the tilted spectrum is the result of the tilted powers.