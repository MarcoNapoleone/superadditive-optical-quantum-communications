%   ex07: exhaustive search of the best sampling threshold in on-off keying
%   (OOK). To detect OOK we need to set-up a threshold, which is not midway
%   of the symbols because of the non-coherent detection mode of OOK.
%
%   We focus on a non-optimal OOK, hence without rootrc pulses, but simpler
%   and chepaer to generate in the real world.

clear
clc

%% Global parameters

Nsymb = 1024;           % number of symbols
Nt = 32;                % number of discrete points per symbol

%% Tx parameters
symbrate = 10;          % symbol rate [Gbaud]
tx.rolloff = 0.2;       % pulse roll-off (note: no rootrc pulses here)
tx.emph = 'asin';       % digital-premphasis type
modfor = 'ook';         % modulation format
PdBm = 0;               % power [dBm]       
lam = 1550;             % carrier wavelength [nm]  

%% Channel parameters
SNRdB = 5:15;           % signal to noise ratio [dB] over symbrate bandwidth

%% Rx parameters
rx.modformat = modfor;   % modulation format
rx.sync.type = 'da';     % time-recovery method
rx.oftype = 'gauss';     % optical filter type
rx.obw = Inf;            % optical filter bandwidth normalized to symbrate
rx.eftype = 'gauss';     % optical filter type
rx.ebw = 0.65;           % electrical filter bandwidth normalized to symbrate
rx.epar = tx.rolloff;    % electrical filter extra parameters (not used with 'gauss')
rx.threshold = 0:.1:2;   % threshold (1 is best in coherent detection + AWGN)
%==========================================================================

%% Init global variables
Nsamp = Nsymb*Nt;        % overall number of samples
fs = symbrate*Nt;        % sampling rate [GHz]
inigstate(Nsamp,fs);     % initialize global variables: Nsamp and fs.

%% Tx side

Plin = 10.^(PdBm/10);   % [mW]
E = lasersource(Plin,lam,struct('pol','single'));  % y-pol does not exist

rng(1);
patx = pattern2(Nsymb,'rand',struct('format',modfor));

[sigx, normx] = digitalmod(patx,modfor,symbrate,'costails',tx);

E   = mzmodulator(E, sigx,struct('norm',normx));

Erx = E;
errors = zeros(length(SNRdB),length(rx.threshold));
%% AWGN channel + Rx
for k=1:length(SNRdB)
    sigma2 = 10^(-SNRdB(k)/10)*fs/symbrate;
    noise = sqrt(sigma2/2)*(randn(Nsamp,1)+1i*randn(Nsamp,1));
    
    Erx.field = E.field + noise; % AWGN
    
    rsigx  = rxfrontend(Erx,lam,symbrate,rx); % front-end
    akxhat = rxdsp(rsigx,symbrate,patx,rx);  % basic dsp
    
    % Now we search for the best threshold for our non-coherent detection.
    % We do an exhausitve search since the optimal thershold is unknown.
    for m=1:length(rx.threshold) % search best threshold minimizing errors
        x.threshold = rx.threshold(m);
        patxhat = samp2pat(akxhat,modfor,x);
        errors(k,m) = sum(patx ~= patxhat); % k: SNR. m: threshold
    end
    % The best threshold for a given SNR will be the one minimizing the 
    % number of errors for that SNR.
end

%% Plot
contourf(rx.threshold,SNRdB,errors)
xlabel('threshold')
ylabel('SNR [dB]')
title('number of errors')
colorbar
