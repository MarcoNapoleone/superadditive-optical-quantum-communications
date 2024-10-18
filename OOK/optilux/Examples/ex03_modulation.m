%   ex03: generating and sampling a linearly modulated digital signal. We
%   use cosine roll-off pulses, which can be sampled directly without the
%   matched filter for the sake of simplicity.

clear
clc

%% Global parameters

Nsymb = 1024;           % number of symbols
Nt = 32;                % number of discrete points per symbol

%% Tx parameters
symbrate = 10;          % symbol rate [Gbaud]
tx.rolloff = 0.01;      % pulse roll-off
modfor = 'qpsk';        % modulation format

%% Init global variables
Nsamp = Nsymb*Nt;       % overall number of samples
fs = symbrate*Nt;       % sampling rate [GHz]
inigstate(Nsamp,fs);    % initialize global variables: Nsamp and fs.

%% Tx side

rng(1);
pat = pattern2(Nsymb,'rand',struct('format',modfor));
% Now create a linearly modulated digital signal
sig = digitalmod(pat,modfor,symbrate,'rc',tx); % rc: raised-cosine pulse

%% Plot
polar(angle(sig(1:Nt:end)),abs(sig(1:Nt:end)),'o')
% we have Nt samples per symbol, hence the sampling 1:Nt:end
