%   ex02: discrete-time channel, thus at one sample per symbol.

clear
clc

%% Global parameters

Nsymb = 1024;           % number of symbols

%% Tx parameters
symbrate = 10;          % symbol rate [Gbaud].
modfor = '16qam';       % modulation format

%% Discrete-time channel

% generate random symbols in modfor alphabet
rng(1); % set the random seed to replicate the simulation
pat = pattern2(Nsymb,'rand',struct('format',modfor));

% pat are labels in decimal form. They will be converted in complex symbols
% or in a binary pattern depending on the operation on them. 

% Let us start converting the symbols into constellation points (ak) and 
% undo (pathat).
% The plot in binary form is not necessary. The output of samp2pat in
% decimal form is necessary since pat is in decimal form in this script.
ak = pat2samp(pat,modfor,struct('plot','bin'));
pathat = samp2pat(ak,modfor,struct('type','dec'));

% % Method 2
% pathat = samp2pat(ak,modfor);
% pathat = bintodec(pathat);

errors = sum(pat ~= pathat);
fprintf('\nErrors: %d\n\n',errors)

