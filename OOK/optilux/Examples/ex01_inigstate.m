%   ex01: understanding GSTATE.
%
%   NOTE: BEFORE STARTING an OptiluX simulation you need to add the Optilux
%   directory to the path by the command addpath.

clear
clc

addpath(['..',filesep]) % PUT HERE THE OPTILUX DIRECTORY! If you just  
%           installed Optilux, such a directory is expected to be one step 
%           above, hence the reason for ['..',filesep] which is translated 
%           in ../ in Linux and ..\ in Windows. 
%           You can use the command addpath JUST ONCE PER SESSION. It is 
%           not necessary to type it within a file like here, just type in 
%           the command window before running any OptiluX file.

%% Global parameters
global GSTATE

Nsymb = 1024;           % number of symbols
Nt = 32;                % number of discrete points per symbol

%% Tx parameters
symbrate = 10;          % symbol rate [Gbaud]

%% Init global variables
Nsamp = Nsymb*Nt;       % overall number of samples
fs = symbrate*Nt;       % sampling rate [GHz]
inigstate(Nsamp,fs);    % initialize global variables: Nsamp and fs. 

% inigstate initializes basic global variables, contained in the global
% struct variable GSTATE. OptiluX operated over signals, and such signals
% have some basic common grounds, such as the number of samples.
%
% This is the most complete call to inigstate. You may call
% inigstate(Nsamp) without fs if you are working with a discrete-time
% channel model with one sample per symbol. See next examples for more
% information.

GSTATE

% GSTATE is a global variable containing the minimum amount of information
% for the signals processed by the OptiluX functions. The fields are:
%
% GSTATE.FN: the discrete frequencies [GHz] used by the Optilux
%   sub-functions.
% GSTATE.NSAMP: the number of samples
% GSTATE.FSAMPLING: the sampling rate [GHz]