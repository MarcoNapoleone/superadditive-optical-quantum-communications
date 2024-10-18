function y = modformatinfo(modformat)

%MODFORMATINFO Get basic info about the digital modulation format
%   Y = MODFORMATINFO(MODFORMAT) returns basic information about the
%   digital modulation format MODFORMAT. MODFORMAT is a string, possibly
%   containing the alphabet size, e.g., '64qam','4pam','ook','qpsk',etc. 
%   Currently supported formats are:
%
%   'nqam': quadrature amplitude modulation of alphabet n, with n an
%       integer power of 2.
%   'bpsk': binary phase shift keying
%   'qpsk': quadrature phase shift keying
%   'npsk': phase-shift keying of alphabet n, with n an integer power of 2.
%   'dpsk': differential-phase shift keying
%   'dqpsk': differential quadrature phase shift keying
%   'ook': on-off keying
%   'npam': pulse amplitude modulation of alphabet n, with n an integer
%       power of 2.
%   '8qam_star': 8qam in star configuration (_ not necessary) [Rios]
%   '8qam_rect': 8qam in rectangular configuration (_ not necessary) [Rios]
%   '8qam_circ': 8qam in circular configuration (_ not necessary) [Rios]
%
%   See the OptiluX examples directory on how to see the constellations.
%
%   On output:
%
%       Y.alphabet: alphabet size, e.g., 64 for '64qam'.
%       Y.family:   general family of the modulation format, e.g., 'qam'
%           for '64qam', 'qam' for '8qam_rect', 'psk' for 'dpsk', etc
%       Y.format:   modulation format type. E.g., 'qam' for '64qam',
%           '8qamrect' for '8qam_rect', 'dpsk' for 'dpsk'.
%           Y.format returns a unique name for different input
%           possibilities, e.g., '8qamcirc' for both, identical, inputs
%           '8qam_circ' and '8qamcirc'.
%       Y.symb.mean: average value of the symbols.
%       Y.symb.var: variance of the (possibly complex) symbols.
%
%   References:
%
%   [Rios] R. R.-Muller et al., "Joint Coding Rate and Modulation Format
%       Optimization for 8QAM Constellations Using BICM Mutual
%       Information," in Proc. OFC 2016, paper W3K.4, 2016
%
%   Author: Paolo Serena, 2021
%   University of Parma, Italy

%    This file is part of Optilux, the optical simulator toolbox.
%    Copyright (C) 2021  Paolo Serena, <serena@tlc.unipr.it>
%
%    Optilux is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 3 of the License, or
%    (at your option) any later version.
%
%    Optilux is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

if ~ischar(modformat)
    error('The modulation format must be a string')
end

y.alphabet = str2num(modformat(isstrprop(modformat,'digit')));
y.format = modformat(isstrprop(modformat,'alpha'));

if strcmpi(modformat,'bpsk') || strcmpi(modformat,'dpsk')
    y.alphabet = 2;
    y.family = 'psk';
    y.symb.mean = 0;
    y.symb.var = 1;
elseif strcmpi(modformat,'ook')
    y.alphabet = 2;
    y.family = 'ook';
    y.symb.mean = 1;
    y.symb.var = 1;
elseif strcmpi(modformat,'qpsk') || strcmpi(modformat,'dqpsk') || ...
        (strcmpi(y.format,'qam') && y.alphabet == 4)
    y.alphabet = 4;
    y.family = 'psk';
    y.symb.mean = 0;
    y.symb.var = 1;
elseif strcmpi(y.format,'psk')
    y.family = 'psk';
    y.symb.mean = 0;
    y.symb.var = 1;
elseif (strcmpi(y.format,'qamcirc') || strcmpi(y.format,'qamstar') || ...
        strcmpi(y.format,'qamrect')) && y.alphabet == 8
    y.family = 'qam';
    y.symb.mean = 0;
    y.symb.var = 1;
elseif strcmpi(y.format,'qam')
    y.family = 'qam';
    y.symb.mean = 0;
    y.symb.var = 1;
elseif strcmpi(y.format,'pam')
    y.family = 'pam';
    y.symb.mean = 0;
    y.symb.var = 1;
elseif strcmpi(y.format,'randn')
    y.family = 'randn';
    y.symb.mean = 0;
    y.symb.var = 1;
    y.alphabet = Inf;
else
    error('Unknown modulation format');
end

if isempty(y.alphabet)
    error('Missing modulation format alphabet size')
end
