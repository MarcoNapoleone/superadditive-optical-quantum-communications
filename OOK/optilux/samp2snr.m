function SNR = samp2snr(pat,sig,modformat,x)

%SAMP2SNR Estimate signal-to-noise ratio from samples
%   SNR = SAMP2SNR(PAT, SIG, MODFORMAT, X) estimates the signal-to-noise
%   ratio (SNR) in [dB] of the samples contained in SIG. MODFOR is a string
%   indicating the modulation format (e.g., '16qam','qpsk', etc). PAT is
%   the symbol pattern of SIG.
%
%   SIG can have up to two columns, i.e., one or two polarizations,
%   respectively. In dual-polarization the first and the second half of PAT
%   columns contain the x and y pattern, respectively.
%
%   The noise is extracted by a difference between SIG and the expected
%   signal with pattern PAT, thus assuming that there is no phase mismatch.
%
%   X is a struct of fields:
%       X.power: if 'mean' the average power of SIG is used. If a number,
%           X.power [dBm] is used as signal power in the SNR. If '', the
%           average power is assumed to be 1. Default: ''.
%
%   SNR is usually called Es/N0 in textbooks.
%
%   Note: if SIG is a signal impaired by a realization of additive Gaussian 
%       noise, as a rule of thumb the accuracy of SNR is about 0.1 dB for 
%       4096 samples, and scales as sqrt(number of samples).
%
%   [1] M . C. Jeruchim et al., Simulation of Communication Systems, 
%       Kluwer, 2nd ed., p. 670.
%   [2] R. E. Blahut, Principles and Practice of Information Theory,
%       Addison Wesley, 1988.
%
%   See Also: BER2Q, Q2BER.
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


p = modformatinfo(modformat);

Ncsig = size(sig,2);
Ncpat = size(pat,2);
if ~(Ncpat == Ncsig || Ncpat == Ncsig*log2(p.alphabet))
    error('The number of columns of the pattern are inconsistent with the number of columns of the data.');
end
if size(sig,2) == 1
    dataid = pat2samp(pat, modformat);
else
    dataid(:,1) = pat2samp(pat(:,1:end/2), modformat);
    dataid(:,2) = pat2samp(pat(:,end/2+1:end), modformat);
end

% real/imaginary point of view
if nargin == 4 && isfield(x,'power')
    if strcmpi(x.power,'mean')
        sigpow = 10*log10(sum(mean(abs(dataid).^2))); % sum: x+y
    else
        sigpow = x.power; % [dBm]
    end
else
    sigpow = 10*log10(size(dataid,2)); % theoretical power of an infinite-length 
    % signal, hence: 3 dBm in dual polarization, 0 dBm in single
    % polarization.
end
noise = dataid - sig; % additive noise
noisepow = 10*log10(sum(var(noise))); % noise power X+Y [dBm]
SNR = sigpow - noisepow; % [dB]
