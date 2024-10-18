function varargout = rxdsp(Iric,symbrate,o1,o2)

%RXDSP digital signal processing
%   AK = RXDSP(IRIC,SYMBRATE,PAT,X) applies digital signal processing to
%   the signal IRIC. SYMBRATE is the symbol rate [Gbaud], PAT is the symbol
%   pattern, either in decimal or binary format. X is a struct with the
%   following fields:
%
%       X.eftype = electrical filter (LPF) type (see MYFILTER). Such a
%           filter is an antialiasing filter before ideal sampling.
%       X.ebw = LPF bandwidth normalized to SYMBRATE.
%       X.epar = electrical filter extra parameters (optional, see
%           MYFILTER).
%       X.sync.type = time recovery method. 'da': data-aided by maximum
%           likelihood estimation. Se also CLOCKREC.
%       X.sync.interp = 'nearest' time recovery with error of +/-1 sample 
%           (fast). 'fine': time recovery by FFT interpolation (default).
%       X.modformat = modulation format, e.g., '16qam','qpsk',etc.
%
%   [AK OUT] = RXDSP(...) returns on output the struct OUT containing some 
%   information about the operation performed by the DSP:
%
%       OUT.sync.rhoratio = ratio of the trace of the autocorrelation
%           matrix along time and the trace of the off-diagonal. It yields
%           information about a polarization exchange.
%       OUT.agc.eta = complex gain applied by the automatic-gain control
%           unit.
%
%   See also RXFRONTEND.
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

global GSTATE

%%% INIT

if nargin == 3
    if isstruct(o1)
        x = o1;
    else
        error('Missing dsp parameters (struct variable).');
    end
elseif nargin == 4
    pat = o1;
    if isstruct(o2)
        x = o2;
    else
        error('Missing dsp parameters (struct variable).');
    end
else
    error('Not enough input parameters')
end

if isfield(GSTATE,'FSAMPLING')
    Nt = GSTATE.FSAMPLING/symbrate; % number of points per symbol
    if floor(Nt) ~= Nt, error('Number of points per symbol is not an integer.'); end
else
    Nt = 1;
end
% Get pattern in decimal form
Nsig = size(Iric,2);
[Nsymb,Ncolp] = size(pat);
if Nsig ~= Ncolp
    p = modformatinfo(x.modformat);
    num_bit = log2(p.alphabet);
    patdec = zeros(Nsymb,Nsig);
    for k=1:Nsig
        patdec(:,k) = bintodec(pat(:,1+(k-1)*num_bit:k*num_bit));
    end
else
    patdec = pat;
end

%%%

% 1: apply the post-detection electrical filter (ADC filter)
Iric = lpfil(Iric,symbrate,x);

% 2: recover time delay
xsync = x.sync; xsync.modformat = x.modformat;
if ~isfield(x.sync,'interp'), xsync.interp='fine';end % default
Iric = clockrec(Iric,patdec,xsync);

% 3: Sample (or downsample and MIMO)
ak = Iric(1:Nt:end,:);

% 4: automatic gain control
for k=1:size(ak,2)
    [ak(:,k) eta(k)] = agc(patdec(:,k),ak(:,k),x.modformat);
end
out.agc.eta = eta;


% output
varargout{1} = ak;
if nargout == 2
    varargout{2} = out;
end

%--------------------------------------------------------------------------
function Iric = lpfil(Iric,symbrate,x)

global GSTATE

if isfield(GSTATE,'FN')
    Fnorm = GSTATE.FN/symbrate;
else
    Nsymb = length(Iric);
    Fnorm = -1/2:1/Nsymb:1/2-1/Nsymb;
end

if ~isfield(x,'epar'), x.epar = NaN; end
Hf = myfilter(x.eftype,Fnorm,x.ebw,x.epar); % lowpass filter
for k=1:size(Iric,2)
    Iric(:,k) = ifft(fft(Iric(:,k)) .* Hf);
end
