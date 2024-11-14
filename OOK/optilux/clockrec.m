function A = clockrec(A,pat,x)

%CLOCKREC clock recovery (synchronization)
%   [A,RHORATIO] = CLOCKREC(A,PAT,X) recovers the time offset of signal A
%   with respect to the corresponding reference signal of pattern PAT. The 
%   signal A can be a column vector or a two-column matrix. PAT follows the 
%   same rules. The two-column case corresponds to dual-polarization
%   transmission.
%
%   On output, A is a delayed version of the input A.
%
%   RHORATIO is the ratio between the trace of the correlation matrix and
%   the trace of the off-diagonal of the same matrix. Such a number has a
%   meaning only in dual-polarization transmission and is an indication of
%   a possible polarization exchange if smaller than 1.
%
%   X is a struct with fields:
%
%       X.sync = Method for clock recovery. valid options are:
%           X.sync='da': maximum-likelihood, data-aided.
%       X.modformat = modulation format, e.g., '16qam','ook', etc.
%       X.interp: (optional) can be 'fine' indicating that a non-integer
%           time-delay is estimated by interpolation. Such a fractional 
%           delay is recovered by FFT, hence this option, although more 
%           accurate, is slower.
%
%   See also PAT2SAMP, RXDSP.
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

Ncol = size(A,2);   % number of signals
if Ncol ~= size(pat,2)
    error('The pattern must be in decimal form.')
end
rho = zeros(Ncol);
ns = zeros(Ncol);
for m1=1:Ncol
    for m2=1:Ncol
        [ns(m1,m2),rho(m1,m2)] = time_estimate(A(:,m1),pat(:,m2),x);
    end
end
[~,imax]=max(rho,[],2); % best candidate. imax(m1)~= m1 is an indication of signal exchange.
if length(unique(imax)) ~= length(imax)
    warning('optilux:clockrec','Clock recovery by the function clockrec failed.')
else
    A = A(:,imax); % signal exchange (if necessary)
    for m1=1:Ncol
        A(:,m1) = fix_delay(A(:,m1),ns(m1,imax(m1)));
    end
end
%     [~,i1] = unique(imax,'first');
%     indexToDupes = find(not(ismember(1:numel(imax),i1)), 1);
%     if ~isempty(indexToDupes)
%         warning('optilux:clockrec','Clock recovery  failed.');
%     end

%--------------------------------------------------------------------------
function [ns,rho] = time_estimate(A,pat,x)

%TIME_ESTIMATE recover timing mismatch
%   NS = TIME_ESTIMATE(A,X) estimates the timing offset NS [samples] of
%   the signal A by non-data aided algorithms (to be implemented).
%
%   [NS,RHO] = TIME_ESTIMATE(A,PAT,X) works in data-aided mode using the
%   data pattern PAT. RHO is the peak of the correlation with the pattern.
%

if nargin == 2
    if ~isstruct(pat)
        error('Non-data aided methods not available.');
    end
    x = pat;
end
if ~isfield(x,'type')
    error('Unknown timing recovery method.');
end

%%% timing recovery
if strcmpi(x.type,'da')
    [ns,rho] = mlda_time(A,pat,x);
elseif strcmpi(x.type,'') || isempty(x.type)
    ns = NaN; rho = NaN;
else
    error('Unknown timing recovery method.');
end

%--------------------------------------------------------------------------
function [ns,rho] = mlda_time(A,pat,x)

%%% Init
Nfft = length(A);
Nsymb = length(pat);
Nt = Nfft/Nsymb;

if rem(Nt,1) ~= 0
    error('fractional number of samples per symbol has not yet done.')
end

%%% go

akid = pat2samp(pat,x.modformat);

aku = upsample(akid,Nt);
Aku = fft(aku);
corr = abs(ifft( conj(Aku) .* fft(A) )); % correlation
[rho,imax] = max(corr);
ns = imax-1;

%%% fine tuning
if isfield(x,'interp')
    if strcmpi(x.interp,'fine') %% fine tuning
        i1 = nmod(imax-1,Nfft); i2 = nmod(imax,Nfft); i3 = nmod(imax+1,Nfft);
        % max of parabola passing for x-coordinates [1 2 3] with values corr.
        ifine = (5*corr(i1)+3*corr(i3)-8*corr(i2))/(corr(i1)+corr(i3)-...
            2*corr(i2))/2;
        ns = nmod(ifine+ns-2,Nfft); % -2: because central point is at coordinate 2
    elseif strcmpi(x.interp,'nearest') %% coarse (but faster) tuning
        % nothing to do
    else
        error('Unknown time interpolation method.')
    end
end

%--------------------------------------------------------------------------
function A = fix_delay(A,ns)

%FIX_DELAY apply delay recovery

global GSTATE

for k=1:length(ns)
    if floor(ns(k)) == ns(k)
        A(:,k) = fastshift(A(:,k),-ns(k));
    else
        omega = 2*pi*GSTATE.FN(:)/GSTATE.FSAMPLING;
        A(:,k) = ifft( fft(A(:,k)) .* fastexp(omega*ns(k)) );
    end
end
