function varargout=evaleye(pat,sig,symbrate,modformat,x)

%EVALEYE Evaluate the eye-opening
%   EYEOP=EVALEYE(PAT,SIG,SYMBRATE,MODFORMAT,X) evaluates the eye
%   opening EYEOP [dB] of signal SIG having symbol rate SYMBRATE [Gbaud].
%   MODFORMAT is a string indicating the modulation format (see
%   MODFORMATINFO).
%
%   SIG and PAT may be matrices with correspondings signal/patterns on
%   columns, respectively.
%
%   For non-binary modulation formats, EYEOP is the worst case among all
%   symbol eye openings.
%
%   X is an optional struct with values:
%
%       X.plot = true: plot the eye in the current figure.
%
%   See also PATTERN, RXFRONTEND
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

if isfield(GSTATE,'FSAMPLING')
    Nt = GSTATE.FSAMPLING/symbrate; % number of points per symbol
    if floor(Nt) ~= Nt, error('Number of points per symbol is not an integer.'); end
    Nsymb = GSTATE.NSAMP/Nt;
else
    error('no eye at one samples per symbol');
end

npol = size(sig,2);

nshift = round(Nt/2);%round(halfbit-delay*GSTATE.NT); % the first bit is centered at index 1
Iricmat = reshape(fastshift(sig,nshift),Nt,Nsymb*npol).';  % Note the transpose!

y = modformatinfo(modformat);

botv = zeros(y.alphabet,Nt); topv = zeros(y.alphabet,Nt);

for k=1:y.alphabet
    eyesig = Iricmat(pat==(k-1),:);
    topv(k,:) =  min(eyesig,[],1); % top of eye
    botv(k,:) =  max(eyesig,[],1); % bottom of eye
end
[~, itop] = sort(max(topv,[],2)); % sort because of pattern
[~, ibot] = sort(min(botv,[],2));

% eye opening: best of the worst among symbols
eyeop = max(topv(itop(2:end),:) - botv(ibot(1:end-1),:),[],2); % among samples
eyeop = min(eyeop,[],1); % among alphabet
eyeop(eyeop<0) = NaN;

eyeop = 10*log10(eyeop); % [dBm]

if nargout == 1
    varargout{1} = eyeop;
end

if nargin == 5 && isfield(x,'plot') && x.plot
    tim = -1/2:1/Nt:1/2-1/Nt;
    if npol == 1
        plot(tim,Iricmat','b')
    else
        plot(tim,Iricmat(1:end/2,:)','b',tim,Iricmat(end/2+1:end,:)','r')
        title('blue: x. red: y')
    end
    xlabel('normalized time [symbols]')
    ylabel('eye')
    grid on
end
