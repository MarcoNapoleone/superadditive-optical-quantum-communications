function pat = diffencoder(pat,modformat)

%DIFFENCODER Differential encoder
%   PAT=DIFFENCODER(PAT,MODFORMAT) given the pattern PAT (either in
%   binary or decimal representation) returns the  pattern differentially
%   encoded. The finite-state machine of the encoder is initialized by 
%   assuming logical zero initial state.
%
%   MODFORMAT is the modulation format (see MODFORMATINFO).
%
%   References:
%       [Seimetz] M. Seimetz, "High-order Modulation for Optical Fiber
%           Transmission", Springer, 2009.
%
%   See also PATTERN, DIFFDECODER.
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
num_bit = log2(p.alphabet);

if size(pat,2)==1 % Not a binary matrix
    pat = dectobin(pat,num_bit);
    isdec = true;
else
    isdec = false;
end

if strcmpi(p.format,'qam') && (p.alphabet == 32 || p.alphabet == 8)
    error('Differential decoding for 8 and 32-QAM not yet implemented.');
end

% Note for the operations below:
%
% y[n] = x[n]*y[n-1]. Product is a sum in log domain -> exp(filter...

switch p.family
    
    case 'ook'
        return
        % nothing to do
        
    case 'psk' % differential coding.
        stars = pat2samp(pat, p.format);
        stars(2:end) = exp(filter(1,[1 1],log(stars(2:end))));
        pat = samp2pat(stars,p.format);
        
    case 'pam'
        stars = pat2samp(pat(:,1), 'bpsk'); % phase ambiguity of pi
        stars(2:end) = exp(filter(1,[1 1],log(stars(2:end))));
        pat(:,1) = samp2pat(stars,'bpsk');
        
    case 'qam'
        pat = dqam(pat,p.alphabet); % prepare pattern ([Seimetz], pp. 43--47)
        stars = pat2samp(pat(:,1:2), 'qpsk'); % phase ambiguity of pi/2
        stars(2:end) = exp(filter(1,[1 1],log(stars(2:end))));
        tmp = samp2pat(stars,'qpsk');
        pat(:,1:2) = tmp;
        
    otherwise
        error('Unknown modulation format.');
end

if isdec
    pat = bintodec(pat);
end

%--------------------------------------------------------------------------
function y = dqam(x,M)

%DQAM prepare pattern to be differentially decoded
%   Y=DQAM(X,M) prepares the binary pattern X of a qam in alphabet M to be
%   differentially decoded. Substantially, the first two bits are equal for
%   each sub-quadrant and differentially decoded, like for a qpsk. The
%   other bits are rotationally invariant. As a result, the pattern in Y is
%   not anymore Gray coded.

n=xor(x(:,1),x(:,end/2+1));
y = ddqam(x,M,n);

%--------------------------------------------------------------------------
function y=ddqam(x,M,n)

b = [x(:,1) x(:,end/2+1)]; % first half: in-phase. Second half: quadrature.
b(n,:) = [b(n,2) b(n,1)];

if M == 4
    y=b; % end of recursion
else
    y=[b,ddqam(x(:,[2:end/2,end/2+2:end]),M/4,n)]; % concatenate recursive calls
end

