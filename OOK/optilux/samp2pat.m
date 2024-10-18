function pat = samp2pat(Iric,modformat,x)

%SAMP2PAT Convert received samples into a pattern (hard decision)
%   PAT = SAMP2PAT(IRIC,MODFORMAT,X) takes the, possibly complex, samples
%   IRIC and takes a hard decision over them. MODFORMAT is the modulation
%   format (e.g., '16qam','qpsk',etc). On output, PAT is the bit-pattern.
%
%   X is a struct for optional parameters:
%
%       X.type: 'dec' or 'bin'. If 'dec' the output PAT is in decimal form.
%           Default: 'bin'.
%       X.threshold: thresholds for the hard-decision. Currently, it is
%           supported by only OOK, BPSK, DPSK, QPSK, and DQPSK. For a real
%           IRIC x.threshold is a single scalar, for complex IRIC it is a
%           vector of size [1,2], the first threshold for the real part,
%           the other for the imaginary.
%           As a reference, the midpoint threshold for OOK is 1, while for 
%           phase-modulated signals is 0.
%
%           Note: for coherent detection, the optimal threshold is
%           automatically set.
%       X.detection: if 'coherent', force optimal threshold for coherent
%           detection. Useful for amplitude modulation schemes such as OOK
%           and PAM that can be detected either coherently or by
%           direct-detection. If not specified, OOK, DPSK, and DQPSK are
%           detected by direct-detection.
%
%   Note 1: the output pattern is Gray when Gray coding is possible.
%
%   Note 2: SAMP2PAT is the inverse operation of PAT2SAMP.
%
%   References:
%
%   [Seimetz] M. Seimetz, "High-order Modulation for Optical Fiber
%           Transmission", Springer, 2009.
%   [Rios] R. R.-Muller, "Joint Coding Rate and Modulation Format
%       Optimization for 8QAM Constellations Using BICM Mutual
%       Information," in Proc. OFC 2016, paper W3K.4, 2016
%
%   See also: RXFRONTEND, RXDSP, PAT2SAMP.
%
%   Author: Paolo Serena, 2021
%   		University of Parma, Italy

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


% Initial checks

% Convert samples to pattern
y = modformatinfo(modformat);
[nr,nc] = size(Iric);
log2M = log2(y.alphabet);
if isinf(log2M), error('It is impossibile to associate a pattern in Inf cardinality');end
pat = false(nr,nc*log2M); % logical array of all 0

switch y.format
    case 'ook' % direct detection
        
        if nargin < 3
            error('Missing threshold');
        else
            if ~isfield(x,'detection') || ~strcmpi(x.detection,'coherent')
                thr = x.threshold;
            else
                thr = 1;
            end
        end
        pat(Iric >  thr) = 1;
        
    case {'bpsk','dpsk'}
        
        if nargin == 4 && isfield(x,'threshold')
            thr = x.threshold;
        else
            thr = 0;
        end
        pat(Iric >  thr) = 1;
        
    case {'qpsk','dqpsk'}
        
        if nargin == 4 && isfield(x,'threshold')
            thr = x.threshold;
        else
            thr = [0 0];
        end
        for k=1:nc
            pat(real(Iric(:,k)) >  thr(1) , 2*k-1) = 1;
            pat(imag(Iric(:,k)) >  thr(2) , 2*k)   = 1;
        end
        
    case 'psk'
        
        phi_s = 3*pi/4;
        samp = mod(round((angle(Iric)+phi_s)*2^(log2M-1)/pi),y.alphabet); % convert in decimal
        %  E.g. 8PSK
        %   ipat (before xor): 0 1 3 2 5 4 6 7
        %   pat              : 0 1 2 3 4 5 6 7
        
        for k=1:nc
            pat1=dectobin(samp(:,k),log2M);
            pat1 = [xor(pat1(:,1),pat1(:,2)),pat1(:,2:end)];
            pat(:,1+(k-1)*log2M:k*log2M) = bintogray(pat1);
        end
        % xor: the first two bits comes from a iq modulator, the
        % others from pm modulator.
        
    case 'qamcirc'
        if y.alphabet ~= 8, error('Alphabet must be 8.'); end
        ALEV = [0; 3; 1; 2; 5; 4; 6];
        LEVEL_data = [cos(2*pi/7*ALEV)+1i*sin(2*pi/7*ALEV); 0]; % 8qam_circ
        Es = 7/8;
        dataid = LEVEL_data/sqrt(nc*Es);
        for k=1:nc
            datahat = nearestsamp(dataid,Iric(:,k));
            pat(:,1+(k-1)*log2M:k*log2M) = dectobin(datahat-1,log2M);
        end
        
    case 'qamrect'
        if y.alphabet ~= 8, error('Alphabet must be 8.'); end
        LEVEL_data = [-3+1i; -3-1i; -1+1i; -1-1i; 3+1i; 3-1i; 1+1i; 1-1i]; % 8-QAM rect
        Es = 6;
        dataid = LEVEL_data/sqrt(nc*Es);
        for k=1:nc
            datahat = nearestsamp(dataid,Iric(:,k));
            pat(:,1+(k-1)*log2M:k*log2M) = dectobin(datahat-1,log2M);
        end
        
    case 'qamstar'
        if y.alphabet ~= 8, error('Alphabet must be 8.'); end
        ps3 = 1+sqrt(3);
        LEVEL_data = [-1-1i; ps3; -1i*ps3; 1-1i; -ps3; 1+1i; -1+1i; 1i*ps3]; % 8qam_star
        Es = (2 + ps3^2)/2;
        pat = getpat(LEVEL_data,nc,Es,log2M,Iric);
        
    case 'pam'
        if ~isfield(x,'detection') || ~strcmpi(x.detection,'coherent')
            warning('optilux:samp2pat','sub-optimal mid-level thresholds')
        end
        Es = (y.alphabet^2-1)/6*2;
        Iric = real(Iric)*sqrt(nc*Es);
        samp_r = pam_hd(y.alphabet,Iric);
        for k=1:nc
            pat1 = dectobin(samp_r(:,k),log2M);
            pat(:,1+(k-1)*log2M:k*log2M) = bintogray(pat1);
        end
        
    case 'qam'
        if y.alphabet == 32 % Cross 32QAM
            LEVEL_data = [-3+5i; -1+5i; -3-5i; -1-5i; -5+3i; -5+1i; -5-3i; -5-1i; ...
                -1+3i; -1+1i; -1-3i; -1-1i; -3+3i; -3+1i; -3-3i; -3-1i; ...
                3+5i;  1+5i;  3-5i;  1-5i;  5+3i;  5+1i;  5-3i;  5-1i; ...
                1+3i;  1+1i;  1-3i;  1-1i;  3+3i;  3+1i;  3-3i;  3-1i]; % 32-QAM
            Es = 20;
            pat = getpat(LEVEL_data,nc,Es,log2M,Iric);
            
        elseif y.alphabet == 128 % Cross 128QAM
            LEVEL_data = [-7+9i; -7+11i; -1+9i; -1+11i; -7-9i; -7-11i; -1-9i; -1-11i; ...
                -5+9i; -5+11i; -3+9i; -3+11i; -5-9i; -5-11i; -3-9i; -3-11i; ...
                -9+7i; -9+5i; -9+1i; -9+3i; -9-7i; -9-5i; -9-1i; -9-3i; ...
                -11+7i; -11+5i; -11+1i; -11+3i; -11-7i; -11-5i; -11-1i; -11-3i; ...
                -1+7i; -1+5i; -1+1i; -1+3i; -1-7i; -1-5i; -1-1i; -1-3i; ...
                -3+7i; -3+5i; -3+1i; -3+3i; -3-7i; -3-5i; -3-1i; -3-3i; ...
                -7+7i; -7+5i; -7+1i; -7+3i; -7-7i; -7-5i; -7-1i; -7-3i; ...
                -5+7i; -5+5i; -5+1i; -5+3i; -5-7i; -5-5i; -5-1i; -5-3i; ...
                7+9i;  7+11i;  1+9i;  1+11i;  7-9i;  7-11i;  1-9i;  1-11i; ...
                5+9i;  5+11i;  3+9i;  3+11i;  5-9i;  5-11i;  3-9i;  3-11i; ...
                9+7i;  9+5i;  9+1i;  9+3i;  9-7i;  9-5i;  9-1i;  9-3i; ...
                11+7i;  11+5i;  11+1i;  11+3i;  11-7i;  11-5i;  11-1i;  11-3i; ...
                1+7i;  1+5i;  1+1i;  1+3i;  1-7i;  1-5i;  1-1i;  1-3i; ...
                3+7i;  3+5i;  3+1i;  3+3i;  3-7i;  3-5i;  3-1i;  3-3i; ...
                7+7i;  7+5i;  7+1i;  7+3i;  7-7i;  7-5i;  7-1i;  7-3i; ...
                5+7i;  5+5i;  5+1i;  5+3i;  5-7i;  5-5i;  5-1i;  5-3i] / 11; % 128-QAM
            Es = mean(abs(LEVEL_data).^2);
            pat = getpat(LEVEL_data,nc,Es,log2M,Iric);
            
        elseif ~mod(log2M,2) % square QAM
            side=sqrt(y.alphabet);
            %             e_r = real(Iric)*sqrt(2); % -1<=e_r<=1
            %             e_i = imag(Iric)*sqrt(2); % -1<=e_i<=1
            Es = (y.alphabet-1)/6*4;
            Iric = Iric*sqrt(nc*Es);
            samp_r = pam_hd(side,real(Iric));
            samp_i = pam_hd(side,imag(Iric));
            
            %   Create the following map for 'samp', e.g.:
            %               ^
            %           s 2 |  0010--0110 | 1110--1010
            %           a 3 |  0011--0111 | 1111--1011
            %           m   |  -----------------------
            %           p 1 |  0001--0101 | 1101--1001
            %           _ 0 |  0000--0100 | 1100--1000
            %           i   |-------------------------->
            %                   0     1      3     2  samp_r
            for k=1:nc
                pat1 = dectobin(samp_i(:,k) +samp_r(:,k)*side,log2M);
                pat(:,1+(k-1)*log2M:k*log2M) = ...
                    [bintogray(pat1(:,1:end/2)),bintogray(pat1(:,end/2+1:end))];
            end
            % E.g., 16QAM
            %
            % ipat: 0  1  3  2  4  5  7  6  12 13 15 14  8  9 11 10
            % pat : 0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15
            
        else
            error('This QAM is not available.');
        end
        
    otherwise
        error('Unknown modulation format.');
end

if nargin == 3 && isfield(x,'type') && strcmpi(x.type,'dec')
    patdec = zeros(nr,nc);
    for k=1:nc
        patdec(:,k) = bintodec(pat(:,1+(k-1)*log2M:k*log2M));
    end
    pat = patdec;
end
%--------------------------------------------------------------------------
function pat = getpat(LEVEL_data,nc,Es,log2M,Iric)

dataid = LEVEL_data/sqrt(nc*Es);
for k=1:nc
    datahat = nearestsamp(dataid,Iric(:,k));
    pat(:,1+(k-1)*log2M:k*log2M) = dectobin(datahat-1,log2M);
end

%--------------------------------------------------------------------------
function samp = pam_hd(Alph,Iric)

Mm1 = Alph-1;
samp = round((Mm1+Iric)/2); % next neighbor
samp(samp<0) = 0;
samp(samp>Alph-1) = Alph-1; % first/last symbols

%--------------------------------------------------------------------------
function n=nearestsamp(aRef,aTest)

%NEARESTSAMP Find the nearest sample
%
%   N=NEAREST(AREF,ATEST) takes a hard-decision on the samples in ATEST and
%   returns in N the index of the nearest sample in the reference set of 
%   values AREF.
%
%   N has the same size of ATEST. AREF has the size of the constellation.
%
%   Note: The decided symbols are thus AREF(N).
%
%   Ref:
%   https://it.mathworks.com/matlabcentral/answers/218033-efficient-method-
%   for-finding-index-of-closest-value-in-very-large-array-for-a-very-
%   large-amount-of-e

x = abs(bsxfun(@minus,aRef(:).',aTest(:)));
[~, n] = min(x,[],2);
