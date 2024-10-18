function varargout=pattern2(Nsymb,ptype,varargin)

%PATTERN sequence pattern
%   PAT=PATTERN(NSYMB,PTYPE,NSEED) returns in PAT a sequence of NSYMB
%   binary symbols representing the pattern of a digital modulation. The
%   alphabet of the corresponding digital constellation can be set by the 
%   optional parameters OPTIONS, see below.
%
%   PTYPE is the type of the pattern and can be one of the following:
%
%       'debruijn': creates a De Bruijn sequence (DBS) (each subsequence of
%               length log2(NSYMB) appears exactly once in a DBS) [1-2]. In
%				the binary case, a DBS is a pseudo-random binary sequence
%			    (PRBS) with an additional zero added to to the longest
%				sequence of 0.
%
%               With PTYPE='debruijn' NSEED is the DBS seed.
%
%               Notes:
%               - 0 <= NSEED < NSYMB/4 yields a unique DBS, i.e., it is not
%               possible to obtain the same DBS with a different NSEED,
%               neither with a circular shift.
%               - For NSEED >= NSYMB/4 the DBS is not unique, hence PATTERN
%               generates a circular random delayed version of a DBS
%               with NSEED < NSYMB/4.
%   	        - NSEED must be < NSYMB/4*(NSYMB-1).
%               - It is not possible to have the same sequence for
%               different NSEED.
%               - Please note hat while a DBS may be best for a single
%               signal, different DBS may not have good correlation
%               properties.
%
%       'rand': create a uniformly distributed random sequence.
%
%               PAT=PATTERN(NSYMB,'rand') does not set the random seed,
%               which can thus be set before calling PATTERN with RNG or
%               RAND with the field 'state'.
%
%               PAT=PATTERN(NSYMB,'rand',OPTIONS) does not set the random
%               seed and uses optional arguments in OPTIONS (see later).
%
%               PAT=PATTERN(NSYMB,'rand',NSEED) sets the state of the
%               random generator to NSEED. Please note that in such a case
%               the original random seed is restored after the call to
%               PATTERN.
%
%               PAT=PATTERN(NSYMB,'rand',NSEED,OPTIONS) is the most
%               general call that sets the random seed and accepts optional
%               parameters.
%
%       'randn': create a normally distributed complex sequence. In this
%               case the sequence is complex with standard normally
%               distributed real and imaginary components (zero mean,
%               complex variance one). The random seed is set with the same
%               rules of PTYPE='rand'.
%
%       'randnreal': same as 'randn', but with real numbers distributed
%               with the standard normal distribution. The random seed is
%               set with the same rules of PTYPE='rand'.
%
%       <sequence of numbers>: a string or a vector of numbers.
%               E.g. '11011100' creates a periodic repetition of the
%               sequence up to length NSYMB, and truncates when necessary.
%               For instance, with NSYMB=16 such a pattern returns 
%               1101110011011100.
%               The sequence can be a string or a vector of double, e.g.
%               '20104101' or [2 0 1 0 4 1 0 1]
%
%               With a sequence of numbers, NSEED makes no sense and must
%               not be provided to the function.
%
%               PAT = PATTERN(NSYMB,<sequence of numbers>,OPTIONS) uses
%               optional arguments in OPTIONS (see below).
%
%       <file>: reads the pattern from 'file' using LOAD. 'file' contains
%               the pattern for the channel.
%               E.g. NSYMB=8, with patterns '30101210', 'file' can have the
%               form, e.g.:
%
%                   3 0 1 0    or        3 0 1 0 1 2 1 0
%                   1 2 1 0
%
%               With this option, NSEED makes no sense and must not be
%               provided to the function.
%
%
%   OPTIONS is an optional parameter containing:
%
%       OPTIONS.alphabet: is the alphabet of the pattern.
%           E.g. PAT=PATTERN(NSYMB'debruijn',0,OPTIONS) with
%               OPTIONS.alphabet=4 returns:
%
%           PAT=[1  1  2  3  0  3  1  3  3  2  2  1  0  2  0  0]
%
%           i.e. a DBS sequence with symbols (0,1,2,3) containing all
%           couples of symbols exactly once.
%
%           E.g. PAT=PATTERN(NSYMB,'rand',OPTIONS) with
%               OPTIONS.alphabet=8 may return:
%
%           PAT=[6  3  7  6  6  2  3  7  3  2  5  4  4  2  0  2]
%
%       OPTIONS.format: alternatively to OPTIONS.alphabet, the modulation
%           format can be directly indicated, e.g., OPTIONS.format='qpsk'.
%           The list of currently supported formats is in MODFORMATINFO,
%
%
%   [PAT,BMAT] = PATTERN(NSYMB,PTYPE,NSEED,OPTIONS) returns in BMAT the
%   binary representation of PAT. BMAT is a matrix of size
%	[NSYMB,ceil(log2(OPTIONS.alphabet))].
%
%   Example: For a quadrature phase-shift keying (QPSK) modulation, the two
%       columns of BMAT represent the bits driving the in-phase and
%       quadrature components.
%
%
%   See also DIGITALMOD, DIFFENCODER, DIFFDECODER.
%
%
%   [1] F. S. Annexstein, "Generating De Bruijn Sequences: An Efficient
%   Implementation," IEEE Transaction on Computers, vol. 48, no. 2,
%   pp. 198-200, Feb. 1997.
%
%   [2] D. v. d. Borne et al, "Bit pattern dependence in optical DQPSK
%   modulation," Electronic Letters, vol. 43, no. 22, Oct. 2007
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



if strcmp(ptype,'debruijn') %%%%% DE BRUIJN SEQUENCE
    
    if nargin < 3
        error('Missing seed');
    else
        nseed = varargin{1};
        if length(nseed) > 1, error('The seed must be a scalar.');end
        if nargin == 4
            options = varargin{2};
            nfnames = fieldnames(options);
            if any(strcmp(nfnames,'alphabet'))
                q = log2(varargin{2}.alphabet);
            elseif any(strcmp(nfnames,'format'))
                mf = modformatinfo(options.format);
                q = log2(mf.alphabet);
            else
                q = 1;
            end
            if rem(q,1)
                error(['De Bruijn sequence can be implemented only for ',...
                    ' power of 2 alphabets']);
            end
            
        else
            q = 1; % default: binary
        end
        maxseed = Nsymb*(Nsymb)/4-1; % largest seed
        if nseed > maxseed, error('nseed must be <= %d',maxseed); end
        lg2 = log2(Nsymb);
        if rem(lg2,q)
            error(['The De Bruijn sequence does not exist! log2(number ',...
                'of symbols)=%d is not a multiple of log2(alphabet length)=%d.'],...
                lg2,q);
        end
        [pat,tmat] = debruijn_seq(log2(Nsymb)/q,nseed,q);
        varargout{2} = tmat';
    end
    
    
elseif strcmp(ptype,'rand') %%%%% RANDOM UNIFORMLY-DISTRIBUTED PATTERN
    setstate = false;
    if nargin == 2
        qq = 2; % alphabet
        
    elseif nargin == 3
        if isnumeric(varargin{1})
            [isrng,oldrandst] = setrandseed(varargin{1});
            qq = 2;
            setstate = true;
        elseif isstruct(varargin{1})
            nfnames = fieldnames(varargin{1});
            if any(strcmp(nfnames,'alphabet'))
                qq = varargin{1}.alphabet;
            elseif any(strcmp(nfnames,'format'))
                mf = modformatinfo(varargin{1}.format);
                qq = mf.alphabet;
                if strcmpi(mf.format,'randn')
                    error('Cannot use uniform-distributed symbols with a randn format.');
                end
            else
                qq = 2;
            end
        else
            error('Second input argument must be an integer or a struct.');
        end
        
    elseif nargin == 4
        nfnames = fieldnames(varargin{2});
        if any(strcmp(nfnames,'alphabet'))
            qq = varargin{2}.alphabet;
        elseif any(strcmp(nfnames,'format'))
            mf = modformatinfo(varargin{2}.format);
            qq = mf.alphabet;
        else
            qq = 2;
        end
        if ~isnumeric(varargin{1})
            error('Random seed must be a double scalar.');
        end
        [isrng,oldrandst] = setrandseed(varargin{1});
        setstate = true;
    end
    pat = floor(qq.*rand(1,Nsymb));
    if setstate
        if isrng
            rng(oldrandst);
        else
            rand('state',oldrandst); % restore old rand state
        end
    end
    varargout{2} = mydec2bin(pat,log2(qq));
    
elseif strcmp(ptype,'randn') % NORMAL-DISTRIBUTED PATTERN
    if nargin > 2 && isnumeric(varargin{1})
        [isrng,oldrandnst] = setrandnseed(varargin{1});
        setstate = true;
    else
        setstate = false;
    end
    pat = 1/sqrt(2)*(randn(1,Nsymb) +1i*randn(1,Nsymb));
    if setstate
        if isrng
            rng(oldrandnst);
        else
            randn('state',oldrandnst); % restore old rand state
        end
    end
    
    varargout{2} = NaN; % bits do not exist
    
elseif strcmp(ptype,'randnreal') % NORMAL-DISTRIBUTED PATTERN
    pat = randn(1,Nsymb);
    varargout{2} = NaN; % bits do not exist
    
elseif isnumeric(ptype)     %%%%% VECTOR OF SYMBOLS
    if length(ptype) > Nsymb
        warning('optilux:pattern','Pattern too long: truncated.');
    end
    
    if min(ptype) < 0, error('Pattern cannot be < 0.');end
    if nargin == 3
        if ~isstruct(varargin{1})
            error('Second input argument must be a struct.');
        end
        nfnames = fieldnames(varargin{1});
        if any(strcmp(nfnames,'alphabet'))
            Alph = varargin{1}.alphabet;
            q = log2(Alph);
        elseif any(strcmp(nfnames,'format'))
            mf = modformatinfo(varargin{1}.format);
            Alph = mf.alphabet;
            q = log2(Alph);
        else
            q = 1;
        end
        if max(ptype) > Alph-1
            error('Pattern cannot be > Alphabet size.');
        end
        if rem(q,1)
            error(['Sequence can be implemented only for ',...
                ' power of 2 alphabets']);
        end
        
    elseif nargin >= 4
        error('Too many input arguments.');
    else
        if max(ptype) > 1
            error('Pattern cannot be > Alphabet size.');
        end
        q = 1; % default: binary
    end
    pat = myseq(ptype,Nsymb); % periodically repeat the pattern
    varargout{2} = mydec2bin(pat,q);
    
elseif (min(double(ptype-48)) >= 0 && max(double(ptype-48)) <= 9)   %%%%% STRING
    ptype2 = double(ptype-48);
    if nargin >= 4
        error('Too many input arguments.');
    elseif nargin == 3
        [pat,patmat] = pattern2(Nsymb,ptype2,varargin{end});
    else
        [pat,patmat] = pattern2(Nsymb,ptype2);
    end
    varargout{2} = patmat;
    
elseif exist(ptype,'file') == 2 %%%%% LOAD FROM FILE
    
    pat = loadpattern(ptype,Nsymb);
    varargout{2} = mydec2bin(pat);
    
else
    
    error('String ptype unknown');
end
varargout{1} = pat(:);

%--------------------------------------------------------------------------
function y=myseq(locpat,Nsymb)

% creates a periodic repetition of the sequence locpat up to Nsymb, and
% truncates if necessary.
%
% E.g.  ptype = [0 1 0] and Nsymb = 8 -> y=[0 1 0 0 1 0 0 1]
%       ptype = [0 1 0 0 1 0] and Nsymb = 4 -> y=[0 1 0 0]
%
%   Author: Paolo Serena, 2009
%   University of Parma, Italy

Lp = length(locpat);
nrem = rem(Nsymb,Lp);
if Lp > Nsymb
    y = locpat(1:Nsymb);
elseif nrem ~= 0
    temp1 = repmat(locpat,1,floor(Nsymb/Lp));   % replicate
    temp2 = locpat(1:nrem);
    y = [temp1,temp2];
else
    y = repmat(locpat,1,Nsymb/Lp);
end


%--------------------------------------------------------------------------
function y=loadpattern(ptype,Nsymb)

% load the pattern from file ptype using the rules of function load
%
%   Author: Paolo Serena, 2007
%   University of Parma, Italy


locpat = load(ptype);
[nr,nc] = size(locpat);
fprintf('%d %d %d\n',nr,nc,Nsymb)
if (nr*nc) ~= (Nsymb)
    error(['File ptype ',...
        'contains a wrong number of symbols']);
end
y = reshape(locpat',1,Nsymb);

%--------------------------------------------------------------------------
function [y,tmat]=debruijn_seq(n,seed,q)

%DEBRUIJN_SEQ De Bruijn sequence generator
%   Y=DEBRUIJN_SEQ(N,SEED) generates a De Bruijn sequence of 2^N bits.
%   A De Bruijn sequence contains all patterns of N bits exactly once.
%
%   SEED is the random seed (integer). For SEED > 2^(N-2) the resulting
%   sequence is a pseudo-random shifted copy of a sequence with
%   SEED <= 2^(N-2).
%
%   Reference:
%
%   [1] F. S. Annexstein, "Generating De Bruijn Sequences: An Efficient
%   Implementation," IEEE Transaction on Computers, vol. 48, no. 2, Feb.
%   1997.

ns = q*n; % log2 length of the DeBruijn sequence
N = 2^(ns-2);
nseed = mod(seed,N);
x = double(dec2bin(nseed,ns-2)-48); % convert to pattern
y = double(generate_debruijn(ns,x));

if q > 1 % multilevel DeBruijn
    mvm = floor(ns/q*((1:q)-1));
    tmat = zeros(q,2^ns); % columns: binary DeBruijn
    for kk=1:q
        tmat(kk,:) = fastshift(y,mvm(kk));
    end
    y = 2.^(q-1:-1:0) * tmat ; % bin2dec conversion
else
    tmat = y;
end
if seed >= N % no more seeds -> apply a random delay shift
    nseed2 = ceil((seed-N+1)/N);
    nshift = mod(97*nseed2,2^ns-1)+1; % congruent random generator
    y = fastshift(y,-nshift); % shift on the right
    if q > 1, tmat = fastshift(tmat.',-nshift).';end
end

%--------------------------------------------------------------------------
function y=generate_debruijn(n,x)

if n == 2
    y=[1,1,0,0] == 1;
elseif n == 3
    if x(1) == 0
        y = [1,0,1,1,1,0,0,0] == 1;
    else
        y = [1,1,1,0,1,0,0,0] == 1;
    end
else
    y=next_debruijn(generate_debruijn(n-1,x),2^(n-2)+(-1)^x(n-3),x(n-2));
end

%--------------------------------------------------------------------------
function y=next_debruijn(w,i,k)

C = xor_iscan(w);
Cbar = not(C);
part1 = C(1:i-k);
part2 = Cbar(i+k:end);
part3 = Cbar(1:i-1+k);
part4 = C(i-k+1:end);
y = [part1,part2,part3,part4];

%--------------------------------------------------------------------------
function y=xor_iscan(x)

% If x=[x1,x2,...,xn] and y=[y1,y2,...,yn] it is:
%       yk = x1 xor x2 xor ... xk

z = cumsum(x);
y = rem(z,2) == 1;

%--------------------------------------------------------------------------
function y=mydec2bin(d,q)

[f,e]=log2(max(d)); %#ok How many digits do we need to represent the numbers?
y=logical(rem(floor(d(:)*pow2(1-max(1,e):0)),2)); % force to logical
cols = size(y,2);
if q > cols
    y(:,q-cols+(1:cols)) = y(:,1:cols);
    y(:,1:(q-cols)) = zeros(size(y,1),q-cols);
end

%--------------------------------------------------------------------------
function [isrng,oldrandst] = setrandseed(seed)

if exist('rng','builtin') == 5
    oldrandst = rng;
    if strcmp(oldrandst.Type,'legacy')
        rand('state',seed);
        isrng = false;
    else
        rng(seed); % set state
        isrng = true;
    end
else % Matlab prior to 7.7
    isrng = false;
    oldrandst = rand('state'); % recovered later
    rand('state',seed); % set state
end

%--------------------------------------------------------------------------
function [isrng,oldrandnst] = setrandnseed(seed)

if exist('rng','builtin') == 5
    oldrandst = rng;
    if strcmp(oldrandnst.Type,'legacy')
        randn('state',seed);
        isrng = false;
    else
        rng(seed); % set state
        isrng = true;
    end
else % Matlab prior to 7.7
    isrng = false;
    oldrandnst = randn('state'); % recovered later
    randn('state',seed); % set state
end
