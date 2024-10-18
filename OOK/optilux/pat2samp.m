function dataid = pat2samp(pat, modformat, x)

%PAT2SAMP convert pattern in physical signals
%   DATAID=PAT2SAMP(PAT,FORMAT) converts the pattern contained in PAT into
%   physical signals, i.e., the constellation symbols of the modulation
%   format. PAT2SAMP works for one signal.
%
%   It turns out that for a linearly modulated digital signal:
%
%   	sum_k ( ak * p(t-k*T) )
%
%   PAT2SAMP returns the levels ak in DATAID.
%
%   Note: while the sample mean square value of DATAID may differ from one,
%       the mean square value is always 1. 
%
%   PAT is a column vector containing the symbols of the signal to be
%   converted. Alternatively, PAT can be a matrix of bits whose number of
%   columns is log2(alphabet size) of the modulation format.
%
%   FORMAT is the modulation format, see MODFORMATINFO.
%
%   DATAID=PAT2SAMP(PAT,FORMAT,X) with X.plot='bin' plots the
%   constellation of the modulation format in the current figure with all
%   the corresponding bits as text markers. If X.plot='dec', the decimal
%   representation is used for the text marker.
%
%   Note: PAT2SAMP is the inverse operation of SAMP2PAT.
%
%   References:
%
%   [Seimetz] M. Seimetz, "High-order Modulation for Optical Fiber
%           Transmission", Springer, 2009.
%
%   See Also: PATTERN, SAMP2PAT, DIGITALMOD.
%
%   Author: P. Serena, 2021
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

%%%%%%% Init and checks

if strcmp(modformat,'nf-dpsk')||strcmp(modformat,'nf-dqpsk')||strcmp(modformat,'psbt')
    error('%s not implemented. Take a look at the old versions of Optilux.',...
        modformat);
end

p = modformatinfo(modformat);
num_bit = log2(p.alphabet);

if (max(max(pat)))>(p.alphabet-1)
    error('Pattern of %s must be in alphabet %d.', modformat,p.alphabet);
end

[nrow,ncol] = size(pat); % ncol == 1: symbol. ncol == 2: bits
if strcmpi(modformat,'randn')
    dataid = pat;
    return
end
if ncol == 1
    pat = dectobin(pat,num_bit);
else
    if nrow == 1
        pat = dectobin(pat(:),num_bit);
    end
    if mod(size(pat,2),num_bit) ~= 0
        error('wrong modulation format');
    end
end
% from now on, pat is a binary matrix
dataid = finalizepat2samp(pat, p);

%%%%%%%%%%%% Plot Constellation
if nargin == 3
    if isfield(x,'plot')
        grid on, hold on
        axis square;
        xlabel('Real'); ylabel('Imag')
        patd = bintodec(pat);
        for k=0:p.alphabet-1
            ik = find(patd==k,1); % search if pattern is present
            plot(real(dataid(ik)),imag(dataid(ik)),'bo','Markerfacecolor','b')
            if strcmpi(x.plot,'bin')
                text(real(dataid(ik))+0.03,imag(dataid(ik)),...
                    num2cell(char(pat(ik,:)+48),2));
            elseif strcmpi(x.plot,'dec')
                text(real(dataid(ik))+0.03,imag(dataid(ik)),...
                    num2str(patd(ik,:)));
            else
                error('field plot must be ''bin'' or ''dec''.');
            end
        end
    end
end

%--------------------------------------------------------------------------
function dataid = finalizepat2samp(pat, p)

% Note: dataid is forced to have average energy 1 with a uniformly
%       distributed pattern.

[nrow,ncol] = size(pat);
M = 2^ncol; % constellation size

if strcmpi(p.family,'ook')
    dataid = 2*pat; % average energy: 1
    
elseif strcmpi(p.format,'bpsk') || strcmpi(p.format,'dpsk') || ...
        strcmpi(p.format,'psbt') || (strcmpi(p.family,'psk') && M==2)
    dataid = 2*pat-1; % average energy: 1
    
elseif strcmpi(p.format,'qpsk') || strcmpi(p.format,'dqpsk')
    level = 2*pat-1; % drive iqmodulator with QPSK
    dataid = complex(level(:,1),level(:,2))/sqrt(2); % average energy: 1
    
elseif strcmpi(p.family,'psk')
    pat = graytobin(pat);
    pat = [pat(:,1),xor(pat(:,1),pat(:,2)),pat(:,3:end)]; % restore Gray for first two (QPSK) symbols
    level = 2*pat - 1; % drive iqmodulator with QPSK
    dataid = complex(level(:,1),level(:,2))/sqrt(2); % QPSK
    for k=3:ncol % >=3 -> PM modulators for M-PSK, M>=8
        level(:,3:ncol) = pat(:,3:ncol); % drive pm_modulator(s)
        dataid=dataid.*fastexp(pi/2^(k-1)*level(:,k)); % M-PSK
    end
    
elseif strcmpi(p.family,'pam')
    pat = graytobin(pat);
    pat_dec = bintodec(pat);
    M = 2^ncol;
    dataid = 2*mod(pat_dec,M)-(M-1); % distance neighboring symbols: 2
    Es = (M^2-1)/6*2;
    dataid = dataid/sqrt(Es); % average energy: 1
    
elseif strcmpi(p.format,'qamrect') && M == 8
    patdec = bintodec(pat); % decimal
    % Look-up Table for 8-QAM rectangular
    LEVEL_data = [-3+1i; -3-1i; -1+1i; -1-1i; 3+1i; 3-1i; 1+1i; 1-1i]; % 8-QAM rect
    % LEVEL_data corresponds to the decimal, Gray, pattern 0:7
    dataid = LEVEL_data(patdec+1); % distance neighboring symbols: 2
    Es = 6;
    dataid = dataid/sqrt(Es); % average energy: 1
    
elseif strcmpi(p.format,'qamstar') && M == 8
    patdec = bintodec(pat); % decimal
    % Look-up Table for 8-QAM star
    ps3 = 1+sqrt(3);
    LEVEL_data = [-1-1i; ps3; -1i*ps3; 1-1i; -ps3; 1+1i; -1+1i; 1i*ps3]; % 8qam_star
    % LEVEL_data corresponds to the decimal, non Gray, pattern 0:7
    dataid = LEVEL_data(patdec+1); % distance neighboring symbols: 2
    Es = (2 + ps3^2)/2;
    dataid = dataid/sqrt(Es); % average energy: 1
    
elseif strcmpi(p.format,'qamcirc') && M == 8
    patdec = bintodec(pat); % decimal
    % Look-up Table for 8-QAM circular
    ALEV = [0; 3; 1; 2; 5; 4; 6];
    LEVEL_data = [cos(2*pi/7*ALEV)+1i*sin(2*pi/7*ALEV); 0]; % 8qam_circ
    % LEVEL_data corresponds to the decimal, non Gray, pattern 0:7
    dataid = LEVEL_data(patdec+1);
    Es = 7/8;
    dataid = dataid/sqrt(Es); % average energy: 1
    
elseif strcmpi(p.family,'qam') && M==32 % Cross-QAM
    patdec = bintodec(pat); % decimal
    % Look-up Table for 32-QAM
    LEVEL_data = [-3+5i; -1+5i; -3-5i; -1-5i; -5+3i; -5+1i; -5-3i; -5-1i; ...
        -1+3i; -1+1i; -1-3i; -1-1i; -3+3i; -3+1i; -3-3i; -3-1i; ...
        3+5i;  1+5i;  3-5i;  1-5i;  5+3i;  5+1i;  5-3i;  5-1i; ...
        1+3i;  1+1i;  1-3i;  1-1i;  3+3i;  3+1i;  3-3i;  3-1i]; % 32-QAM
    % LEVEL_data corresponds to the decimal, non Gray, pattern 0:31
    dataid = LEVEL_data(patdec+1); % distance neighboring symbols: 2
    Es = 20;
    dataid = dataid/sqrt(Es); % average energy: 1
    
elseif strcmpi(p.family,'qam') && M==128 % Cross-QAM
    patdec = bintodec(pat); % decimal
    % Look-up Table for 128-QAM
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
    % LEVEL_data corresponds to the decimal, non Gray, pattern 0:127
    dataid = LEVEL_data(patdec+1); % distance neighboring symbols: 2
    Es = mean(abs(LEVEL_data).^2);
    dataid = dataid/sqrt(Es); % average energy: 1
    
elseif strcmpi(p.family,'qam') % Square QAM
    pat = [graytobin(pat(:,1:end/2)),graytobin(pat(:,end/2+1:end))];
    % first half: real. Second half: imag.
    sqrtM = sqrt(M); % (square) QAM side length
    if rem(sqrtM,2) ~= 0, error('Unknown QAM');end
    ncolh = ncol/2;
    level = zeros(nrow,2);
    for m=1:2 % real/imag components
        pat_dec = bintodec(pat(:,1+(m-1)*ncolh:m*ncolh));
        level(:,m) = 2*mod(pat_dec,sqrtM)-(sqrtM-1);
    end
    dataid = complex(level(:,1),level(:,2)); % distance neighboring symbols: 2
    Es = (M-1)/6*4;
    dataid = dataid/sqrt(Es); % average energy: 1
    
else
    error('Unknown modulation format');
end


