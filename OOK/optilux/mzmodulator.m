function E=mzmodulator( E, modsig, options)
%MZMODULATOR Mach-Zehnder modulator
%   E = MZMODULATOR(E, MODSIG, OPTIONS) modulates the optical field E with  
%   the electric signal MODSIG by a Mach-Zehnder interferometer [1].
%   E is a struct of fields lambda, field (see LASERSOURCE).
%
%   MODSIG is the electrical driving signal (e.g., produced by DIGITALMOD).
%
%   OPTIONS is an optional struct variable with fields:
%
%       - exratio : extinction ratio [dB] (default = inf)
%       - biasl: bias, normalized to Vpi, of the lower arm only. It is
%           ignored if OPTIONS.bias exists. Default = -1
%       - biasu: similar meaning of bials, but for the upper arm only.
%       - mode: operating mode. It can be 'push-push' where the same MODSIG
%           is applied to both arms, or 'push-pull' where a phase shift of
%           pi is induced in the lower arm. Default: 'push-pull'.
%       - nch: channel (wavelength) index. This option must be used when 
%           E.field is a matrix and represents the column on which the 
%           function operates. In dual-polarization the function operates 
%           on columns [nch (nch+1)] representing x and y polarizations.
%       - Vpi: voltage yielding a phase shift of pi in each MZ branch.
%           Default: pi/2 in push-pull mode.
%
%
%   See also LASERSOURCE, DIGITALMOD, IQMODULATOR.
%
%   [1] S. Walklin and J. Conradi, "Effect of Mach-Zehnder modulator DC
%   extinction ratio on residual chirp-induced dispersion in 10-Gb/s binary
%   and AM-PSK duobinary lightwave systems," IEEE Photon. Technol. Lett.,
%   vol. 9, pp. 1400-1402, Oct. 1997.
%   [2] K. Hoon and A. H. Gnauck, "Chirp characteristics of dual-drive
%   Mach-Zehnder modulator with a finite DC extinction ratio," IEEE Photon.
%   Technol. Lett., vol. 14, pp. 298-300, Mar. 2002.
%
%   Author: Paolo Serena, 2021
%			Massimiliano Salsi, 2009
%			Marco Bertolini, 2009
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

% Default values
biasl      = -1;    % bias of lower arm
biasu      = -1;    % bias of upper arm
exratio    = inf;   % extinction ratio
pp = 0;             % 0: push-pull. 1: push-push
Vpi = pi/2*(pp+1);  % voltage of phase shift pi
normf = 1;          % normalization factor

if length(E.lambda) > 1 % need to know which channel to modulate
    if nargin < 3 || ~isfield(options,'nch')
        error('Missing channel index');
    end
else
    nch = 1;  % number of channels (wavelengths)
end

if nargin > 2
    checkfields(options,{'Vpi','exratio','biasu','biasl','bias',...
        'mode','nch','norm'});
    if isfield(options,'Vpi')
        Vpi = options.Vpi;
    end
    if isfield(options,'exratio')
        exratio = options.exratio;
    end
    if isfield(options,'biasu')
        biasu = options.biasu;
    end
    if isfield(options,'biasl')
        biasl = options.biasl;
    end
    if isfield(options,'bias')
        biasu = options.bias;
        biasl = options.bias;
    end
    if isfield(options,'mode')
        if strcmpi(options.mode,'push-push')
            pp = 1;
        elseif strcmpi(options.mode,'push-pull')
            pp = 0;
        else
            error('Mode can be ''push-push'' or ''push-pull''');
        end
    end
    if isfield(options,'nch')
        nch = options.nch;
        if nch > length(E.lambda)
            error('Channel index does not exist');
        end
    end
    if isfield(options,'norm')
        normf = options.norm;
        if normf == 1
            warning('optilux:mzmodulator','normalization factor not used.')
        end
    end
end

% Set-up extinction ratio
invexr_lin = 10^(-exratio/10);
gamma = (1-sqrt(invexr_lin))/(1+sqrt(invexr_lin)); % [1]

modsig = real(modsig); % signal must be real
if pp == 1 % push-push
    Phi_U = pi * (modsig + biasu*Vpi)/Vpi;
    Phi_L = pi * (modsig + biasl*Vpi)/Vpi;
else % push-pull
    Phi_U = pi/2 * (modsig + biasu*Vpi)/Vpi;
    Phi_L = -pi/2 * (modsig + biasl*Vpi)/Vpi;
end

if size(E.field,2) == 2*length(E.lambda)% dual polarization
    Npol = 2;
else
    Npol = 1;
end

% Now set polarizations in alternate way (if they exist)
ncols = (Npol*(nch-1)+1):(Npol*nch);
% E.g. Npol = 2, nch = 3 -> ncols = [5 6]
% E.g. Npol = 1, nch = 3 -> ncols = 3


% now modulation only on the existing polarizations
for m=1:Npol
    np = ncols(m);
    if any(E.field(:,np))
        E.field(:,np)  = normf*E.field(:,np) .* (fastexp(Phi_U) + ...
            gamma*fastexp(Phi_L))/(1+gamma);
        % with default values -> E.field = E.field .* sin(modsig)
    end
end
