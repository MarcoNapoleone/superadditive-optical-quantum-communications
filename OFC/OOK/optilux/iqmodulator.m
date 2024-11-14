function E=iqmodulator(E, modsig, options)

%IQMODULATOR in-phase/quadrature optical modulator
%   EOUT=IQMODULATOR(E,MODSIG, OPTIONS) modulates the in-phase and 
%   quadrature components of the field E (see LASERSOURCE) with
%   Mach-Zehnder modulators [1].
%
%           ____ MZ (in-phase) ____
%         /                         \
%   ------                           ------
%         \ ____ MZ (quad) _ pi/2 _ /
%
%   MODSIG is a complex vector containing in the real and in the imaginary
%   parts the in-phase and the quadrature components, respectively. Such
%   driving signals are created, e.g., by DIGITALMOD.
%
%   OPTIONS is an optional structure whose fields can be:
%
%       - iqratio : power ratio [dB] between I and Q arms (default = 0)
%       - biasc : normalized bias between I and Q arms. The bias is
%           normalized to the voltage yielding a phase shift of pi in
%           both arms, usually called Vpi. (default = 0)
%       - exratio : extinction ratio of the two nested modulators [dB]
%           (default = [inf inf])
%       - bias  : bias of the two nested modulators. See MZMODULATOR.
%           (default = [-1 -1]).
%       - biasl: similar meaning of bias, but for the lower arm only of
%           each Mach-Zehnder modulator. It is ignored if OPTIONS.bias
%           exists.
%       - biasu: similar meaning of bias, but for the upper arm only of
%           each Mach-Zehnder modulator. It is ignored if OPTIONS.bias
%           exists.
%       - Vpi: voltage of pi phase shift in the MZ modulators.
%       - nch: channel index. This option must be used when E.field is a
%           matrix and represents the column on which the function
%           operates.
%
%   See also PATTERN, DIGITALMOD, LASERSOURCE, MZMODULATOR.
%
%   [1] R. A. Griffin, "Integrated DQPSK transmitters," in Proc. OFC 2005,
%   2005 Anaheim, CA, USA
%
%   [2] K. Hoon and A. H. Gnauck, "Chirp characteristics of dual-drive
%   Mach-Zehnder modulator with a finite DC extinction ratio," IEEE Photon.
%   Technol. Lett., vol. 14, pp. 298-300, Mar. 2002.
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

iqratio     = 0; % [dB]
biasc       = 0;
biasl       = [-1 -1];
biasu       = [-1 -1];
Vpi         = pi/2*[1 1];
exratio     = [inf inf]; % [dB]
mode        = 'push-pull';
normf       = 1;

if length(E.lambda) > 1
    if nargin < 3 || ~isfield(options,'nch')
        error('Missing channel index');
    end
else
    nch = 1;         % number of channels
end

if size(E.field,2) == 2*length(E.lambda)% dual polarization
    Npol = 2;
else
    Npol = 1;
end

if nargin > 2
    if isfield(options,'iqratio')
        iqratio = options.iqratio;
    end
    if isfield(options,'biasc')
        biasc = options.biasc;
    end
    if isfield(options,'biasl')
        biasl = options.biasl;
        if length(biasl) ~= 2
            error('biasl must be a vector of length 2')
        end
    end
    if isfield(options,'biasu')
        biasu = options.biasu;
        if length(biasu) ~= 2
            error('biasu must be a vector of length 2')
        end
    end
    if isfield(options,'bias')
        biasu = options.bias;
        biasl = options.bias;
        if length(bias) ~= 2
            error('bias must be a vector of length 2')
        end
    end
    if isfield(options,'Vpi')
        Vpi = options.Vpi;
        if length(Vpi) ~= 2
            error('Vpi must be a vector of length 2')
        end
    end
    if isfield(options,'exratio')
        exratio = options.exratio;
        if length(exratio) ~= 2
            exratio = [exratio exratio];
        end
    end
    if isfield(options,'nch')
        nch = options.nch;
    end
    if isfield(options,'norm')
        normf = options.norm;
    end
end

iqratio = 10^(iqratio/20);
sr = iqratio/(1+iqratio);

% Now set polarizations in alternate way (if they exist)
ncols = (Npol*(nch-1)+1):(Npol*nch);
% E.g. Npol = 2, nch = 3 -> ncols = [5 6]
% E.g. Npol = 1, nch = 3 -> ncols = 3

%   The main function

Ei = E; Eq = E;
Ei.field(:,ncols) = E.field(:,ncols) * sr;         % In Phase Field
Eq.field(:,ncols) = E.field(:,ncols) * (1-sr);     % Quadrature Field

Ei = mzmodulator(Ei,real(modsig),struct('exratio',exratio(1),'biasl',...
    biasl(1),'biasu',biasu(1),'Vpi',Vpi(1),'mode',mode,'nch',nch,'norm',normf));
Eq = mzmodulator(Eq,imag(modsig),struct('exratio',exratio(2),'biasl',...
    biasl(2),'biasu',biasu(2),'Vpi',Vpi(2),'mode',mode,'nch',nch,'norm',normf));

% *2: undo the two 3dB coupler loss
E.field(:,ncols) = 2*(Ei.field(:,ncols) + Eq.field(:,ncols)*fastexp(pi/2 + biasc));

