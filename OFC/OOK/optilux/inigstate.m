function inigstate(Nsamp,Fs)

%INIGSTATE Initializes the global variable GSTATE of OptiluX.
%
%   INIGSTATE(NSAMP) initializes the global struct variable GSTATE of
%   OptiluX by setting the number of samples equal to NSAMP. NSAMP will be
%   the reference values used by the OptiluX signals. This call should be 
%   used for digital systems working at one sample per symbol.
%
%   When more than one sample per symbol is used, the initialization is by
%   INIGSTATE(NSAMP,FS), FS being the sampling frequency FS [GHz], also
%   called Nyquist frequency. Similarly, FS will be the reference
%   sampling frequency for all OptiluX signals.
%
%   Such variables become fields of the global struct variable GSTATE:
%
%       GSTATE.NSAMP = NSAMP
%       GSTATE.FSAMPLING = FS   % [GHz]
%
%   INIGSTATE initializes also the discrete frequencies in GSTATE.FN:
%
%       GSTATE.FN = fftshift( -FS/2 : FS/NSAMP : FS/2 - FS/NSAMP); % [GHz]
%
%   Some important notes:
%
%       * GSTATE.FN is in [GHz]. The normalized frequencies as in digital
%           signal processing textbooks are GSTATE.FN/FS.
%       * The LOWEST discrete frequency (resolution) is FS/NSAMP.
%       * The LARGEST discrete frequency (Nyquist frequency) in absolute
%           value is FS/2.
%
%   Please remember that OptiluX works with FFT, hence the discrete
%   frequency axis wraps at multiples of the sampling rate FS. For
%   instance, the next discrete frequency on the right of (FS/2 - FS/NSAMP)
%   is -FS/2.
%
%   INIGSTATE initializes also the global struct variable CONSTANTS
%   containing fundamentals physical constants (units of the international
%   system):
%
%       CONSTANTS.HPLANCK:      Planck's constant [J*s]
%       CONSTANTS.CLIGHT:       speed of light in vacuum [m/s]
%       CONSTANTS.ECHARGE:      electron charge [C]
%       CONSTANTS.KBOLTZMANN:   Boltzmann's constant [J/K]
%
%   See also: DIGITALMOD, PATTERN
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

global GSTATE;  % GSTATE is a global structure variable.
global CONSTANTS;  % CONSTANTS is a global structure variable.

CONSTANTS.CLIGHT= 299792458;        % speed of light in vacuum [m/s]
CONSTANTS.HPLANCK= 6.62606896e-34;  % Planck's constant [J*s] (CODATA value,
% std. uncertainty= 3.3e-43). NOTE: in October 2005, the National Physical
% Laboratory reports 6.62607095e-34.
CONSTANTS.ECHARGE= 1.602176487e-19;  % Electron's charge [C] (CODATA value, year 2006)
CONSTANTS.KBOLTZMANN= 1.3806504E-23; % Boltzmann's constant [J/oK] (CODATA value,
% std. uncertainty= 2.4e-31)

if ~(floor(Nsamp) == Nsamp) % not integer
    error('The number of samples must be an integer');
end

GSTATE.NSAMP = Nsamp;           % number of samples
if nargin == 2
    stepf = Fs/Nsamp; % minimum frequency [GHz]
    GSTATE.FN = fftshift(-Fs/2:stepf:Fs/2-stepf);   % frequencies [GHz]
    GSTATE.FSAMPLING = Fs;        % sampling frequency [GHz]
end
