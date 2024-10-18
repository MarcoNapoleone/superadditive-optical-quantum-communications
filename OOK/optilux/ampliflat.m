function [E,out] = ampliflat(E,x)

%AMPLIFLAT optical amplifier with flat gain in wavelength
%   E = AMPLIFLAT(E,X) amplifies the optical field E. E is a struct of
%   fields E.lambda, E.field, see LASERSOURCE.
%
%   X contains information about the amplification mode:
%
%       X.gain: amplifier gain [dB] in constant gain mode. This option is
%           ignored if X.outpower exists.
%       X.f:    noise figure [dB]. If absent, the amplifier is noiseless.
%               Note: by convention, the noise figure is related to the
%               one-sided power spectral density N0 of amplified 
%               spontaenous emission (ASE)  by N0=h*nu*X.f*G, with h
%               Planck's constant, nu carrier frequency, G gain.
%       X.outpower: output power [dBm]. If present, the amplifier works in
%           constant-output power mode, with or without ASE. The main
%           difference with constant gain mode (default) is in the 
%           presence of amplified spontaneous emission (ASE) noise. In 
%           constant-output power mode the input signal experiences a droop 
%           to make room for the incoming ASE power. The gain is estimated 
%           by using an estimation of the average power of the input 
%           signal. Alternatively, the power can be provided by the user 
%           (see X.inpower).
%       X.inpower: input power [dBm]. By default the input power is
%           estimated from E. However, the estimation may be unrealistic,
%           especially with few samples. The user can manually set the
%           value of the average power in X.inpower. X.inpower is necessary
%           only in constant output power mode.
%       X.asebandwidth: If X.f exists, the ASE, by default, is generated
%           over the simulation bandwidth, set in INIGSTATE. The user can
%           locally set the ASE bandwidth by X.asebandwidth [GHz]: the ASE
%           is still white over the simulation bandwidth as in the default
%           case, but with a variance scaled to meet the requirements of 
%           X.asebandwidth.
%
%   [E,OUT]=AMPLIFLAT(E,X) also returns on output information about the
%   amplification:
%
%       OUT.gain: gain [dB] used by the amplifier (possibly with droop).
%       OUT.droop: droop factor (if used). linear scale.
%       OUT.asebandwidth: bandwidth [GHz] of ASE (if ASE is present).
%
%   See Also: FIBER.
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

global CONSTANTS;  % physical constants
global GSTATE;

HPLANCK = CONSTANTS.HPLANCK;    % Planck's constant [J*s]
CLIGHT = CONSTANTS.CLIGHT;      % speed of light [m/s]

%%% Init
if nargin < 2, error('Missing amplifier parameters'); end

if isfield(x,'outpower') && isfield(x,'gain')
    warning('optilux:ampliflat','Gain value ignored');
    x.gain = NaN;
end

%%% Noiseless Gain

if isfield(x,'outpower') % constant (average) output power
    if isfield(x,'inpower')
        Pi = 10^(x.inpower*0.1); % [mW]
    else
        Pi = sum(mean(abs(E.field).^2),2); % [mW]
    end
    Po = 10^(x.outpower*0.1); % [mW]
    gain = Po/Pi; % noiseless, estimated, gain. Updated with droop below.
else
    gain = 10^(x.gain*0.1); % noiseless gain
end

%%% ASE
if isfield(x,'f') % noisy amplifier
    if isfield(x,'asebandwidth')
        if x.asebandwidth < GSTATE.FSAMPLING
            warning('optilux:ampliflat','ASE bandwidth smaller than simulation bandwidth');
        end
        Bw = x.asebandwidth; % ASE bandwidth [GHz]
    else
        Bw = GSTATE.FSAMPLING; % bandwidth [GHz] of numerical simulation
    end
    Flin = 10^(x.f*0.1); % noise figure
    sigma = sqrt(Flin/4*HPLANCK*CLIGHT./E.lambda*gain.*Bw*1e21);
else
    sigma = 0; % sigma^2: ASE variance [mW]
end

%%% Droop
if isfield(x,'outpower')
    sigma2ASE = 4*sigma.^2; % Rex,Imx,Rey,Imy
    droop = Po/(Po + sigma2ASE); % droop factor
    gain = gain*droop;
    sigma = sigma*sqrt(droop);
end

%%% Amplification
E.field = E.field*sqrt(gain); % amplification

%%% Add ASE
if sigma ~= 0
    [nrow,ncol] = size(E.field);
    noise = sigma*(randn(nrow,ncol) + 1i*randn(nrow,ncol));
    E.field = E.field + noise;
end

%%% Info on output

if nargout == 2
    if isfield(x,'f'), out.asebandwidth = Bw; end % [GHz]
    out.gain = 10*log10(gain); % [dB]
    if isfield(x,'outpower'), out.droop = droop; end
end
