function Iric = rxfrontend(E,lam,symbrate,x)

%RXFRONTEND receiver front-end.
%   IRIC=RXFRONTEND(E,LAM,SYMBRATE,X) returns the received
%   current of a direct-detection or coherent receiver. Direct-detection is
%   used with on-off keying (OOK), differential-phase shift keying (DPSK)
%   and differential-quadrature phase shift keying (DQPSK). In all the
%   other cases coherent detection is used.
%
%   The block implements the following operations:
%
%          __
%         /  \
%        |    |      -------        ----------
%  E      \  /      |       |      | OPTICAL  |
%     --------------| OBPF  |----- |   TO     |-----IRIC
%         post      |       |      |ELECTRICAL|
%         fiber      -------        ----------
%
%
%
%   A more detailed description for each modulation format can be found
%   below.
%
%   E is a struct of fields:
%
%       E.lambda: central wavelength [nm] of the electric field, i.e.,
%           wavelength that relates the lowpass equivalent signal to the
%           corresponding bandpass signal.
%       E.field: time samples of the electric field, with polarizations (if
%           existing) alternated on columns
%
%   SYMBRATE is the symbol rate [Gbaud] of the signal. LAM [nm] is the
%   central wavelength of the optical bandpass filter. LAM may differ from
%   E.lambda for instance with wavelength division multiplexing (WDM).
%
%   X is a structure of fields:
%
%       X.oftype = optical filter (OBPF) type (see MYFILTER)
%       X.obw = OBPF bandwidth normalized to SYMBRATE.
%       X.opar = optical filter extra parameters (optional, see MYFILTER)
%
%   Other optional parameters of X:
%
%       X.dcum = post compensating fiber accumulated dispersion [ps/nm].
%           The fiber is inserted before the optical filter.
%       X.slopecum = post compensating fiber accumulated slope [ps/nm^2].
%       X.lambda = wavelength [nm] at which the post compensating fiber
%              has an accumulated dispersion equal to X.dcum.
%
%   The post-compensating fiber is assumed to be a purely linear lossless
%   fiber, while the photodiodes are ideal (abs(.)^2). For this reason, it
%   is not important the fiber dispersion coefficient and its length, but
%   only their product X.dcum.
%
%   IRIC is a vector containing the received current [mA].
%
%
%   The receiver for OOK modulation is the following:
%          __
%         /  \                      photodiode
%        |    |        -------          |
%  E      \  /        |       |       __|__
%     ----------------| OBPF  |------  / \ --------IRIC
%         post        |       |       /   \
%         fiber        -------        -----
%                                       |
%                                       |
%
%   The receiver for DPSK modulation is the following:
%
%          __                 Mach Zehnder (MZ)   __|__
%         /  \                    ---------        / \
%        |    |     -------     /          \   /   ---
%  E      \  /     |       |   /            \ /     |
%     -------------| OBPF  |---              -      |-----IRIC
%         post     |       |   \   ------   / \   __|__
%         fiber     -------     \_| 1 bit|_/   \   / \
%                                 | delay|         ---
%                                  ------           |
%
%   The receiver for DQPSK modulation is the following:
%
%
%                             Mach Zehnder (MZ)   __|__
%                                 ---------        / \
%                                /  -pi/4   \   /  ---
%                               /            \ /    |
%                             -               -     |-----real(IRIC)
%          __                /  \   ------   / \  __|__
%         /  \              /    \_| 1 bit|_/   \  / \
%        |    |  -------   /       | delay|        ---
%  E      \  /  |       | /         ------          |
%    -----------| OBPF  |/
%         post  |       |\
%         fiber  -------  \   Mach Zehnder (MZ)   __|__
%                          \      ---------        / \
%                           \    /  +pi/4   \   /  ---
%               		     \  /            \ /    |
%                              -              -     |-----imag(IRIC)
%                               \   ------   / \  __|__
%                                \_| 1 bit|_/   \  / \
%                                  | delay|        ---
%                                   ------          |
%
%   The receiver for coherent detection is the following:
%
%          __
%         /  \
%        |    |        -------      -----
%  E      \  /        |       |    |     |
%     ----------------| OBPF  |----|  *  |-----IRIC
%         post        |       |     -----
%         fiber        -------        |
%                                     |
%                              oscillator @ lambda
%
%   See also MYFILTER, RXDSP.
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



%%% INIT AND CHECKS
if nargin < 4
    error('Not enough input arguments');
end

global GSTATE   % GSTATE is a structure whose fields are defined in inigstate.m

if ~isfield(x,'opar')
    x.opar = 0;     % dummy value
end

if ~isfield(x,'mzdel')
    x.mzdel = 1;     % default interferometer delay: 1 symbol
elseif (x.mzdel <= 0) || (x.mzdel > 1)
    error('The delay of the interferometer must be  0 < mzdel <= 1')
end

if ~isfield(x,'modformat')
    error('Missing modulation format');
end

% create linear optical filters: OBPF (+fiber)
Fnorm = GSTATE.FN/symbrate;
if isfield(x,'dcum')
    Hf = postfiber(x.dcum,x.slopecum,x.lambda,lam,E.lambda);
    Hf = Hf .* myfilter(x.oftype,Fnorm,0.5*x.obw,x.opar);
else
    Hf = myfilter(x.oftype,Fnorm,0.5*x.obw,x.opar);
end

%%% 1: apply optical filter

E = filterenv(E,lam,Hf);

%%% 2: optical to electrical conversion

Nt = GSTATE.FSAMPLING/symbrate; % number of points per symbol
Iric = o2e(E,Nt,x);

%--------------------------------------------------------------------------
function Hf = postfiber(dcum,slopecum,lam0,lam,lamc)

%POSTFIBER post-compensating single-mode fiber
%   Hf = POSTFIBER(DCUM,SLOPEZ,LAM0,LAM,LAMC) returns the frequency
%   response of a lossless and purely linear single-mode fiber that
%   accumulates a dispersion DCUM [ps/nm] and a third-order dispersion
%   SLOPEZ [ps^2/nm] at the wavelength LAM0 [nm]. LAMC [nm] is the central
%   wavelength of the bandpass electric field propagating in such a fiber.
%   LAM [nm] is the central wavelength of the optical filter.

global CONSTANTS
global GSTATE

CLIGHT = CONSTANTS.CLIGHT; % [m/s]
b20z = -lam0^2/2/pi/CLIGHT*dcum*1e-3; % beta2*z [ns^2] @ lam0
b30z = (lam0/2/pi/CLIGHT)^2*(2*lam0*dcum+...
    lam0^2*slopecum)*1e-3; % beta3*z [ns^3] @ lam0

%   lam:  wavelength of the channel's carrier under detection
%   lam0: wavelength @ fiber parameters
%   lamc: wavelength of central frequency of bandpass signal

% Domega_ik: [1/ns]. "i" -> at ch. i, "0" -> at lam0
Domega_i0 = 2*pi*CLIGHT*(1./lam-1/lam0);
Domega_ic = 2*pi*CLIGHT*(1./lam-1/lamc);
Domega_c0 = 2*pi*CLIGHT*(1./lamc-1/lam0);
beta1z = b20z*Domega_ic+0.5*b30z*(Domega_i0^2-Domega_c0^2);    %[ns]
beta2z = b20z+b30z*Domega_i0;  % beta2*z [ns^2]@ lam
% dispersion of the channels
omega = 2*pi*GSTATE.FN'; % angular frequency [rad/ns]
betat = omega*beta1z+0.5*omega.^2*beta2z+omega.^3*b30z/6;

Hf = fastexp(-betat);

%--------------------------------------------------------------------------
function E = filterenv(E,lam,Hf)

%FILTERENV filter around a given wavelength
%   E = FILTERENV(E,LAM,HF) first aligns the channel's carrier
%   wavelength LAM to the central wavelength of the optical filter of
%   frequency response HF, then applies the filter.

global CONSTANTS
global GSTATE

freqc = CONSTANTS.CLIGHT./E.lambda; % central frequency [GHz] (corresponding to the zero
%   frequency of the lowpass equivalent signal by convention
freq = CONSTANTS.CLIGHT./lam; % carrier frequency [GHz]
deltafn = freqc - freq;   % carrier frequency spacing [GHz]
minfreq = GSTATE.FN(2)-GSTATE.FN(1);    % resolution [GHz]
ndfn = round(deltafn./minfreq);  % spacing in points

E.field = fft(E.field);
E.field = fastshift(E.field,ndfn); % undo what did in MULTIPLEXER
E.field = ifft( bsxfun(@times, Hf , E.field) );

%--------------------------------------------------------------------------
function Iric = o2e(E,Nt,x)

Nfft = length(E.field);

if strcmpi(x.modformat,'ook')
    Iric = sum(abs(E.field).^2,2); % PD. sum is over polarizations
elseif strcmpi(x.modformat,'dpsk')
    ndel = nmod((1:Nfft)-round(x.mzdel*Nt),Nfft); % interferometer delay
    
    Iric = sum( real(E.field .* conj(E.field(ndel))) ,2); % MZI + PD
elseif strcmpi(x.modformat,'dqpsk')
    ndel = nmod((1:Nfft)-round(x.mzdel*Nt),Nfft); % interferometer delay
    
    Iric = sum( fastexp(-pi/4)*real(E.field .* conj(E.field(ndel))) ,2) + ...
        1i*sum( fastexp(pi/4)*real(E.field .* conj(E.field(ndel))) ,2);
else % coherent detection
    Iric = E.field;
end

