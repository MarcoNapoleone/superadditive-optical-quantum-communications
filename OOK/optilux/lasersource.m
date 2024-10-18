function E=lasersource(Ptx,lam,o1,o2,o3)

%LASERSOURCE Multi-channel laser source
%   E = LASERSOURCE(PTX,LAM) creates a comb of constant waves of powers
%   PTX [mW] at wavelengths LAM [nm]. Such waves are stored in the struct E
%   of fields:
%
%       E.lambda = LAM
%       E.field  = time samples of the electric field along rows
%
%   PTX can be a scalar or a vector of the same length of LAM. In the first
%   case, the same power is used for all channels.
%
%   In dual-polarization mode (default), E.field has 2*length(LAM) columns, 
%   where the odd columns represent x-polarizations samples, while 
%   the even columns the y-polarizations samples. In single-polarization 
%   mode (see options), the length is halved and only the x-polarizations 
%   exist along columns.
%
%   For instance, the call E = LASERSOURCE(5,[1550 1551]) yields:
%
%       E.lambda = [1550 1551]
%       E.field(k,:) = [sqrt(5), 0, sqrt(5), 0] 
%
%   where E.field(k,1) refers to carrier 1550, x-polarization, E.field(k,2)
%   to carrier 1550, y-polarization, E.field(k,3) to carrier 1551,
%   x-polarization, E.field(k,4) to carrier 1551, y-polarization.
%
%   E = LASERSOURCE(PTX,LAM,SPAC,NLAMBDA) works with a scalar LAM and the
%   additional input SPAC representing the channel spacing [nm]. In such a
%   case, a comb of uniformly spaced constant waves is created up to
%   NLAMBDA carriers, with central wavelength LAM.
%
%   Note: if SPAC is present, the channel spacing is uniform in the
%       frequency domain, according to the ITU-T recommendations. Hence,
%       the spacing is slightly non-uniform in the wavelength domain. For
%       instance, with LAM = 1550, SPAC = 0.4 and NLAMBDA = 5, the 
%       following wavelengths, and the corresponding carrier frequencies, 
%       are created (only four digits are shown):
%
%       [1549.2004 1549.6001 1550.0000 1550.4001 1550.8004]     [nm]
%       [193.3147  193.3646  193.4145  193.4644  193.5143]      [THz]
%
%
%   E = LASERSOURCE(PTX,LAM,SPAC,NLAMBDA,OPTIONS) or
%   E = LASERSOURCE(PTX,LAM,OPTIONS) has the additional struct variable
%   OPTIONS:
%
%      OPTIONS.pol: if 'single' only the x-polarization is created,
%           otherwise even the y-polarization is created and set equal to
%           zero in absence of noise. This option is useful to work in a
%           purely scalar regime.
%      OPTIONS.linewidth: the 3 dB two-sided linewidth [GHz] of the laser.
%           It can be a scalar or a vector of the same length of the
%           wavelengths. A Wiener phase noise with such a linewidth is
%           added to the phase of  E. The function uses the Brownian bridge
%           trick [Wiki].
%      OPTIONS.n0: the one-sided spectral density [dB/GHz] of a Gaussian
%           complex noise added to the laser field. The noise variance,
%           i.e., the sum of the real and the imaginary components
%           variance, is:
%
%               variance = 10^(OPTIONS.n0/10)*GSTATE.FSAMPLING
%
%           Under a small-signal assumption, such a noise is related to the
%           laser relative-intensity noise (RIN), on a dB scale, by:
%
%               RIN = 0.5*(3 + OPTIONS.n0 - 10*log10(PTX))
%
%
%   See also MZMODULATOR, DIGITALMOD, PATTERN, INIGSTATE
%
%   References:
%
%   [Wiki] https://en.wikipedia.org/wiki/Brownian_bridge
%
%   Author: Paolo Serena, 2021
%           Massimiliano Salsi, 2009
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


global GSTATE   % GSTATE is a structure whose fields are defined in inigstate.m

if isempty(GSTATE)
    error('The simulation must start with inigstate.m');
end
Nsamp = GSTATE.NSAMP;
Npow = length(Ptx);

spac = 0;
if nargin == 2 % (power, lambda)
    Nch = length(lam);   % number of carriers
    
elseif nargin == 3 % (power, lambda, options)
    Nch = length(lam);   % number of carriers
    if isstruct(o1)
        options = o1;
    else
        error('wrong input variables type. Missing number of channels?');
    end
    
else % nargin > 3
    if ~isstruct(o1) % (power, lambda, spac, Nch)
        if length(lam) > 1
            error('LAM is the central wavelength of the WDM comb.');
        end
        spac = o1; % wavelength spacing [nm]
        if ~(floor(o2) == o2)
            error('The number of channels must be an integer.');
        else
            Nch = o2; % number of carriers
        end
    else
        error('wrong input variables type');
    end
    if nargin > 4 % (power, lambda, spac, Nch, options)
        if ~isstruct(o3)
            error('The options must be a struct variable');
        end
        options = o3;
    end
end

if Npow > 1 && Npow ~= Nch
    error('Powers and wavelengths must have the same length');
end

if exist('options','var')
    checkfields(options,{'n0','linewidth','pol'});
    if isfield(options,'n0')
        N0 = options.n0;
        if length(N0) > 1, error('N0 length must be 1'); end
    else
        N0 = -inf;
    end
    if isfield(options,'pol') && strcmpi(options.pol,'single')
        Npol = 1; % scalar propagation
    else
        Npol = 2; % dual-polarization propagation
    end    
    if isfield(options,'linewidth')
        linewidth = options.linewidth;
        Llnw = length(linewidth);
        if Llnw == 1
            linewidth = linewidth*ones(1,Nch);
        elseif Llnw ~= Nch
            error('The linewidth must have the same length of the number of wavelengths.')
        end
    else
        linewidth=0;
    end
else
    
    linewidth = 0;
    N0 = -Inf;
    Npol = 2;
end     % end of options checks

if Npow == 1
    power = repmat(Ptx(:).',1,Nch);
else
    power = Ptx(:).';
end

if spac ~= 0 % uniformly spaced carriers
    E.lambda = getlambda(lam,spac,Nch);
else
    E.lambda = lam(:).'; % force row
end

E.field = zeros(Nsamp,Nch*Npol);
% by default, fully polarized on x (odd columns):
E.field(:,1:Npol:Nch*Npol) = repmat(sqrt(power),Nsamp,1);

%%% Add phase noise
if any(linewidth)
    freq_noise  = (ones(Nsamp,1) * sqrt(2*pi*linewidth./GSTATE.FSAMPLING)) .* ...
        randn( Nsamp, Nch);
    freq_noise(1) = 0;
    phase_noise = cumsum(freq_noise,1);
    % Brownian bridge [Wiki]
    tim = 0:Nsamp-1;
    phase_noise = phase_noise - tim(:)/(Nsamp-1)*phase_noise(end,:);
    %     E1(:,:,k) = repmat(E(k,:),Nsamp,1) .* fastexp(phase_noise);
    E.field(:,1:Npol:Nch*Npol) = fastexp(phase_noise) .* E.field(:,1:Npol:Nch*Npol);
end

% Add Gaussian complex white noise (ASE)
if ~isinf(N0)
    N0_lin=10^(N0/10);
    sigma = sqrt(N0_lin/2*GSTATE.FSAMPLING); % noise std
    E.field = E.field + sigma*complex(randn(Nsamp,Nch*Npol),randn(Nsamp,Nch*Npol));
end
