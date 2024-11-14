function Eo = multiplexer(E,o1,o2)

%MULTIPLEXER optical multiplexer
%   E=MULTIPLEXER(E) multiplexes the electric fields contained in E. in
%   particular, E is a struct variable with the following fields:
%
%       E.lambda: carrier wavelengths [nm] of the input fields, identical
%           for both x and y polarizations.
%       E.field:  electric field time samples. E.field is a
%           matrix of size [Nsamp,Ncol], with Nsamp the number of samples
%           (see INIGSTATE), while Ncol=2*length(E.lambda) in dual
%           polarization or length(E.lambda) in scalar propagation. Such
%           columns will be multiplexed together.
%
%   For instance, E.field columns may experience the following mixing:
%
%       E.field(:,1)   ---------------------
%                                            \
%                                             \
%                                              -------- E.field
%                                             /
%                                            /
%       E.field(:,end) ---------------------
%
%   Since E.field is the lowpass equivalent signal of the electric field,
%   each column is first modulated at the proper wavelength before
%   combining the signals. The zero frequency of the resulting signal is,
%   by convention, equal to the barycenter of the carrier frequencies of
%   the input signals.
%
%   E=MULTIPLEXER(E,OPTIONS) accepts the optional struct input OPTIONS,
%   whose fields are:
%
%       OPTIONS.lambda = central wavelength [nm] of the resulting
%           multiplexed signal. If not specified, the central frequency is
%           equal to the barycenter of the two extreme carriers of the
%           fields, hence the central wavelength is its inverse times the
%           speed of light.
%       OPTIONS.mux: if equal to 'sepfields', the wavelengths are not
%           multiplexed together but just concatenated along columns of
%           E.field.
%       OPTIONS.delay: can be 'rand', which adds a random delay to each
%           column of E.field before multiplexing. An independent delay is
%           added to each polarization.
%           Otherwise, OPTIONS.DELAY can be a vector of integers, of length
%           equal to the number of polarizations, containing the delays in
%           samples.
%
%   E=MULTIPLEXER(E,E1) or E=MULTIPLEXER(E,E1,OPTIONS) operates as before
%   but concatenates the electric fields contained in the struct E and E1
%   (E1 has similar fields of E). Such calls are useful to multiplex
%   several signals by iterative calls to MUTIPLEXER.
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

%%% Initial checks
if ~isstruct(E)
    error('The input electric field must be a struct variable.');
end
options = [];
checkfields(E,{'lambda','field'});

[Elambda,Npol] = expandlambda(E);
if nargin >= 2
    if ~isstruct(o1)
        error('All input parameters must be struct variables.');
    end
    if isfield(o1,'field') % Concatenate the structs E and o1
        checkfields(o1,{'lambda','field'});
        [o1lambda,Npolo1] = expandlambda(o1);
        if Npol ~= Npolo1, error('cannot multiplex a scalar field with a dual-polarization field.');end
        E.field = [E.field,o1.field];
        E.lambda = [Elambda(:)',o1lambda(:)'];
        if nargin > 2 % optional inputs
            if ~isstruct(o2)
                error('All input parameters must be struct variables.');
            end
            options = o2;
        end
    else
        options = o1;
    end
else
    E.lambda = Elambda(:)';
end

% Is separate fields in wavelength?
issep = ~isempty(options) && isfield(options,'mux') && strcmpi(options.mux,'sepfields');

% multiplex together the wavelengths
if issep
    Eo = E; % separate fields: do not mux wavelengths
else
    Eo = multiplexcarriers(E,options);
end
% If only one wavelength
if all(Eo.lambda == Eo.lambda(1)), Eo.lambda = Eo.lambda(1);end

%--------------------------------------------------------------------------
function E = multiplexcarriers(E,options)

global GSTATE
global CONSTANTS

Nfft = GSTATE.NSAMP; % number of samples

% % Search for (exactly) duplicated lambda
% if length(unique(E.lambda)) ~= length(E.lambda)
%     warning('optilux:multiplexer','Duplicated wavelengths.');
% end

if size(E.field,1) ~= Nfft
    error('wrong number of samples');
end
if size(E.field,2) == 2*length(E.lambda)% dual polarization
    Npol = 2;
else
    Npol = 1;
end

lambda = E.lambda;

% if length(lambda) ~= length(unique(lambda))
%     error('Duplicated wavelengths');
% end

% Sort wavelengths in ascending order
lambda = sort(lambda);

if ~isempty(options)
    %%% Add a delay
    if isfield(options,'delay')
        if strcmp(options.delay,'rand')
            tau = floor((Nfft-1)*rand(1,size(E.field,2)));
        elseif all(floor(options.delay) == options.delay)
            if length(options.delay) ~= size(E.field,2)
                error('the number of elements of the delay must match the total number of polarizations');
            end
        else
            error('the delay must be made of integers');
        end
        
        for k=1:size(E.field,2)
            E.field(:,k) = fastshift(E.field(:,k),tau(k));
        end
    end
end

%%% Multiplex

freq = CONSTANTS.CLIGHT./lambda; % carrier frequencies [GHz]
maxf = max(freq); % [GHz]
minf = min(freq); % [GHz]
if (maxf - minf) > GSTATE.FSAMPLING
    reply = input(['The sampling frequency set in INIGSTATE is too small',...
        ' to avoid aliasing. Do you want to continue? Y/N [enter: no]: '], 's');
    if ~strcmpi(reply(1),'Y')
        error('Too small sampling frequency');
    end
end
if ~isempty(options) && isfield(options,'lambda')
    freqc = CONSTANTS.CLIGHT/options.lambda;
else
    freqc = (maxf+minf)/2; % central frequency [GHz] (corresponding to the zero
    %   frequency of the lowpass equivalent signal by convention)
end
deltafn = freqc - freq;   % carrier frequency spacing [GHz]
minfreq = GSTATE.FN(2)-GSTATE.FN(1);    % resolution [GHz]
ndfn = round(deltafn./minfreq);  % spacing in points
E.lambda = CONSTANTS.CLIGHT./freqc; % new central wavelength [nm]

E.field = fft(E.field);
Eout = zeros(Nfft,Npol);
for kch=1:length(lambda)    % multiplexing
    if Npol == 1
        Eout = Eout + fastshift(E.field(:,kch),-ndfn(kch));
    else
        Eout(:,1) = Eout(:,1) + fastshift(E.field(:,2*kch-1),-ndfn(kch));
        Eout(:,2) = Eout(:,2) + fastshift(E.field(:,2*kch),-ndfn(kch));
    end
end
E.field = ifft(Eout); % overwrite

%--------------------------------------------------------------------------
function [lambda,Npol] = expandlambda(A)
% Expand all lambda, and return number of polarizations

if (size(A.field,2) == 2*length(A.lambda))
    Npol = 2; % dual-polarization
else
    Npol = 1;
end

lambda = A.lambda;
if Npol*length(lambda) ~= size(A.field,2)
    error('The number of columns of the field does not match the number of wavelengths')
end

