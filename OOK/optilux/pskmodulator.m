function E = pskmodulator(E,elec,par)
%PSKMODULATOR optical phase-shift keying (PSK) modulator
%   PSKMODULATOR(E,ELEC) modulates the laser source E with the electrical
%   signal ELEC returned by DIGITALMOD.
%
%   For PSK of order higher than 4, the psk modulation can be done by a IQ
%   modulator (default) or by the cascade of a QPSK modulation and PSK
%   modulators [1]. Such an option can be set by PAR:
%
%   PAR.mode = 'iq. IQ modulation (default). 'parallel-psk': cascade
%       QPSK+PSKs.
%
%   See also DIGITALMOD, MZMODULATOR, IQMODULATOR, LASERSOURCE
%
%   References:
%       [1] M. Seimetz, "High-order Modulation for Optical Fiber
%           Transmission", Springer, 2009, p.28.
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

ncol = size(elec,2);
if nargin == 2
    if ncol > 2
        error('IQ modulation expects just two columns on input signal.')
    end
end

if ncol == 1 && all(isreal(E.field))    % bpsk
    E = mzmodulator(E,elec,par);
else            % qpsk or psk
    E = iqmodulator(E, elec(:,1),par);
    for kk=2:ncol % m-psk with m>4
        
        % '0': do not rotate
        % '1': rotate by pi/(2^(kk+1))
        E = pm_modulator(E,pi/2^kk,elec(:,kk));
    end
end


%--------------------------------------------------------------------------
function E=pm_modulator( E, phi, modsig)

modsig = (modsig+1)/2; % remove average of normalized BPSK signal

if size(E.field,2) == 2*length(E.lambda)% dual polarization
    Npol = 2;
else
    Npol = 1;
end

for m=1:Npol
    if any(any(E.field(:,m)))
        E.field = E.field.*fastexp(phi*modsig);
    end
end