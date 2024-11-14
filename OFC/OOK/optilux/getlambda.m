function lambda = getlambda(lamc,spac,Nch)

%GETLAMBDA generate wavelengths evenly spaced in the frequency domain
%   LAMBDA = GETLAMBDA(LAMC,SPAC,NCH) generates a vector of wavelengths
%   LAMBDA [nm] evenly spaced in the frequency domain. LAM [nm] is the
%   central wavelength, SPAC the spacing [nm], NCH the number of
%   wavelengths.
%
%   Note: the channel spacing is uniform in the frequency domain, according
%       to the ITU-T recommendations. Hence, the spacing is slightly
%       non-uniform in the wavelength domain. For instance, with LAM = 1550
%       , SPAC = 0.4 and NCH = 5, the following wavelengths, and the
%       corresponding carrier frequencies, are created (only four digits
%       are shown):
%
%       [1549.2004 1549.6001 1550.0000 1550.4001 1550.8004]     [nm]
%       [193.3147  193.3646  193.4145  193.4644  193.5143]      [THz]
%
%   See Also: LASERSOURCE
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

global CONSTANTS

freq = CONSTANTS.CLIGHT./lamc; % [GHz]
DF = spac./lamc.^2*CONSTANTS.CLIGHT; % frequency spacing [GHz]
lambda = zeros(1,Nch);
for chan = 1:Nch
    freqt = freq + DF*(chan-(Nch+1)/2);
    lambda(Nch-chan+1) = CONSTANTS.CLIGHT./freqt;
end
