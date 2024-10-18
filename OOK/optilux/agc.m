function [ak, eta] = agc(pat,ak,modformat)

%AGC automatic gain control
%   IRIC = AGC(PAT,AK,Y) performs data-aided automatic (complex) gain
%   control of the digital signal AK having pattern PAT and modulation
%   format MODFORMAT (e.g., '16qam', see MODFORMATINFO). Hence, the 
%   function adjusts the constant gain and the phase rotation of AK to meet 
%   the requirements of the decision regions.
%
%   The data-aided operation is the result of a maximum-likelihood (ML)
%   estimation.
%
%   [AK ETA] = AGC(...) also returns on output the gain factor ETA such
%   that ETA*AK is maximally similar in ML sense to the expected sample.
%
%   See also: SAMP2PAT.
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

Npol = size(ak,2); % number of polarizations
Ns = size(pat,2)/Npol; % number of bits per polarization
eta = ones(1,Npol);
for k=1:Npol
    dataid = pat2samp(pat(:,1+(k-1)*Ns:k*Ns), modformat);
    eta(k) = sum(ak(:,k) .* conj(dataid))/sum(dataid.*conj(dataid));
    ak(:,k) = ak(:,k)/eta(k);
end
