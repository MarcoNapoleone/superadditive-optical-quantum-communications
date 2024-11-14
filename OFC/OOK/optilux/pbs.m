function [Ex,Ey] = pbs(E)

%PBS polarization beam splitter
%   [EX,EY] = PBS(E) is a polarization beam splitter. E contains
%   information about the electric field. It is a a struct of fields:
%
%       E.lambda:   carrier wavelength [nm]
%       E.field:    time samples of the x-polarization
%
%   The block returns on output two electric fields EX and EY with 
%   orthogonal polarizations at 45 degrees with respect to E.
%   For instance, if E is fully polarized along x:
%
%   EX.field(:,1) = E.field(:,1)/sqrt(2);
%   EX.field(:,2) = 0;
%   EY.field(:,1) = 0;
%   EY.field(:,2) = E.field(:,1)/sqrt(2);
%
%   Hence, an input field polarized along the x axis is equally split
%   between the two output polarizations. The same occurs for a
%   polarization along y, but a phase shift of pi is added on EX because of
%   the unitarity of the 45 degree transformation.
%
%   See also: LASERSOURCE, PBC.
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

if size(E.field,2) ~= 2*length(E.lambda)
    error('A beam splitter can be used only in dual-polarization mode.')
end
Ex = E;
Ex.field(:,1:2:end) = (E.field(:,1:2:end) - E.field(:,2:2:end))/sqrt(2);
Ex.field(:,2:2:end) = 0;

Ey = E;
Ey.field(:,1:2:end) = 0;
Ey.field(:,2:2:end) = (E.field(:,1:2:end) + E.field(:,2:2:end))/sqrt(2);


