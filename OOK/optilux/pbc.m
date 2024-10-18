function E = pbc(Ex,Ey)

%PBS polarization beam combiner
%   E=PBC(EX,EY) is a polarization beam combiner. EX, EY, and E are struct
%   describing electric fields (see PBS). On output, the electric field in
%   E is the superposition of the x-polarization of EX and the
%   y-polarization of EY, alternated on columns.
%
%   Note: the combination of a beam splitter and a beam combiner is not the
%   identity operation since it rotates the reference system by 45 degrees.
%
%   See also: LASERSOURCE, PBS.
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

if Ex.lambda ~= Ey.lambda
    error('different wavelengths: use a multiplexer')
end
if size(Ex.field,2) == 2*length(Ex.lambda)% dual polarization
    Npol = 2;
else
    Npol = 1;
end
E.lambda = Ex.lambda;
E.field = zeros(size(Ex.field));
E.field(:,1:Npol:end) = Ex.field(:,1:Npol:end);
if Npol == 2, E.field(:,2:Npol:end) = Ey.field(:,2:Npol:end); end
