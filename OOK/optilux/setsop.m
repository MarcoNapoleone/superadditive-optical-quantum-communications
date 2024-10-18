function E = setsop(E,x)

%SETSOP set the state of polarization
%   E=SETSOP(E,X) sets the state of polarization of the field contained in
%   the struct E according to the method described in X.
%
%   Valid strings for X are:
%
%       'rand': Each channel of the field E.field is rotated by a random
%           unitary matrix. The matrices are statistically independent.
%
%   Author: Paolo Serena, 2021
%   		University of Parma, Italy

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


Ncol = size(E.field,2);
if nargin == 1
    error('Missing methodology')
end
if ischar(x) && strcmpi(x,'rand')
    if Ncol == 2*length(E.lambda)% dual polarization
        for k=1:Ncol/2
            E.field(:,[2*k-1,2*k]) = E.field(:,[2*k-1,2*k])*randU(2);
        end
    end
else
    error('Methodology unknown.')
end