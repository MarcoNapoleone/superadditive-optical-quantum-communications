function E = sep2uni(E,k)

%SEP2UNI extract field from separate-field propagation
%   E = SEP2UNI(E,K) extract the Kth field from the matrix representing a
%   separate-field propagation. It thus extracts columns (2K-1,2K) in
%   dual-polarization, or column K in single-polarization.
%
%   K can be a vector of channel indexes.
%
%   In unique field propagation the function clones the input field on
%   output.
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

if length(E.lambda) == 1 % unique field
    warning('optilux:sep2uni','Nothing to do: already unique field');
else% separate fields
    if size(E.field,2) == 2*length(E.lambda) % dual polarization
        ind = [2*k-1;2*k]; % x and y polarizations are contiguous
        E.field = E.field(:,ind(:));
    else % scalar propagation
        E.field = E.field(:,k);
    end
    E.lambda = E.lambda(k);
end