function dec = bintodec(bin)
%BINTODEC Convert binary patterns to decimal integers.
%   DEC = BINTODEC(BIN) converts the binary matrix BIN in the equaivalent
%   decimal representation in DEC. The function operates by rows.
%
%
%   Example
%       bintodec([0 1 0 1 1 1;1 1 1 1 1 1]) returns [23; 63]
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

nc = size(bin,2);
pow = (2.^((nc-1):-1:0))';
dec = bin*pow;

