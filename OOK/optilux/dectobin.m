function bin=dectobin(dec,nbit)

%DECTOBIN Convert decimal integers to binary patterns.
%   BIN = DECTOBIN(DEC,NBIT) converts the array of decimal numbers DEC in 
%   the equivalent binary pattern BIN of NBIT bits.
%
%   DEC must be a column vector.
%
%   Example
%       dectobin([1; 7],3) returns
%
%       BIN = [0 0 1; 1 1 1]
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

bin = logical(rem(floor(dec*2.^((1-nbit):0)),2));
% % substituting the previous row with the following let dectobin work
% % for a matrix
% bin = rem(floor(dec(:)*2.^((1-nbit):0)),2);
% [m,ncol]=size(dec);
% n = ncol*nbit;
% q = nbit;
% % see P. J. Acklam, "Matlab array manipulation tips and tricks", p.17-18
% % URL: http://home.online.no/~pjacklam
% X = reshape(bin,[m n/q q]);
% X = permute(X,[1 3 2]);
% bin = reshape(X,[m n]);
