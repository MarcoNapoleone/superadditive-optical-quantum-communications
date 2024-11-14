function gray = bintogray(bin)
%BINTOGRAY Gray encoder
%   GRAY = BINTOGRAY(BIN) generates a Gray encoded output with the same
%   size as its input parameter BIN. BIN is a binary matrix.
%
%   Example:
%
%      0     0     0                 0     0     0
%      0     0     1                 0     0     1
%      0     1     0    BINTOGRAY    0     1     1
%      0     1     1  -------------> 0     1     0
%      1     0     0                 1     1     0
%      1     0     1                 1     1     1
%      1     1     0                 1     0     1
%      1     1     1                 1     0     0
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

if any(any(~(bin == 0 | bin == 1)))
    error('The input must be binary');
end

temp =  [zeros(size(bin,1),1),bin(:,1:end-1)];
gray = xor(bin,temp);
