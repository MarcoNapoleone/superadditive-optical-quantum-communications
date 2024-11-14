function bin = graytobin(gray)

%GRAYTOBIN Gray decoder
%   BIN = GRAYTOBIN(GRAY) generates a Gray decoded output with the same
%   size as its input GRAY. BIN and GRAY are binary vector or matrices.
%
%   Example:
%
%      0     0     0                 0     0     0
%      0     0     1                 0     0     1
%      0     1     1    GRAYTOBIN    0     1     0
%      0     1     0  -------------> 0     1     1
%      1     1     0                 1     0     0
%      1     1     1                 1     0     1
%      1     0     1                 1     1     0
%      1     0     0                 1     1     1
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

if max(gray)>1
    error('vector or matrix must be binary');
end

num_bit = size(gray,2);

bin(:,1) = gray(:,1);
for ii=2:num_bit
    bin(:,ii) = xor(gray(:,ii),bin(:,ii-1));
end

