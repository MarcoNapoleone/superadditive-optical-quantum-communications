function y=circfilter(h,x)

%CIRCFILTER Fast circular convolution using filter
%   Y=CIRCFILTER(H,X) returns in Y the circular convolution between the 
%   filter H and the signal X. This operation is faster than
%   ifft(fft(X) .* fft(H)) for length(H) << length(X).
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

M = length(h);
y = filter(h,1,x);
if size(x,1) == 1
    extx = [x(end-M+2:end),x(1:M-1)];
else
    extx = [x(end-M+2:end);x(1:M-1)];
end
yfirst = filter(h,1,extx);
y(1:M-1) = yfirst(M:end);
