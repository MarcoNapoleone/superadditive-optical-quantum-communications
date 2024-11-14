function q=ber2q(ber)
%Q=BER2Q(BER) bit-error rate to Q-factor
%	Q=BER2Q(BER) Converts the bit-error rate in Q-factor [dB]
%	using the formula:
%
%	Q = 20*log10( sqrt(2) * erfcinv(2*ber));
%
%	As a reference (BER -> Q [dB}): 
%						1e-3 -> 9.80		8e-4 ~= 10
%						1e-5 -> 12.60		2e-4 ~= 11
%						1e-9 -> 15.56		
%
%   See Also: Q2BER.
%
%   Author: Massimiliano Salsi, 2009
%   Comments: Armando Vannucci, 2009
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

q=20*log10( sqrt(2) * erfcinv(2*ber));

