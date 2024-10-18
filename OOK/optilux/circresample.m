function yout = circresample(x, p, q)

% CIRCRESAMPLE circular resample
%   CIRCRESAMPLE(X,P,Q) works as RESAMPLE(X,P,Q) but with a circular, i.e.,
%   periodic, sequence X.
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

N=10; % default Neighbor term number of resample.m
[p,q] = rat( p/q, 1e-12 );
pqmax = max(p,q);
L = 2*N*pqmax + 1;
Lx = length(x);
Ltot = ceil(Lx*p/q);

if Lx >= L
    xext=[x(end-L+2:end,:);x;x(1:L-1,:)]; % extended signal
    y=resample(xext,p,q);
    L2 = ceil((L-1)*p/q);
    yout = y(L2+1:Ltot+L2,:);
else % the resample filter is longer than x
    xext=[];
    n=0;
    while length(xext) < L
        xext=[xext;x]; % cat until xext is longer than the resample filter
        n=n+1;
    end
    xext=[xext;xext]; % two-sided: replication on left and right.
    y=resample(xext,p,q);
    L2 = ceil(((n-1)*Lx+1)*p/q);
    yout = y(L2+1:Ltot+L2,:);
end