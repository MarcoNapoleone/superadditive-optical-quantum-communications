function [cond,avg,nruns,stdber,wrn]=meanestimate(sig,x,nind)

%MEANESTIMATE Mean estimate by Monte Carlo simulation
%   [COND,AVG,NSAMP,STDMEAN]=MEANESTIMATE(SIG,X) estimates the
%   mean value of the random variable SIG by Monte Carlo (MC) simulation.
%
%   SIG can be a vector or a matrix. In the last case, all values are used.
%   Otherwise, with the call:
%
%   MEANESTIMATE(SIG,X,NIND)
%
%   only column NIND is used. This option is useful to run parallel Monte
%   Carlo estimations. In such a case, the struct X needs also the field
%   DIM (see below).
%
%   The Monte Carlo estimation proceeds until a certain accuracy is 
%   reached. Such accuracy is controlled by the struct X:
%
%       X.stop = vector [t1 t2]. t1 is the relative accuracy of MEAN, while
%           t2 the Gaussian confidence of the accuracy. Hence, if E[SIG] is
%           the unknown exact average value and MEAN the estimated one,
%           with Gaussian confidence the following probability holds:
%
%           Pr( E[sig]*(1-t1) < MEAN < E[SIG]*(1+t1) ) = t2/100
%
%           In such a case the MC simulation tests a sufficient number of
%           runs as soon as the accuracy is reached. For instance, with
%           t1=0.01 and t2=95, the relative error on the estimation is 0.01
%           with confidence 95%. Either small values of t1 or large values
%           of t2 increase the accuracy of the estimation, however with
%           longer computational times.
%           Beware that the concept of confidence works for Gaussian
%           distributed random variables.
%
%           A Gaussian confidence of 68% means that a Gaussian distributed
%           MEAN is within +/- one standard deviation with probability 68%.
%
%       X.nmin = minimum number of samples observed. This option can be
%           used with SIG of binary values. Default value: 1.
%
%   X can also have the following parameters:
%
%       X.minsamp: minimum number of samples analyzed. COND cannot be true
%           before observing X.MINSAMP samples. The default value is 1.
%           Note: to run MEANESTIMATE at least n times, set:
%
%               X.minsamp = n*size(SIG,1)*size(SIG,2).
%
%       X.maxsamp: maximum number of samples. If MEANESTIMATE counts more
%           than X.maxsamp samples, COND is set to false and WRN on output 
%           is set to true. The default value is Inf.
%
%       X.dim = total number of parallel estimations (only when NIND
%           exists).
%
%   On output NSAMP is the number of samples observed during the
%   estimation. STDMEAN is the standard deviation of the estimated MEAN.
%   COND is true during the MC cycle, and false when the convergence stop
%   criterion has been reached.
%
%   [COND,AVG,NSAMP,STDMEAN,WRN]=MEANESTIMATE(SIG,X) also returns
%   on output a warning, equal to true when the number of samples overcames
%   X.maxsamp.
%
%   See also EVALEYE.
%
%   Authors: Paolo Serena, 2021
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

%%%%%%%%%%%%%%%%%%%%%% INITIALIZATION

fnames=fieldnames(x);

if ~any(strcmp(fnames,'minsamp')), x.minsamp = 1; end
if ~any(strcmp(fnames,'maxsamp')), x.maxsamp = Inf; end
if ~any(strcmp(fnames,'nmin')), x.nmin = 1; end
if ~any(strcmp(fnames,'dim'))
    if nargin < 4
        x.dim = 1;
        nind = 1;
    else
        error('missing variable nind in ber_estimate');
    end
end
if any(strcmp(fnames,'stop'))
    if (x.stop(2) > 100) || (x.stop(2) < 0)
        error('the Gaussian confidence must be < 100 and > 0');
    end
    isepsilon=true;
else
    isepsilon=false;
end

[cond,avg,nruns,stdber,wrn]=mc_run(sig,x,isepsilon,nind);


%--------------------------------------------------------------------------
function [cond2,avg2,n2,stdber,wrn2]=mc_run(sig,x,isepsilon,nind)

%%%%%%%%%%%%%%%%%%%%%% MC SIMULATION
% The simulation collects M samples each iteration, hence it is a block
% iteration.

persistent n avg varavg first epsilon cond wrn; % keep memory

M = size(sig,1)*size(sig,2); % block dimension
if isempty(first)
    n = zeros(1,x.dim);  % first time only: initialization
    avg = zeros(1,x.dim);
    varavg = zeros(1,x.dim);
    first = 1; % remember for next time (persistent variable)
    cond = true(1,x.dim);  % vector of trues
    wrn = false(1,x.dim);
    if isepsilon
        epsilon = sqrt(2)*erfcinv(1-x.stop(2)/100);
    end
    if x.minsamp > x.maxsamp
        error('MIN number of samples > MAX number of samples');
    end
end
n(nind) = n(nind) + 1;
nnew = n(nind)*M;

err = sum(sum(sig));

N = (n-1)*M; % running number of samples
varerr = (err - err^2/M)/(M-1); % variance of the current block of M samples
avgerr = err/M; % mean of the current block
% varavg(nind) is the cumulative variance
varavg(nind) = ((N(nind)-1)*varavg(nind) + (M-1)*varerr + ...
    N(nind)*M/(N(nind)+M)*(avg(nind)-avgerr)^2)/(N(nind)+M-1);
avg(nind) = ((n(nind)-1)*avg(nind) + avgerr)/n(nind); % cumulative mean
stdber = sqrt(varavg./(N+M)); % standard deviation.
if isepsilon
    if nnew  >= x.minsamp
        % Note: the N+M is the MC decreasing variance factor
        if (epsilon*stdber(nind) < x.stop(1)*avg(nind)) && (avg(nind)*nnew >= x.nmin)
            cond(nind) = false; % stop MC
            if ~any(cond), first = []; end
            wrn(nind) = false;
        elseif nnew >= x.maxsamp
            cond(nind) = false; % stop MC
            if ~any(cond), first = []; end
            wrn(nind) = true;
            warning('optilux:meanestimate',['Too many samples.',...
                ' Early exit with low accuracy.']);
        end
    end
else
    if nnew  >= x.minsamp
        if avg(nind)*nnew >= x.nmin
            cond(nind) = false;
            if ~any(cond), first = []; end
            wrn(nind) = false;
        elseif nnew >= x.maxsamp
            cond(nind) = false; % stop MC
            if ~any(cond), first = []; end
            wrn(nind) = true;
            warning('optilux:meanestimate',['Too many samples.',...
                ' Early exit with low accuracy.']);
        end
    end
end

n2 = n*M;
avg2 = avg;
cond2 = cond;
wrn2 = wrn;
