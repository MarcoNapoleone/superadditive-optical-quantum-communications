function plotfield(E,flag,o1,o2)

%PLOTFIELD plot the electric field
%   PLOTFIELD(E,FLAG) plots the electric field contained in the struct
%   variable E in the current figure. E must have the fields lambda, field. 
%   See LASERSOURCE for more information.
%
%   FLAG determines the type of plot. It is a 4 char string such that:
%
%       FLAG(1) = 'p': plot the power in the time domain. '-': do not plot.
%       FLAG(2) = 'a': plot the phase ('a': angle) in the time domain. 
%           '-': do not plot.
%       FLAG(3) = 'p': plot the power spectral density (PSD) in the
%           frequency domain. '-': do not plot.
%       FLAG(4) = 'a': plot the phase of the Fourier transform in the
%           frequency domain. '-': do not plot.
%
%   For instance, FLAG='p---' plots only the power in the time domain.
%
%   PLOTFIELD(E,FLAG,X) accepts the optional parameter X of fields:
%
%       X.norm: normalization frequency. The frequency axis is normalized
%           to X.norm and the time axis to 1/X.norm. For instance, if
%           X.norm is equal to the symbol rate, the time axis is 1 after 
%           one symbol time.
%
%           Note: By default, the time axis is in [ns] and the frequency
%               axis in [GHz].
%
%       X.power: if 'xy', both the x and y polarizations are shown in
%           separate lines. By default, the total power is plotted.
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

global GSTATE

%%% Init
ts = 1/GSTATE.FSAMPLING; % sampling time [ns]
tim = 0:ts:(GSTATE.NSAMP-1)*ts; % time axis [ns]

if size(E.field,2) == 2*length(E.lambda)% dual polarization
    isy = true;
else
    isy = false;
end
if ~isy
    E.fieldx = E.field;
else
    E.fieldx = E.field(:,1:2:end);
    E.fieldy = E.field(:,2:2:end);
end
E = rmfield(E,'field');

if nargin == 3
    if isstruct(o1)
        x = o1; % options
        col = '-'; % color
    else
        x = [];
        col = o1; % color
    end
elseif nargin == 4
    col = o1; % color
    x = o2; % options
else
    x = [];
    col = '-';
end
imm = strfind(col,'--');
imd = strfind(col,'-.');
idm = strfind(col,'.-');
idd = strfind(col,'..');
if ~isempty(imm) || ~isempty(imd) || ~isempty(idm) || ~isempty(idd)
    col1 = [col(isstrprop(col,'alpha')),'-'];
else
    col1 = [col(isstrprop(col,'alpha')),'--'];
end

if ~isempty(x) && isfield(x,'norm')
    tim = tim*x.norm;
    freq = fftshift(GSTATE.FN/x.norm);
    fs = GSTATE.FSAMPLING/x.norm; % normalized sampling frequency
    timlab = 'time [symbols]';
    freqlab = 'frequency [1/symbols]';
    psdlab = 'PSD [dBm*symbols]';
else
    freq = fftshift(GSTATE.FN);
    fs = GSTATE.FSAMPLING; % sampling frequency [GHz]
    timlab = 'time [ns]';
    freqlab = 'frequency [GHz]';
    psdlab = 'PSD [dBm/GHz]';
end

%%% Check
Lnot = length(findstr(flag,'-'));
switch Lnot
    case {0,1}
        nrow = 2; ncol = 2;
    case 2
        nrow = 2; ncol = 1;
    case 3
        nrow = 1; ncol = 1;
    case 4
        return
end
for k=[1 3]
    if ~(strcmpi(flag(k),'p')  || strcmpi(flag(k),'-'))
        error('wrong flag. Valid options are, e.g., ''papa'',''p-p-'',etc');
    end
end
for k=[2 4]
    if ~(strcmpi(flag(k),'a')  || strcmpi(flag(k),'-'))
        error('wrong flag. Valid options are, e.g., ''papa'',''p-p-'',etc');
    end
end
if ~isempty(x) && isfield(x,'power') && strcmpi(x.power,'xy')
    istot = false;
else
    if ~isy, E.fieldy = zeros(size(E.fieldx)); end
    istot = true;
end

%%% Plot

npl = 1;
if strcmpi(flag(1),'p') %%% power, time-domain
    axt(npl)=subplot(nrow,ncol,npl); npl = npl + 1;
    if istot
        plot(tim,abs(E.fieldx).^2+abs(E.fieldy).^2,col);
    else
        if isy
            plot(tim,abs(E.fieldx).^2,col,tim,abs(E.fieldy).^2,col1);
            legend('x','y');
        else
            plot(tim,abs(E.fieldx).^2,col);
        end
    end
    grid on, hold on
    xlabel(timlab)
    ylabel('power [mW]')
end

if strcmpi(flag(2),'a') %%% phase, time-domain
    axt(npl)=subplot(nrow,ncol,npl); npl = npl + 1;
    if isy
        plot(tim,angle(E.fieldx),col,tim,angle(E.fieldy),col1);
        legend('x','y');
    else
        plot(tim,angle(E.fieldx),col);
    end
    grid on, hold on
    xlabel(timlab)
    ylabel('phase [rad]')
    if length(axt) > 1, linkaxes(axt,'x'); end
end

if strcmpi(flag(3),'p') || strcmpi(flag(4),'a')
    nf = 1;
    Efx = fftshift(fft(E.fieldx));
    if isy, Efy = fftshift(fft(E.fieldy)); end
    
    if strcmpi(flag(3),'p') %%% power, frequency-domain
        Nsamp = length(E.fieldx);
        axf(nf)=subplot(nrow,ncol,npl); npl = npl + 1; nf = nf + 1;
        if istot
            if isy
                plot(freq,10*log10((abs(Efx).^2+abs(Efy).^2)/Nsamp/fs),col);
            else
                plot(freq,10*log10(abs(Efx).^2/Nsamp/fs),col);
            end
        else
            if isy
                plot(freq,10*log10(abs(Efx).^2/Nsamp/fs),col,...
                    freq,10*log10(abs(Efy).^2/Nsamp/fs),col1);
                legend('x','y');
            else
                plot(freq,10*log10(abs(Efx).^2/Nsamp/fs),col);
            end
        end
        grid on, hold on
        xlabel(freqlab)
        ylabel(psdlab)
        %         Note: For instance, let isy=false. Then:
        %         mean(abs(E.fieldx).^2
        %         % coincides with:
        %         sum(abs(Efx).^2/Nsamp/fs)*(freq(2)-freq(1))
    end
    
    if strcmpi(flag(4),'a') %%% phase, frequency-domain
        axf(nf)=subplot(nrow,ncol,npl); npl = npl + 1;
        if isy
            plot(freq,angle(Efx),col,freq,angle(Efy),col1);
            legend('x','y');
        else
            plot(freq,angle(Efx),col);
        end
        grid on, hold on
        xlabel(freqlab)
        ylabel('Spectrum phase [rad]')
    end
    if length(axf) > 1, linkaxes(axf,'x'); end
end

