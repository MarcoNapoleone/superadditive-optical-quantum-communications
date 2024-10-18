function [elec,normf] = digitalmod(pat, modformat, symbrate, ptype, par)

%DIGITALMOD linearly modulated digital signal
%   ELEC = DIGITALMOD(PAT, MODFORMAT, SYMBRATE, PTYPE, PAR)
%   creates the linearly modulated digital signal ELEC by:
%
%       ELEC = sum( ak*p( t -k*T) )
%
%   where ak are the digital symbols according to the modulation format
%   MODFORMAT and the transmitted pattern PAT.
%
%   MODFORMAT is a string (see MODFORMATINFO). Some examples (in square
%   brackets an example of the corresponding optical modulator) :
%
%   'ook': on-off keying [Mach Zehnder modulator]
%   'bpsk': binary phase-shift keying (PSK) [Mach Zehnder modulator]
%   'qpsk': quadrature phase-shift keying (QPSK) [iq modulator]
%   'psk': PSK [iq modulator + (possibly) phase modulator(s)]
%   '[num]qam': quadrature amplitude modulation. [num] is the alphabet size
%       [iq modulator]
%   '[num]pam': M-ary (power of two) pulse amplitude modulation. [num] is
%       the alphabet size. [Mach Zehnder modulator]
%
%   PAT is a matrix containing the pattern (see PATTERN). PAT can be a
%   matrix of bits, of size number of bits-by-log2(alphabet size), or a
%   column vector containing the decimal representation.
%
%   SYMBRATE is the signal symbol-rate in [Gbaud].
%
%   Note: if possible, use a sampling rate (see INIGSTATE) that is a
%       multiple of the symbol rate. If not, the function approximate their
%       ratio by a rational number, thus by first upsampling ELEC to the
%       number of samples equal to the numerator and finally, after all
%       computations, by downsampling to the denominator.
%       The upsampling part may take long times and/or occupy a large
%       amount of memory.
%
%   PTYPE is the pulse type and it is one of the following strings:
%
%   <string used in MYFILTER>:      the supporting pulse is the impulse
%                                   response of the corresponding filter in
%                                   MYFILTER. See also PAR.bw and PAR.par
%                                   and the note below for Gaussian pulses.
%
%   'rc':       raised cosine pulse. The elementary pulse takes the
%               expression (roll is the roll-off) [1]:
%                _
%               |    1       0 <= abs(t) <= (1-roll)/2*duty
%               |
%       p(t) = .     1/2*{1+cos[pi/roll/duty*...
%               |        (abs(t)-(1-roll)*duty/2)]},
%               !         if  (1-roll)*duty/2 <= abs(t) <= (1+roll)*duty/2
%               |
%                _   0      abs(t) > (1+roll)/2*duty
%
%   'rootrc':   root raised-cosine pulse. The Fourier transform of such a
%       pulse is sqrt(P(f)), with P(f) the Fourier transform of a
%       raised-cosine p(t).
%
%   'userfir': the shape of the pulse is given by the user as a column
%       vector (see PAR.firtaps).
%
%   'costails': a non-square pulse with cosine-shaped tails. Such a pulse
%       has temporal expression equal to P(t), P being the Fourier
%       transform of the raised cosine pulse (see option 'rc' above). The
%       pulse is something like the following:
%
%                       ________
%                      /        \
%                     /          \
%                     |          |
%                    /            \
%                   /              \
%       --------------------------------------------> time
%
%       This pulse is similar to real pulses not finely tuned by digital
%       signal processing thus showing some rise and fall time.
%
%   'dirac': the pulse is a Dirac delta, i.e., a Kronecker delta in the 
%       digital domain.
%
%   'rect': rectangular pulse. The duty cycle can be set by PAR.duty.
%
% PAR is a struct with the following fields (for the variable GSTATE, 
%  please see INIGSTATE):
%
%       PAR.rolloff = roll-off (0<=PAR.rolloff<=1) of raised cosine and
%           root-raised cosine pulses.
%       PAR.duty = pulse duty cycle (0 < PAR.duty <= 1). The time axis is
%           scaled by PAR.duty. Default: PAR.duty = 1.
%       PAR.norm. power normalization.
%           'iid': the output signal is normalized to unit average power
%               under the assumption that the symbols are independent and
%               identically distributed symbols with uniform distribution.
%               [Default]. See also the note below.
%           'mean': the power is normalized to the average power of the
%               current signal.
%           'no': no normalization is applied.
%
%           Note: 'iid' normalizes to the average power of an infinitely
%               long signal. hence, mean(abs(ELEC).^2) may differ from 1
%               because ELEC is of finite length. This is the correct
%               normalization, hence 'iid' is recommended compared to
%               'mean'.
%
%       PAR.nsps = number of samples per symbol. This option allows using
%           a different sampling rate for internal operations before the
%           final setting to match the global GSTATE.FSAMPLING.
%           The transmitted sequence is upsampled to PAR.Nsps samples, and 
%           then filtered with a FIR filter with impulse response equal to 
%           the desired pulse shape.
%           Default value: GSTATE.FSAMPLING/SYMBRATE. If not an integer, a
%           rational approximation with RAT is used.
%       PAR.firtaps = taps of the FIR filter creating the pulse. The
%           corresponding sampling rate is set by PAR.nsps.
%           Default: PAR.firtaps = length(PAT). This field is used only
%           with PTYPE='userfir'.
%       PAR.emph = emphasis type. After creation, the signal may experience
%           an emphasis. Valid options are:
%           PAR.emph = 'asin': an arcsin law is applied to the signal:
%
%               ELEC = asin( real(ELEC) ) + 1i*asin( imag(ELEC) )
%
%           This is useful to perfectly compensate for the sin(.) transfer
%           function of a Mach-Zehnder modulator.
%       PAR.bw = bandwidth, normalized to the symbol rate SYMBRATE, of the 
%           filter used by MYFILTER when PTYPE is a valid string supported
%           by MYFILTER.
%       PAR.par = optional parameters of MYFILTER when PTYPE is a valid 
%           string supported by MYFILTER.
%
%   [ELEC NORM] = DIGITALMOD(PAT, MODFORMAT, SYMBRATE, PTYPE, PAR) also
%   returns on output the normalization factor NORM that can be used by the
%   Mach-Zehnder modulator to undo the normalization performed by
%   PAR.emph='asin'. This way, the argument of the asin has abs value
%   smaller than 1, such that the Mach Zehnder modulator performs a
%   transparent ideal operation from the electrical to the optical domain.
%
%   Note: to generate Gaussian-shaped pulses with standard deviation sigma
%       [symbols] in the time domain, use:
%
%           PTYPE = 'gauss';
%           PAR.bw = sqrt(log(2))/(2*pi*sigma).
%
%       where sigma=1 means one symbol time (symbol rate: SYMBRATE).
%
%   See also INIGSTATE, PATTERN, MZMODULATOR, IQMODULATOR, LASERSOURCE.
%
%   References:
%       [1] M. Seimetz, "High-order Modulation for Optical Fiber
%           Transmission", Springer, 2009.
%       [2] M. Seimetz, "Multi-format Transmitters for Coherent Optical
%           M-PSK and M-QAM Transmission," in Proc. ICTON05, Th.B1.5, 2005.
%       [3] J. G. Proakis, Digital Communications, Mc Graw Hill, 4ed, 2001
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

global GSTATE;

%%% Init and checks

Nfft = GSTATE.NSAMP;
Ntini = GSTATE.FSAMPLING/symbrate; % wished samples per symbol
Nsymb = length(pat);
Nsymbupdw = ceil((Nfft-1)/Ntini);
if Nsymb < Nsymbupdw
    error('Too few symbols to fill the required number of samples.');
elseif Nsymb > Nsymbupdw
    Nsymb = Nsymbupdw;
    warning('optilux:digitalmod','Too many symbols. Pattern truncated to %d symbols.',...
        Nsymb);
end

[Nt,Nd] = rat(Ntini); % oversample, then downsample at the end
if isfield(par,'nsps')
    if floor(par.nsps) ~= par.nsps
        error('the number of samples per symbol must be an integer.');
    end
    
    Nsps = par.nsps; % the function works with Nsps samples per symbol till the final
    % resampling emulating the DAC
else
    Nsps = Nt;
end

if Nt/Ntini > 10
    warning('optilux:digitalmod',...
        ['resampling may take big times/memory. Consider using sampling ',...
        'or symbol rates such that \n[N,D]=rat(GSTATE.FSAMPLING/symbrate) contains ',...
        'a small integer N (now N=%d, D=%d). See also inigstate.m.'],Nt,Nd);
end

if ~isfield(par,'norm')
    par.norm = 'iid'; % default
end

if strcmpi(ptype,'rc') || strcmpi(ptype,'rootrc') || strcmpi(ptype,'userfir') || ...
        strcmpi(ptype,'dirac') || strcmpi(ptype,'costails') || strcmpi(ptype,'rect')
    % ok
else
    if ~isfield(par,'par'), par.par = 0; end
    if ~isfield(par,'bw'), error('Missing filter bandwidth for pulse shaping'); end
    try
        myfilter(ptype,1,par.bw,par.par);
    catch
        error('Unknown pulse type.');
    end
end

par.modformat = modformat;

%%% 1: convert the pattern into stars of the constellations
level = pat2samp(pat,modformat);

%%% 2: create a linearly modulated digital signal
elec = elecsrc(level,ptype,par,Nsymb,Nsps,Nd,Nfft);

%%% 3: resample if necessary

if Nt ~= Nsps
    % Apply DAC (Note: resample.m generates border effects. interp is better)
    if mod(Nt,Nsps) == 0
        elec = interp(elec,Nt/Nsps);
    else
        elec = circresample(elec,Nt,Nsps); % (circular) resample
    end
end

%%% 4: Perform pre-emphasis
if isfield(par,'emph')
    normf = max(max(abs(real(elec))),max(abs(imag(elec))));
    if strcmpi(par.emph,'asin')
        % thanks to normf, the result of each asin is for sure real, as
        % required by a M.ach-Zehnder modulator. However, to preserve the
        % energy, the Mach-Zehnder must know such a normalization factor.
        elec = asin(real(elec)/normf) + 1i*asin(imag(elec)/normf);
    end
else
    normf = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBFUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function elec = elecsrc(ak,ptype,par,Nsymb,Nsps,Nd,Nfft)

%   The idea is the following: the pattern is first upsampled to par.nsps
%   samples per symbol, and then filtered to create the PAM signal.

if strcmpi(ptype,'rc') || strcmpi(ptype,'rootrc') || strcmpi(ptype,'costails')
    if ~isfield(par,'rolloff'), error('Undefined rolloff'); end
    if (par.rolloff > 1 || par.rolloff < 0), error('The roll-off must be 0<= roll-off <= 1'); end
end

if ~isfield(par,'duty')
    par.duty = 1;
elseif par.duty <= 0 || par.duty > 1
    error('must be 0 < duty <= 1 ');
end

if ~(strcmpi(ptype,'rc') || strcmpi(ptype,'rootrc') ||  strcmpi(ptype,'userfir') ...
        || strcmpi(ptype,'dirac') || strcmpi(ptype,'rect') || strcmpi(ptype,'costails'))
    flag = true; % the pulse is the impulse response of myfilter
    if ~isfield(par,'par'), par.par = 0; end
    if ~isfield(par,'bw'), error('missing bandwidth for myfilter'); end
else
    flag = false;
end

%%% Modulate

if strcmpi(ptype,'userfir')
    if ~isfield(par,'firtaps')
        error('Missing FIR filter taps');
    else
        if length(par.firtaps) > Nsymb*Nsps, error('too many taps for FIR'), end
        if ~isreal(par.firtaps), error('FIR taps must be real.'); end
        par.firtaps = par.firtaps(:);
    end
end


aku = upsample(ak,Nsps);
Aku = fft(aku(1:Nsps*Nsymb)); % truncate if necessary
% Aku = repmat(fft(ak),Nsps,1); % same as fft(upsample(ak,Nsps), but faster

if flag == 1    % filter the signal
    ffir = par.duty*fftshift(-Nsps/2:1/Nsymb:Nsps/2-1/Nsymb);
    Hfir = myfilter(ptype,ffir(:),par.bw,par.par);
else
    elpulse = pulsedesign(ptype,Nsps,Nsymb,par);  % single pulse
    Hfir = fft(fftshift(elpulse));
    if strcmpi(ptype,'rootrc') % square-root raised cosine
        Hfir = sqrt(Hfir*Nsps);
        % *Nt: because I'm using filters normalized in peak spectrum (as if
        % symbol time was 1)
    end
end

elec = ifft( Aku .* Hfir ); % create PAM signal

if length(elec) < Nfft
    error('It is impossible to get the desired number of samples with the given pattern and sampling rate');
elseif length(elec) > Nfft
    elec = elec(ceil(1:Nd:Nfft*Nd));
end

%%% normalize to unit power
if strcmpi(par.norm,'iid') %i.i.d. symbols
    % See [3], power spectra of linearly modulated signals
    y = modformatinfo(par.modformat);
    varak = y.symb.var; % expected variance
    meanak = y.symb.mean; % expected value or mean
    avge = (varak*sum(abs(Hfir).^2)/Nsymb + ...
        abs(meanak)^2*sum(abs(Hfir(1:Nsymb:end)).^2))/Nsps^2;
elseif strcmpi(par.norm,'mean')
    avge = mean(abs(elec).^2);
elseif strcmpi(par.norm,'no')
    avge = 1;
else
    error('Unknwon normalization method');
end
elec = elec/sqrt(avge);

%--------------------------------------------------------------------------
function elpulse=pulsedesign(ptype,Nsps,Nsymb,par)

%PULSEDESIGN Creates the fundamental pulse
%   Y=PULSEDESIGN(PTYPE,NSPS,NSYMB,PAR) returns in the vector [NSYMB*NSPS,1]
%   the fundamental pulse whose type is defined in PTYPE.
%
%   PAR is a struct (see main help of DIGITALMOD).

elpulse = zeros(Nsps*Nsymb,1);

switch ptype
    case {'rc','rootrc'}
        tfir = (1/par.duty)*(-Nsymb/2:1/Nsps:Nsymb/2-1/Nsps)';
        elpulse = sinc(tfir).*cos(pi*tfir*par.rolloff)./(1-(2*par.rolloff*tfir).^2);
        elpulse(~isfinite(elpulse)) = par.rolloff/2*sin(pi/2/par.rolloff); % Avoid NaN
        
    case 'costails'
        nl = round(0.5*(1-par.rolloff)*par.duty*Nsps);    % start index of cos roll-off
        nr = par.duty*Nsps-nl-1;                       % end index of cos roll-off
        
        nmark = 1:nl;                       % indices where the pulse is 1
        ncos  = nl:nr;                      % transition region of cos roll-off
        
        elpulse(Nsps*Nsymb/2+nmark) = 1;
        hperiod = par.duty*Nsps-2*nl;
        if hperiod ~= 0
            elpulse(ncos+Nsps*Nsymb/2+1) = 0.5*(1+cos(pi/(hperiod)*(ncos-nl+0.5)));
        end
        elpulse(1:Nsps*Nsymb/2) = flipud(elpulse(Nsps*Nsymb/2+1:Nsps*Nsymb/2*2)); % first half of the pulse
        
    case 'userfir'
        Ntaps = length(par.firtaps);
        Ntapshalf = ceil(Ntaps/2);
        par.firtaps = ifftshift(par.firtaps); % zero padding will be done in the middle
        Npulse = length(elpulse);
        Npad = Npulse-Ntaps;
        elpulse = [par.firtaps(1:Ntapshalf);zeros(Npad,1);par.firtaps(Ntapshalf+1:end)];
        elpulse = fftshift(elpulse); % fftshift(ifftshift(x)) = x
        
    case 'dirac'
        elpulse(Nsps*Nsymb/2+1) = 1;
        
    case 'rect'
        nl = round(0.5*par.duty*Nsps);
        if nl == 0
            elpulse(Nsps*Nsymb/2+1) = 1;  % same as Dirac's delta.
            % Why this? Because in this way you can create a Gaussian pulse, for
            % instance, by filtering this train of delta.
        else
            elpulse(Nsps*Nsymb/2-nl+1:Nsps*Nsymb/2+nl) = 1;
        end
        
    otherwise
        error('The pulse ptype does not exist');
end
