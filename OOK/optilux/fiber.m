function [E,zbrf]=fiber(E,x)

%FIBER Single-mode optical fiber in the nonlinear regime
%   FIBER(E,X) solves the nonlinear Schroedinger equation (NLSE). E is the
%   propagating electric field, such that E.fieldx are the time-samples of
%   the x-polarization,  E.fieldy the samples of the y-polarization (if
%   existing), and E.lambda is the carrier wavelength.
%
%   X contains fiber parameters. It is a structure of fields:
%
%       X.length:     fiber length [m]
%       X.alphadB:    fiber attenuation [dB/km]
%       X.aeff:       fiber effective area [um^2] of the fundamental mode
%       X.n2:         fiber nonlinear index [m^2/W]
%       X.lambda:     wavelength [nm] of fiber parameters
%       X.disp:       fiber chromatic dispersion coefficient [ps/nm/km]
%                     @ X.lambda.
%       X.slope:      fiber slope [ps/nm^2/km], i.e. derivative of X.disp
%                     @ X.lambda.
%
%                     Note: If you want the beta3 coefficient = 0 set:
%                       X.slope = -2*X.disp/X.lambda.
%
%   The propagation along the fiber is simulated by the split-step Fourier
%   method (SSFM). Different choices are possible for the step setup:
%
%       X.steptype:   step type. It can be 'asymm' (asymmetric step) or
%                     'symm' (symmetric step) [Sha14]. Default: 'asymm'.
%       X.dzmax:      max. step [m]. Default: x.length.
%       X.dphi1fwm:   true: set the first step according to the maximum
%                     four-wave mixing (FWM) phase criterion [Mus18].
%                     false: set the first step according to the maximum
%                     nonlinear phase criterion [Sin03,Mus18].
%                     Default: true.
%
%                     Note: Setting the first step according to the maximum
%                           FWM is a wise choice since it shrinks the step
%                       	for increasing bandwidth and dispersion.
%       X.stepupd:    Step update rule. Sets the updating rule for all
%                     steps following the first one. It can be 'nlp'
%                     (nonlinear phase criterion, [Sin03]) or 'cle'
%                     (constant local error criterion, [Zha05,Zha08]).
%                     Default: 'cle'.
%       X.dphimax:    Accuracy parameter. If X.dphi1fwm is true, X.dphimax
%                     is the maximum FWM phase [rad] in the first step.
%                     In such a case, the reference bandwidth may be
%                     provided, see X.bandwidth.
%                     If X.dphi1fwm=false, X.stepupd is forced to 'nlp'
%                     and X.dphimax is the maximum nonlinear phase [rad]
%                     per step [Sin03].
%                     If X.dphimax does not exist, the following values
%                     [rad] are used [Mus18] over the signal bandwidth:
%
%                        X.dphimax = 20  if X.dphi1fwm & X.stepupd='cle'
%                        X.dphimax = 4   if X.dphi1fwm & X.stepupd='nlp'
%                        X.dphimax = 0.003 if ~X.dphi1fwm & X.stepupd='nlp'
%
%                     Note: if X.stepud='nlp' and X.dphi1fwm is true
%                       do not confuse X.dphimax with the maximum nonlinear
%                       phase rotation per step. Such a value is calculated
%                       internally to FIBER to match the first step length
%                       given by X.dphimax.
%       X.bandwidth:  Bandwidth [GHz] for the step size set-up. Usually the
%                     bandwidth of the propagating wavelength-division
%                     multiplexing (WDM) signal. Default: the simulation
%                     bandwidth, GSTATE.FSAMPLING (See INIGSTATE).
%       X.trace:      true: print SSFM information on the screen. Default:
%                     false.
%
%   The asymmetric/symmetric step discretize the fiber in the following way
%
%      Asymmetric                          Symmetric
%   -----------------       ------------------------------------------
%   | N | L | N | L |...    | L/2 | N | L/2 | L/2 | N | L/2| N | L/2 | ...
%   --------*--------       ----------*-------------------------------
%
%   where N is the nonlinear block (including attenuation, [Sha14]), while
%   L is the linear block. The next step choice is made at the end of the
%   step in the asymmetric case or at the end of the nonlinear block in the
%   symmetric case (see * in the above figure). In the symmetric case, the
%   concatenation L/2 + L/2 is done in a single linear propagation.
%
%   In dual-polarization, the fiber is described by the concatenation of
%   waveplates, such that within each waveplate the fiber has constant
%   parameters.
%   In such a case X has the following additional fields:
%
%       X.nplates:   Number of waveplates, also called trunks. In each
%                    waveplate, the propagation is modeled as for a
%                    polarization-maintaining fiber [Hobook]. Default: 100.
%       X.coupling:  Polarization coupling mode. It can be 'none' for no
%                    coupling, or 'pol' for strong polarization coupling.
%                    The default value depends on the PMD parameter and the
%                    Kerr effect.
%       X.pmdpar:    Polarization-mode dispersion (PMD) parameter
%                    [ps/sqrt(km)] in strong-coupling mode. The 
%                    differential group delay (DGD) is a random variable 
%                    whose average value accumulated along X.length [m] is
%                    X.pmdpar*sqrt(X.length*1e-3) [ps].
%                    The root mean square (rms) value is related to the
%                    average value by:
%
%                       DGDrms = DGDmean * sqrt(3*pi/8)
%
%                    Note: within a waveplate, the DGD grows linear with
%                    distance. The square root law [Hobook] is created by 
%                    the concatenation of many waveplates.
%		X.beatlength:beat length [m] between polarizations. It is related
%                    to the average differential propagation constant Dbeta
%                    between polarizations by:
%
%                       X.beatlength = 2*pi/<Dbeta>
%
%                    Default: 20 [m].
%       X.ismanakov: true: nonlinear Kerr effect is modeled by the
%                    Manakov equation [Mar97, Ant16, Mum14]. false: solve
%                    the coupled-NLSE (CNLSE). Default: false.
%
%   [E,OUT]=FIBER(E,X) returns in OUT a struct containing some useful
%   parameters used by the fiber:
%
%   OUT.firstdz:      first step length [m] used by the SSFM.
%   OUT.ncycle:       number of SSFM iterations.
%   OUT.time:         elapsed time [s] within the function FIBER.
%
%   If X.exportpar=true the following additional parameters are available
%   on output:
%
%   OUT.lin:          birefringence coupling matrices.
%   OUT.betat:        beta(omega) used by FIBER, with omega=2*pi*GSTATE.FN,
%                     see INIGSTATE.
%   OUT.lcorr:        correlation length [m]: X.length/X.nplates.
%   OUT.npol:         number of polarizations.
%
%   See also AMPLIFLAT.
%
%   References:
%
%      [Sin03]  O. V Sinkin, R. Holzlöhner, J. Zweck, and C. R. Menyuk,
%               “Optimization of the split-step Fourier method in modeling
%               optical-fiber communications systems," J. Lightw. Technol.,
%               vol. 21, no. 1, pp. 61–68, 2003.
%      [Zha05]  Q. Zhang and M. Hayee, “An SSF scheme to achieve
%               comparable global simulation accuracy in WDM systems," IEEE
%               Photon. Technol. Lett., vol. 17, no. 9, pp. 1869–1871, 2005
%      [Zha08]  Q. Zhang and M. Hayee, “Symmetrized split-step Fourier
%               scheme to control global simulation accuracy in fiber-optic
%               communication systems," J. Lightw. Technol., vol. 26, no. 2
%               pp. 302–316, 2008.
%      [Mus18]  S. Musetti, P. Serena and A. Bononi, "On the Accuracy of
%               split-step Fourier simulations for wideband nonlinear
%               optical communications," J. Lightw. Technol., vol. 36,
%               no. 23, pp. 5669–5677, Dec 2018.
%      [Sha14]  J. Shao, X. Liang, and S. Kumar, “Comparison of Split-Step
%               Fourier Schemes for Simulating Fiber Optic Communication
%               Systems," IEEE Photonics J., vol. 6, no. 4, 2014.
%      [Df00]   A. O. Dal Forno, A. Paradisi, R. Passy, and J. P.
%               Von Der Weid, “Experimental and theoretical modeling of
%               polarization-mode dispersion in single-mode fibers," IEEE
%               Photon. Technol. Lett., vol. 12, no. 3, pp. 296–298, 2000.
%      [Mar97]  D. Marcuse, C. R. Menyuk, and P. K. A. Wai, “Application of
%               the Manakov-PMD equation to studies of signal propagation
%               in optical fibers with randomly varying birefringence,"
%               J. Lightw. Technol., vol. 15, no. 9, pp. 1735–1746, 1997.
%      [Ant16]  C. Antonelli, M. Shtaif, and A. Mecozzi, "Modeling of
%               nonlinear propagation in Space-Division Multiplexed
%               fiber-optic transmission," J. Lightw. Technol.,vol.34,no.1,
%               pp. 366-554, Jan 2016
%      [Hobook] K.-Po Ho and J. M. Kahn, "Mode Coupling and its Impact on
%               Spatially Multiplexed Systems," in Optical Fiber
%               Telecommunications VIB, Systems and Networks, I. Kaminov,
%               T. Li and A. Willner eds, Academic Press, 2013.
%      [Ack00]  P. J. Acklam, "MATLAB array manipulation tips and tricks,"
%               tech. report, University of Oslo, Norway, 5 May 2000.
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

global CONSTANTS;               % CONSTANTS is a global structure variable.
CLIGHT = CONSTANTS.CLIGHT;      % speed of lightin vacuum [m/s]
global GSTATE                   % GSTATE is a structure whose fields are
% defined in inigstate.m.

DEF_PLATES = 100;       % number of waveplates when birefringence is on
DEF_BEATLENGTH = 20;    % Default beat length [m]
DEF_STEPUPD = 'cle';    % Default step-updating rule
DEF_ISSYMM = false;     % Default step computation
DEF_PHIFWMCLE = 20;     % Default X.dphimax [rad] for CLE stepupd rule
DEF_PHIFWMNLP = 4;      % Default X.dphimax [rad] for NLP & X.dphi1fwm
DEF_PHINLP = 0.003;     % Default X.dphimax [rad] for NLP & ~X.dphi1fwm

timeini = tic;

%%% Init

checkfields(x,{'length','alphadB','n2','lambda','slope','dzmax','dphimax',...
    'disp','aeff','ismanakov','nplates','coupling','dphi1fwm','bandwidth',...
    'stepupd','trace','steptype','db0','dgd','lin',...
    'pmdpar','beatlength','exportpar'});
% check that the fields of X really exist

nfnames = fieldnames(x);

if size(E.field,2) == 2*length(E.lambda)% dual polarization
    x.isy = true;
else
    x.isy = false;
end

Nfft = size(E.field,1);

if length(E.lambda) > 1
    x.isuni = false; % separate-field propagation
else
    x.isuni = true; % unique-field propagation
end

%%% Check mandatory parameters
if ~isfield(x,'length'), error('Missing fiber length [m].'); end
if ~isfield(x,'lambda'), error('Missing wavelength [nm] of fiber parameters.'); end
if ~isfield(x,'alphadB'), error('Missing fiber attenuation [dB/km].'); end
if ~isfield(x,'disp'), error('Missing fiber dispersion [ps/nm/km].'); end
if ~isfield(x,'slope'), error('Missing fiber slope [ps/nm^2/km].'); end
% if ~isfield(x,'pmdpar'), error('Missing fiber PMD parameter [ps/sqrt(km)].'); end % later
if ~isfield(x,'n2'), error('Missing fiber nonlinear index n2 [m^2/W].'); end
if ~isfield(x,'aeff'), error('Missing fiber effective area [um^2].'); end

%%% Check if Kerr effect is active
if isinf(x.aeff) || (x.n2 == 0) % no Kerr
    x.iskerr = false;
else
    x.iskerr = true;
end

%%% Set-up coupling
if isfield(x,'coupling')
    if ~(strcmpi(x.coupling,'none') || strcmpi(x.coupling,'pol'))
        error('Coupling must be ''none'' or ''pol''.');
    end
end

if ~x.isy % scalar propagation
    
    if isfield(x,'ismanakov') && x.ismanakov
        warning('optilux:fiber','Cannot use Manakov equation in scalar propagation: forced to NLSE.')
    end
    x.ismanakov = false;   % because of no birefringence in scalar case
    if isfield(x,'coupling') && ~strcmpi(x.coupling,'none')
        warning('optilux:fiber','coupling not possible in scalar propagation');
    end
    if isfield(x,'pmdpar') && x.pmdpar ~= 0
        warning('optilux:fiber','PMD does not exist in scalar propagation: set to 0.');
    end
    x.pmdpar = 0;
    x.beatlength = 0;
    x.coupling = 'none'; % coupling impossible with one pol
    
else % dual-polarization
    
    if ~isfield(x,'pmdpar')
        error('Missing fiber PMD parameter [ps/sqrt(km)].');
    end
    if isfield(x,'ismanakov') && x.ismanakov
        if x.pmdpar == 0
            % Note: coupling indeed makes a difference for the Kerr effect of CNLSE
            x.coupling = 'none';
        else
            if isfield(x,'coupling') && strcmpi(x.coupling,'none')
                warning('optilux:fiber','No coupling but PMD ~= 0. ');
            else
                x.coupling = 'pol';
            end
        end
        
    else % CNLSE
        
        x.ismanakov = false;
        if ~isfield(x,'coupling')
            if ~x.iskerr && x.pmdpar == 0
                x.coupling = 'none';
            else
                x.coupling = 'pol';
            end
        end
    end
end

if strcmpi(x.coupling,'pol')  % X-Y coupling
    if ~isfield(x,'beatlength')
        x.beatlength = DEF_BEATLENGTH; % default value
    end
    if ~any(strcmp(nfnames,'nplates'))
        x.nplates = DEF_PLATES;    % default value
    end
    lcorr = x.length/x.nplates; % waveplate length [m] (aka correlation length)
    dgd1 = x.pmdpar/sqrt(lcorr)*sqrt(3*pi/8)/sqrt(1000)*1e-3; % Differential
    % group delay (DGD) per unit length [ns/m] @ x.lambda within a
    % waveplate. To get [Df00] remember that within a waveplate the delay is
    % dgd1*lcorr.
    x1.coupling = x.coupling;
    x1.beatlength = x.beatlength;
    if ~isfield(x,'lin')
        lin.matin = zeros(2,2,x.nplates);
        lin.db0 = zeros(x.nplates,2);
        for l=1:x.nplates % SVD, hence different with the old FIBER version,
            [lin.matin(:,:,l),lin.db0(l,:)] =  eigendec(x1);
        end
        % lin.db0 extended later
    else
        lin = x.lin;
    end
else
    dgd1 = 0; % turn off all polarization and birefringence effects
    x.nplates = 1;
    x.coupling = 'none';
    lin.matin =  1; % scalar linear effects
    lin.db0 = 0;
end
lin.isuni = x.isuni;
%%% SSFM checks

% maximum step
if ~any(strcmp(nfnames,'dzmax')) || (x.dzmax > x.length)
    x.dzmax = x.length;
end
% trace on screen
if ~any(strcmp(nfnames,'trace'))
    x.trace = false;
end
% first step criterion
if ~any(strcmp(nfnames,'dphi1fwm'))
    x.dphi1fwm = true; % FWM criterion for the first step setup.
end
% step update
if x.iskerr
    if any(strcmp(nfnames,'stepupd'))
        if ~strcmpi(x.stepupd,'cle') && ~strcmpi(x.stepupd,'nlp')
            error('Unknown step update rule')
        end
    else
        x.stepupd = DEF_STEPUPD; % default
    end
    if strcmpi(x.stepupd,'cle')
        if ~x.dphi1fwm
            error('the combination X.dphi1fwm=false and X.stepupd=''cle'' is not possible');
        end
        x.iscle = true;
    else
        x.iscle = false;
    end
    
    % first step parameter value
    if ~any(strcmp(nfnames,'dphimax'))
        if x.dphi1fwm
            if strcmpi(x.stepupd,'cle')
                x.dphimax = DEF_PHIFWMCLE; % CLE
            else
                x.dphimax = DEF_PHIFWMNLP; % NLP + PhiFWM in 1st step
            end
            if ~isfield(x,'bandwidth')
                x.dphimax = x.dphimax*(1.5)^2; % because in [Mus18] they used
                % the signal bandwidth, x1.5 smaller than the simulation
                % bandwidth.
            end
        else
            x.dphimax = DEF_PHINLP; % NLP
        end
    end
end
% step type
if any(strcmp(nfnames,'steptype'))
    if strcmpi(x.steptype,'asymm')
        x.issym = false;
    elseif strcmpi(x.steptype,'symm')
        x.issym = true;
    else
        error('Wrong step type')
    end
else
    x.issym = DEF_ISSYMM;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% SET-UP PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x.alphalin = (log(10)*1e-4)*x.alphadB; % [m^-1]

%%%%%%%%%% Linear Parameters

omega = 2*pi*GSTATE.FN(:);  % angular frequency [rad/ns]

% b00 = 0; % phase reference of propagation constant
b10 = 0; % retarded time frame
b20 = -x.lambda^2/2/pi/CLIGHT*x.disp*1e-6; % beta2 [ns^2/m] @ x.lambda
b30 = (x.lambda/2/pi/CLIGHT)^2*(2*x.lambda*x.disp+x.lambda^2*x.slope)*1e-6;
% beta3 [ns^3/m] @ x.lambda

betat = zeros(size(E.field));
% Domega_ik: [1/ns]. "i" -> at ch. i, "0" -> at x.lambda, "c" -> at lambda of frequency 0
Domega_i0 = 2*pi*CLIGHT*(1./E.lambda - 1/x.lambda); % [1/ns]
ipm = 1;
if x.isuni % unique field
    beta1 = b10 + dgd1/2;    % [ns/m] @ E.lambda
    beta2 = b20 + b30*Domega_i0;  % beta2 [ns^2/m] @ E.lambda
    beta3 = b30; % [ns^3/m] @ E.lambda
    betat(:,ipm) = omega.*(beta1 + omega.*(beta2/2 + omega*beta3/6));
    % betat: deterministic beta coefficient [1/m]
    if x.isy % Add DGD on polarizations
        betat(:,ipm+1) = betat(:,ipm) - dgd1/2*omega;
    end
else % separate-field
    freq = CLIGHT./E.lambda; % carrier frequencies [GHz]
    maxf = max(freq); % [GHz]
    minf = min(freq); % [GHz]
    freqc = (maxf+minf)/2; % central frequency [GHz] (corresponding to the zero
    %   frequency of the lowpass equivalent signal by convention)
    lamc = CLIGHT/freqc; % central wavelength [nm]

    Domega_ic = 2*pi*CLIGHT*(1./E.lambda - 1/lamc); % Domega [1/ns] channel i vs central lambda
    Domega_c0 = 2*pi*CLIGHT*(1./lamc - 1/x.lambda); % Domega [1/ns] central lambda vs lambda fiber parameters
%     beta0 = b10*Domega_ic + 0.5*b20*(Domega_i0.^2-Domega_c0.^2) + ...
%         b30/6*(Domega_i0.^3-Domega_c0.^3); % [ns]
    b1 = b10 + b20*Domega_ic+0.5*b30*(Domega_i0.^2-Domega_c0.^2);  %ch's beta1 [ns/m]
    beta1 = b1 + dgd1/2;   % [ns/m] @ GSTATE.LAMBDA
    beta2 = b20 + b30*Domega_i0;  % beta2 [ns^2/m]@ E.lambda
    
    ich = 1;
    for kch=1:length(E.lambda)
        betat(:,ich) = omega*beta1(kch)+0.5*omega.^2*beta2(kch)+...
            omega.^3*b30/6;   % beta coefficient [1/m]
        x.b1(:,ich) = omega*beta1(kch);% needed by walk-off filter
        ich = ich + 1;
        if x.isy % Add DGD on polarizations
            betat(:,ich) = betat(:,ich-1) - dgd1*omega;
            x.b1(:,ich) = x.b1(:,ich-1) - dgd1*omega;
            ich = ich + 1;
        end
    end
    if ~isfield(x,'lin') && ~x.isuni && strcmpi(x.coupling,'pol')
        % https://it.mathworks.com/matlabcentral/answers/417541-how-do-i-duplicate-a-certain-column-of-a-matrix
        lin.db0 = repmat(lin.db0,1,length(E.lambda)); % duplicate columns
        lin.db0 = bsxfun(@plus,lin.db0,kron(dgd1/2*Domega_ic,[1 1]));
    end

end

if isa(E.field,'gpuArray')
    betat = gpuArray(betat); % put in GPU even betat to speed up
end

%%%%%%%%%% Nonlinear Parameters
if ~x.iskerr
    gam = 0;
else
    gam = 2*pi*x.n2./(E.lambda*x.aeff(1))*1e18; % nonlinear coeff [1/mW/m]
    if ~isfinite(gam)
        error('Cannot continue: not finite nonlinear Kerr coefficient.');
    end
end
% x.aeff(1) is correct, because gam is normalized to fundamental
% mode [Ant16,Mum13].

if x.ismanakov % Manakov
    if x.isy
        nl_mat.coeff = 8/9;
    else
        nl_mat.coeff = 1; % Not a Manakov actually, but same model
    end
else % CNLSE
    nl_mat.coeff = 1;
end
x.gam =  nl_mat.coeff*gam; % [1/mW/m], including Manakov correction, if active

if (x.disp == 0 && x.slope == 0) || any(x.gam == 0)  % only GVD or only Kerr
    x.dphimax = Inf;
    x.dzmax = x.length;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SSFM PROPAGATION  %%%%%%%%%%%%%%%%%%%%%%%%%%%

[firstdz,ncycle,E]=ssfm(E,lin,betat,x);

if nargout
    if any(strcmp(nfnames,'exportpar')) && x.exportpar
        zbrf.lin = lin;
        zbrf.betat = betat;
        zbrf.lcorr = x.length/x.nplates;
        zbrf.npol = 1+x.isy;
    end
    zbrf.time = toc(timeini);
    zbrf.ncycle = ncycle;
    zbrf.firstdz = firstdz;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUB-FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [firstdz,ncycle,E]=ssfm(E,lin,betat,x)

%SSFM Split-step Fourier method
%   The function solves the coupled nonlinear Schroedinger equation or the
%   Manakov equation [Mar97,Ant16] through the split-step Fourier method
%   (SSFM).
%
%   [FIRSTDZ,NCYCLE,E]=SSFM(E,LIN,BETAT,X)
%
%   E contains the electric field and its information. Hence E.fieldx are
%   the time samples of the x-polarization, E.fieldy of the y-polarization
%   (if exists), and E.lambda the wavelength [nm].
%
%   BETAT is the deterministic beta(omega) coefficient.
%
%   LIN is a struct containing information about the waveplates. In
%   particular, LIN.MATIN contains the unitary matrix for change of basis,
%   while LIN.DB0 contains the differential propagation constant between
%   polarizations.
%
%   X contains parameters regarding the fiber and the SSFM implementation.

%%% Init

if x.trace, fprintf('\n%-12s %-8s %-9s\n\n','Stepupd','step #','z [m]');end

% Concatenate polarizations by columns. It works for both
% unique and separate fields case. We have:
%
%   unique field   : u = (Nfft,Npol)
%   separate fields: u = (Nfft,Npol,Nchannel)

x.chlambda = E.lambda; % save for later, because E will be canceled
if isa(E.field,'gpuArray')
    E.field = gpuArray(E.field);
end

%%% SSFM

ncycle = 1;                 % number of steps
lcorr = x.length/x.nplates; % waveplate length [m]

[dz,x.dphimax] = firststep(E.field,x);
% alldz = dz;

halfalpha = 0.5*x.alphalin; % [1/m]
firstdz = dz;
zprop = dz;             % running distance [m]
if x.issym % first half LIN step outside the cycle
    [dzb,nindex] = checkstep(dz/2,dz/2,lcorr);
    E.field = linear_step(lin,betat,dzb,nindex,E.field);% Linear step
    % 	fprintf('dz in FIRST LIN = %.2f.\n',dz/2);
else
    dzs = 0; % symmetric step contribution.
end

while zprop < x.length   % all steps except the last
    
    %%% Nonlinear step
    E.field = nonlinear_step(E.field,x,dz);
    %     fprintf('\n%d. dz in NL = %.2f. zprop = %.2f\n',ncycle,dz,zprop);
    
    %%% Linear step 1/2: attenuation (scalar)
    E.field = E.field*exp(-halfalpha*dz);
    
    if x.issym
        dzs = nextstep(E.field,x,dz);
        if zprop+dzs>x.length, dzs=x.length-zprop; end % needed in case of last step
        hlin = (dz+dzs)/2;
    else
        hlin = dz;
    end
    
    %%% Linear step 2/2: GVD + birefringence
    [dzb,nindex] = checkstep(zprop+dzs/2,hlin,lcorr); % zprop+dzs/2: end of step
    
    E.field = linear_step(lin,betat,dzb,nindex,E.field);% Linear step
    %     fprintf('dz in LIN = %.2f = (dz/2 = %.2f) + (dzs/2 = %.2f).\n',hlin,dz/2,dzs/2);
    
    if x.issym
        dtmp = dz; dz = dzs; dzs = dtmp; % exchange dz and dzs
    else
        dz = nextstep(E.field,x,dz);
    end
    
    if x.trace,fprintf('%-12s %-8d %-9.2f\n',x.stepupd,ncycle,zprop);end
    
    zprop = zprop+dz;
    ncycle = ncycle+1;
    %     alldz(ncycle) = dz;
end
last_step = x.length-zprop+dz; % last step
zprop = zprop -dz + last_step;
% alldz(ncycle) = last_step;

if x.trace,fprintf('%-12s %-8d %-9.2f\n',x.stepupd,ncycle,zprop);end

if x.issym
    hlin = last_step/2;
else
    hlin = last_step;
end

%%% Last Nonlinear step
if x.gam ~= 0, E.field = nonlinear_step(E.field,x,last_step); end % last NL
% fprintf('dz in LAST1 NL = %.2f. zprop = %.2f\n',last_step,zprop-dz+last_step);

E.field = E.field*exp(-halfalpha*last_step);

%%% Last Linear step: GVD + birefringence
[dzb,nindex] = checkstep(x.length,hlin,lcorr);

E.field = linear_step(lin,betat,dzb,nindex,E.field); % last LIN


% fprintf('dz in LAST1 LIN = %.2f.\n',hlin);
% fprintf('issym=%d. LAST coordinate = %.2f\n\n',x.issym,sumdz)

E.lambda = x.chlambda;

% brf.betat = betat;
% brf.db1 = db1;

%--------------------------------------------------------------------------
function u = linear_step(lin,betat,dzb,nindex,u)

%LINEAR_STEP linear step of the optical fiber
%   U=LINEAR_STEP(LIN,BETAT,DZB,NINDEX,U,ISUNI) applies the linear
%   vectorial step.
%
%   U is a matrix containing the electric field temporal samples, whose
%   size depends on the propagation type ISUNI:
%
%       unique field:   size(U) = [@ samples, number of pol]
%       separate-field: size(U) = [@ samples, number of pol, number channels]
%
%   BETAT is a matrix of the same size of U and contains the
%   beta(omega) coefficient, i.e.
%   BETAT= beta0 + (beta1+dgd)*omega + beta2/2*omega^2+beta3/6*omega^3
%   for all frequencies, channels, modes, with the same convention
%   adopted for U.
%
%   DZB [m] is the step and may be a vector. In such a case, it means that
%   the step crosses more waveplates, and the distance traveled in each
%   waveplate is an element of the vector.
%
%   NINDEX is a vector containing the waveplate indexes corresponding to
%   the substeps of DZB.
%
%   LIN is a struct containing the unitary matrix of waveplate change of
%   basis (LIN.MATIN) and the eigenvalues (LIN.DB0), according to
%   the following diagonalization (the waveplate is Hermitian):
%
%                           Hdd = MATIN * B(w) * MATIN'
%
%   where B(w) is a diagonal matrix of eigenvalues LIN.DB0. MATIN is
%   independent of frequency.
%

u = fft(u);

for nt=1:length(dzb) % the step is made of multi-waveplates
    if lin.isuni
        u = u*conj(lin.matin(:,:,nindex(nt))); %% apply unitary matrix
        u = u.*fastexp(-bsxfun(@plus,betat,lin.db0(nindex(nt),:))*dzb(nt));
        u = u*lin.matin(:,:,nindex(nt)).'; %% return in the reference system
    else
        u = fastmultmat(u,conj(lin.matin(:,:,nindex(nt)))); %% apply unitary matrix
        u = u.*fastexp(-bsxfun(@plus,betat,lin.db0(nindex(nt),:))*dzb(nt));
        u = fastmultmat(u,lin.matin(:,:,nindex(nt)).'); %% return in the reference system
    end
end

u = ifft(u);

%--------------------------------------------------------------------------
function u = nonlinear_step(u,x,dz)

%NONLINEAR_STEP nonlinear step of the optical fiber
%   U = NONLINEAR_STEP(U,X,DZ) propagates the electric field U into the
%   nonlinear step of the NLSE.
%
%   X is a struct containing the step parameters, see FIBER.
%
%   If X.ismanakov is true the propagation is according to the Manakov
%   equation, otherwise, the CNLSE is used. DZ is the step length [m].
%
%   For instance, in single-mode and ISMANAKOV=false this function solves
%   the CNLSE with the following solution:
%
%   U = expm(-i*gam*leff*(|U|^2*[1 0;0 1] - 1/3*U'*sig3*U*sig3))
%
%   where U = [ux(k),uy(k)], k is the kth sample, and sig3 = [0 -i;i 0].
%   leff is the step effective length [m].
%

if x.alphalin == 0
    leff = dz;
else
    leff = (1-exp(-x.alphalin*dz))/x.alphalin;  % effective length [m] of dz
end
gamleff = x.gam*leff; % [1/mW]

if x.isuni %%% UNIQUE FIELD
    phi = sum(abs(u).^2,2)*gamleff; % nl phase [rad].
    expiphi = fastexp(-phi);
    u = bsxfun(@times,expiphi,u); % expiphi .* u
    
    if ~x.ismanakov && size(u,2) > 1 % CNLSE, only in dual-polarization
        s3 = 2*(real(u(:,1)).*imag(u(:,2))-imag(u(:,1)).*real(u(:,2)))*gamleff(1); % stokes comp. #3
        expigamleffs3d3 = fastexp(s3/3);
        cosphi = real(expigamleffs3d3); % (fast) cos(gamleff*s3/3)
        sinphi = imag(expigamleffs3d3); % (fast) sin(gamleff*s3/3)
        uxx = cosphi.*u(:,1) + sinphi.*u(:,2);
        uyy = -sinphi.*u(:,1) + cosphi.*u(:,2);
        u(:,1) = uxx;
        u(:,2) = uyy;
        
    end % if ismanakov
    
else % separate-field
    
        [Nfft,nfc] = size(u);
        Npol = x.isy+1;
        
        Pow = fft(real(u).^2 + imag(u).^2); % same size of u
        if x.isy
            Cxy  = fft(u(:,1:2:end).*conj(u(:,2:2:end))); 
        else
            Cxy = 0;
        end
        for nch=1:Npol:nfc % forward
            f11 = zeros(Nfft,1); f12=f11; f22=f11;
            inch = ceil(nch/(x.isy+1)); % counter single/dual-pol
            for kch=1:Npol:nfc
                ikch = ceil(kch/(x.isy+1)); % counter single/dual-pol
                if kch == nch
                    f11 = f11+gamleff(inch)*Pow(:,kch); % SPM
                    f22 = f22+gamleff(inch)*Pow(:,kch);
                    if x.isy
                        f11 = f11+gamleff(inch)*Pow(:,kch+1); % SPM of y
                        f22 = f22+gamleff(inch)*Pow(:,kch+1);
                    end
                else
                    %                     cxy  = ux(:,kch).*conj(uy(:,kch));
                    %                     Ddsp = -dsp(kch)+dsp(nch);
                    wDdsp = x.b1(:,nch)-x.b1(:,kch);
                    Hwo = x.gam(inch)*(1 - exp(-(x.alphalin - 1i*wDdsp)*dz))./...
                        (x.alphalin - 1i*wDdsp); % walk-off filter
                    
                    powx = Pow(:,kch) .* Hwo;
                    
                    f11 = f11 + 2*powx;
                    
                    if x.isy
                        powy = Pow(:,kch+1) .* Hwo;
                        f11 = f11 + powy;
                        cxy = Cxy(:,ikch) .* Hwo; % NL polarization crosstalk
                        f12 = f12 + cxy;
                        %                     f21 = conj(f12);
                        f22 = f22 + 2*powy+powx;
                    end
                    % NL cross-channel matrix between two channels:
                    % [f11, f12; conj(f12), f22]
                end
            end
            f11 = real(ifft(f11)); % f11 and f22 are for sure real
            uux = u(:,nch);
            if x.isy
                f12 = ifft(f12);
                f22 = real(ifft(f22));
                uuy = u(:,nch+1);
                %             ux(:,nch) = uux.*fastexp(-f11) - 1i*(uuy.*f12);
                %             uy(:,nch) = uuy.*fastexp(-f22) - 1i*(uux.*conj(f12));
                
                % Cayley-Hamilton expansion of expm(.) for each time
                % E.g., at time 1, if B = [f11(1) f12(1);conj(f12(1)), f22(1)],
                %   expm(-1i*B) coincides with a0(1)*eye(2)-a1(1)*B.
                % This approach is more numerical stable than pure Euler step
                sqrt12 = sqrt(4*abs(f12).^2+(f11-f22).^2);
                lam1 = 0.5*(-(f11+f22) - sqrt12); % lam1 and lam2 are the
                lam2 = 0.5*(-(f11+f22) + sqrt12); % eigenvalues of matrix exponent
                explam1 = fastexp(lam1); explam2 = fastexp(lam2);
                a0 = (lam1.*explam2-lam2.*explam1)./(lam1-lam2);
                a1 = (explam1 - explam2)./(lam1-lam2);
                u(:,nch) = (a0 -a1.*f11).*uux -a1.*f12.*uuy;
                u(:,nch+1) = (a0 -a1.*f22).*uuy -a1.*conj(f12).*uux;
            end
        end % for nch=1:nfc    
end % if isuni

%--------------------------------------------------------------------------
function [step,phimax] = firststep(u,x)

%FIRSTSTEP first step setup of the SSFM
%   [STEP,PHIMAX}=FIRSTSTEP(U,X) evaluates the first STEP [m] of the SSFM.
%   U is the electric field.
%
%   X is a struct containing several options. If X.dphi1fwm is true, the
%   first step is chosen according to the maximum FWM criterion [Mus18]. If
%   X.dphi1fwm is false, the first step is chosen according to the
%   nonlinear phase criterion [Sin03].
%
%   The step cannot, in any case, be greater than DZMAX.

global CONSTANTS
global GSTATE

if x.length == x.dzmax
    step = x.dzmax;
    phimax = Inf;
else
    if x.dphi1fwm
        if ~isfield(x,'bandwidth')
            x.bandwidth = GSTATE.FSAMPLING;
        end
        spac = x.bandwidth.*x.chlambda.^2/CONSTANTS.CLIGHT; % bandwidth in [nm]
        if x.isuni % min: worst case among spatial modes
            step = min(x.dphimax./abs(x.disp)./(2*pi*spac.*x.bandwidth*1e-3)*1e3); % [m]
        else % separate fields: FWM bandwidth is substituted by walk-off bandwidth
            step = min(x.dphimax./abs(x.disp)./(2*pi*spac*x.bandwidth*1e-3)*1e3); % [m]
        end
        if step > x.dzmax
            step = x.dzmax;
        end
        if strcmp(x.stepupd,'nlp') % nonlinear phase criterion
            if size(x.gam,1) == 1  % CNLSE
                maxgam = max(x.gam); % worst case
            end
            invLnl = max(max(abs(u).^2))*maxgam; % max of 1/Lnl [1/m]
            
            if alphalin == 0
                leff = step;
            else
                leff = (1 - exp(-x.alphalin*step))/x.alphalin;
            end
            phimax = invLnl*leff; % recalculate max nonlinear phase [rad] per step
        else
            phimax = x.dphimax;
        end
    else % nonlinear phase criterion
        step = nextstep(u,x,NaN);
        phimax = x.dphimax;
    end
end
%--------------------------------------------------------------------------
function dz=nextstep(u,x,dz_old)

%NEXTSTEP step size setup for the SSFM algorithm
%   DZ=NEXTSTEP(U,X,DZ_OLD) evaluates the step size of the SSFM.
%   U is the electric field, DZ_OLD is the step used in the
%   previous SSFM cycle.
%

if x.iscle % constant local error (CLE)
    if x.issym, q=3; else q=2; end % See [Zha05,Zha08]
    step = dz_old*exp(x.alphalin/q*dz_old); % [m]
    
else % nonlinear phase criterion
    if x.isy % max over time
        Pmax = max(abs(u(:,1:2:end)).^2 + abs(u(:,2:2:end)).^2);
    else
        Pmax = max(abs(u).^2);
    end
    invLnl = max(Pmax.*x.gam); % max over channels
    
    leff = x.dphimax/max(invLnl);  % effective length [m] of the step
    dl = x.alphalin*leff;  % ratio effective length/attenuation length
    
    if dl >= 1
        step = x.dzmax; % [m]
    else
        if x.alphalin == 0
            step = leff;
        else
            step = -1/x.alphalin*log(1-dl);
        end
    end
end
if step > x.dzmax
    dz = x.dzmax;
else
    dz = step;
end

%--------------------------------------------------------------------------
function [dzb,nindex] = checkstep(zprop,dz,lcorr)

%CHECKSTEP split the step in multi-step within a correlation length
%   [DZB,NINOUT,NINDEX]=CHECKSTEP(ZPROP,DZ,LCORR) split the step DZ [m]
%   into a  sequence of sub-steps, each of maximum length equal to LCORR.
%
%   LCORR is the waveplate (trunk) length [m], i.e., the length over which
%   the birefringence is supposed to be constant.
%
%   On output, DZB is a vector containing all the sub-steps within DZ. It
%   is sum(DZB) = DZ.
%
%   NINDEX is a vector containing the waveplate indexes corresponding to the
%   substeps of DZB.
%
%   NINOUT is a vector of length 2*length(DZB). DZB contains the indexes of
%   the steps crossed along DZ. The factor 2 accounts for input/output
%   matrices at each waveplate boundary.
%
%   For instance, take LCORR = 1 [m] and a step starting at coordinate
%   8.5 [m] with DZ=1.8 [m]. The step DZ crosses 3 waveplates and
%   DZB = [0.5, 1, 0.3].
%
%   ---------------------------------------------------------------
%            |      *     |           |    *      |
%  ZPROP:    8     8.5    9           10  10.3    11    [m]
%  TRUNK:    9            10          11          12
% NINOUT:           9    9 10       10 11  11
%
%   NINOUT in this example is [0,9,10,10,11,0]. First/last zero: the step
%   is not starting/ending over a multiple of the correlation length. For
%   the other numbers remember that trunk numbered k is within coordinates
%   (k-1)*LCORR and k*LCORR.
%
%   All lengths are in [m].

zini = zprop - dz;  % starting coordinate of the step [m]
zend = zprop;       % ending coordinate of the step [m]
% first waveplate index is 1
nini = floor(zini/lcorr)+1; % waveplate of starting coordinate
nend = ceil(zend/lcorr); % waveplate of ending coordinate

nmid = nini:nend; % waveplate indexes
if nini == nend     % start/end of the step within a waveplate
    dzsplit = dz;
    %     ntrunk = [nini nend];
    ntrunk = [0 0];
    if isa(nini,'gpuArray')
        ntrunk = gpuArray(ntrunk);
    end
else                % multi-waveplate step
    zmid = lcorr*nmid(1:end-1); % waveplate mid-coordinates
    dzsplit = [zmid(1)-zini,diff(zmid),zend-zmid(end)];
    ntmp = [nmid;nmid];
    ntrunk = ntmp(:);
end

if (nini-1)*lcorr ~= zini, ntrunk(1) = 0; end % starting coord. over a trunk
if (nend-1)*lcorr ~= zend, ntrunk(end) = 0; end % ending coord. over a trunk

dzb = dzsplit(dzsplit ~= 0); % remove zero-length steps
nindex = nmid(dzsplit ~= 0); % remove zero-length steps
% ninout = ntrunk(1:2*length(dzb)); % trunk coordinates. 0: do not apply
% if ~mod(zini,lcorr), ninout(1) = nini; end % only at beginning of fiber
% if ~mod(zend,lcorr), ninout(end) = nend; end % only at beginning of fiber

% fprintf('zini = %.2f zend = %.2f lcorr = %.1f\n',zini,zend,lcorr);
% fprintf('nini = %d. nend = %d\n',nini,nend);
% fprintf('dzb = ');fprintf('%.2f ',dzb);fprintf('\n')
% fprintf('dzsplit = ');fprintf('%.2f ',dzsplit);fprintf('\n')
% fprintf('ntrunk = ');fprintf('%d ',ntrunk); fprintf('\n');
% fprintf('ninout = ');fprintf('%d ',ninout); fprintf('\n\n');

%--------------------------------------------------------------------------
function u=pol2vec(E)

%POL2VEC concatenate polarizations
%   U=POL2VEC(E) creates matrix U by concatenating the polarizations
%   contained in E. The size of the output is the following:
%
%   unique field   : U = (Nfft,Npol)
%   separate fields: U = (Nfft,Npol,Nchannel)
%
%   with Nfft the number of samples, Npol the number of polarizations,
%   Nchannel the number of channels.
%
%   POL2VEC is the opposite of VEC2POL.

u(:,1,:) = E.fieldx;
if isfield(E,'fieldy')
    u(:,2,:) = E.fieldy;
end

%--------------------------------------------------------------------------
function E = vec2pol(u)

%VEC2POL restore polarizations
%   E=VEC2POL(U) is the undo of POL2VEC.

E.fieldx = u(:,1,:);
if size(u,2) > 1
    E.fieldy = u(:,2,:);
end

%--------------------------------------------------------------------------
function [U,S] = eigendec(x)
%EIGENDEC - spectral decomposition of a random waveplate
%   [U,S] = EIGENDEC(X) returns in S a vector containing the eigenvalues
%   of the 2 x 2 matrix describing a random waveplate. U is the matrix
%   with the corresponding eigenvectors on columns.
%
%   The matrix describing a waveplate of length h is Hermitian,
%   hence with spectral decomposition:
%
%       U * diag(fastexp(S*h)) * U'
%
%   X is a struct whose fields are the following:
%
%       X.coupling:   linear coupling type. Valid options are 'pol', i.e.,
%                     coupling occurs only between polarizations, or
%                     'none'.
%
%		X.beatlength: beat-length [m] between polarizations. If dbeta is
%                     the differential beta between polarizations we have:
%
%                       X.beatlength = 2*pi/<Dbeta>
%
%					  with <.> ensemble average. When the coupling is
%					  active, polarizations couples according to the
%					  Langevin model of [Wai96]. The Stokes vector
%   			      describing the coupling matrix is described by zero
%					  average Gaussian random variables, independent
%					  waveplate-to-waveplate.
%                     The variance of such variables is set to have an
%					  average beat length equal to x.beatlength.
%
%   References:
%
%   [Wai96] P. K. A. Wai and C. R. Menyuk, “Polarization mode dispersion,
%       decorrelation, and diffusion in optical fibers with randomly
%       varying birefringence," J. Lightw. Technol., vol. 14, no. 2, pp.
%       148–157, 1996.
%
%   [Ant13] C. Antonelli, A. Mecozzi, M. Shtaif, and P. J. Winzer, “Random
%       coupling between groups of degenerate fiber modes in mode
%       multiplexed transmission," Opt. Express, vol. 21, no. 8, p. 9484,
%       2013.
%
%   Author: Paolo Serena, 2018
%   University of Parma, Italy

if strcmpi(x.coupling, '') || strcmpi(x.coupling, 'none')
    
    U = eye(2);
    S = zeros(1,2); % eigenvalues
    
elseif strcmpi(x.coupling, 'pol') % couple only polarizations
    
    U = randU(2); % random unitary matrix (Haar matrix)
    q = sqrt(pi^3/2)/x.beatlength*randn(1,3); % [1/m]
    w = sqrt(sum(q.^2,2)); % Maxwellian distribution. deltabeta0 [1/m]
    S = [w -w]/2;
else
    error('Unknown coupling method.');
end

%--------------------------------------------------------------------------
function u=fastmultmat(u,M)

% https://stackoverflow.com/questions/40508500/how-to-dynamically-reshape-matrix-block-wise
if isscalar(M)
    % nothing to do
else
    [nr,nc]=size(u);
    u=reshape(permute(reshape(u,nr,2,[]),[1,3,2]),[],2); % 2: polarizations
    u = u * M;
    u=reshape(permute(reshape(u,nr,nc/2,[]),[1 3 2]),nr,[]);
end