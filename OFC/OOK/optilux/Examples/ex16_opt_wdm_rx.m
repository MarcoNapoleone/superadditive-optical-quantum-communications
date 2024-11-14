%   ex16: Generation and detection of a wavelength-division multiplexing
%   (WDM) signal in back-to-back. Moreover, a practice with saving to file
%   and print info on the screen.
%

clear
clc

%% Global parameters

Nsymb = 1024;           % number of symbols
Nt = 32;                % number of discrete points per symbol
issave = false;         % save results to file?

%% Tx parameters

Nch = 5;                % number of WDM channels
Ncut = ceil(Nch/2);     % channel under test (CUT) index
symbrate = 32;          % symbol rate [Gbaud]
tx.rolloff = 0.1;       % pulse roll-off
tx.emph = 'asin';       % digital-premphasis type
modfor = '16qam';       % modulation format
PdBm = 3;               % power [dBm]
lam = 1550;             % central carrier wavelength [nm]
deltaf = 37.5;          % channel spacing [GHz]

%% Rx parameters
rx.modformat = modfor;      % modulation format
rx.sync.type = 'da';        % time-recovery method
rx.oftype = 'gauss';        % optical filter type
rx.obw = Inf;               % optical filter bandwidth normalized to symbrate
rx.eftype = 'rootrc';       % optical filter type
rx.ebw = 0.5;               % electrical filter bandwidth normalized to symbrate
rx.epar = tx.rolloff;
rx.type = 'bin';            % binary pattern

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Init

% Conversions
spac = deltaf*lam^2/299792458;  % channel spacing [nm]
Plin = 10.^(PdBm/10);           % [mW]

% init global variables
Nsamp = Nsymb*Nt;        % overall number of samples
fs = symbrate*Nt;        % sampling rate [GHz]
inigstate(Nsamp,fs);     % initialize global variables: Nsamp and fs.

% init for saving
thisfile = textread([mfilename,'.m'], '%s', 'delimiter', '\n', 'whitespace','');
try % for LinuX
    [~,thisnohup]=system(['ls -drt ',mfilename,'*nohup | tail -1']);
catch
    thisnohup = '';
end

file_name = ['Nch=',num2str(Nch),...
    '_',modfor,...
    '_PdBm=',num2str(PdBm),...
    '_Df=',num2str(deltaf)]; % save info in this file (.mat or .txt)

outdir = mfilename; % same name of this file
[outfile,fid] = preparefile(outdir,file_name,issave); % get output file
% fid = 1: do not write on file, otherwise write on outfile.txt

%% Print info on screen
% Personalize to fit your needs. Printing info is important, do not
% underrate it!
hostn = char(java.net.InetAddress.getLocalHost.getHostName); % host name
for m=unique([1 fid]) % if output file exists, write both on file and command window
    fprintf(m,'Start (%s):\t%s.\n\n',hostn, datestr(now));
    [~,outname] = fileparts(outfile);
    fprintf(m,'Matlab file:\t%s\nOutput dir:\t%s\n',mfilename,outdir);
    if fid ~=1, fprintf(m,'Output file:\t%s\n\n',outname); else, fprintf(m,'\n\n'); end
    fprintf(m,'-------Tx-------\n');
    fprintf(m,'Nsymb = %d. Samp x symb = %d. Sampling freq = %.1f [GHz]\n',...
        Nsymb, Nt, fs);
    if Nch == 1, chstr='1ch'; else, chstr=sprintf('WDM: %dx%.1f=%.1f [GHz]',Nch,deltaf,Nch*deltaf);end
    fprintf(m,'%s @ %.1f [Gbd], %s, P = %.2f [dBm]\n',modfor,symbrate,chstr,PdBm);
end

%% Tx side

lamv = getlambda(lam,spac,Nch); % carrier wavelengths [nm]

for k=1:Nch
    rng(1e5+k-ceil(Nch/2)); % ceil(..): keep seed for increasing Nch
    E = lasersource(Plin,lamv(k));
    
    [Ex,Ey] = pbs(E); % split in two orthogonal polarizations
    patx(:,k) = pattern2(Nsymb,'rand',struct('format',modfor));
    paty(:,k) = pattern2(Nsymb,'rand',struct('format',modfor));
    
    [elecx, normx] = digitalmod(patx(:,k),modfor,symbrate,'rootrc',tx);
    [elecy, normy] = digitalmod(paty(:,k),modfor,symbrate,'rootrc',tx);

    Ex  = iqmodulator(Ex, elecx,struct('norm',normx));
    Ey  = iqmodulator(Ey, elecy,struct('norm',normy));
    
    E = pbc(Ex,Ey); % combine creating a PDM signal
    
    if k == 1
        Em = E; % nothing to multiplex together
    else
        Em = multiplexer(Em,E,struct('lambda',lamv(Ncut)));
    end
end

%% Receiver

for k=1:Nch % detect channel-by-channel
    rsig = rxfrontend(Em,lamv(k),symbrate,rx); % front-end kth-channel
    akhat = rxdsp(rsig,symbrate,[patx(:,k) paty(:,k)],rx);
    pathat = samp2pat(akhat,modfor,rx); % Rx bits
    SNRdBb2bhat(k) = samp2snr([patx(:,k) paty(:,k)],akhat,modfor);
end

% The SNR is not Inf because of numerical noise
fprintf('\nCh:\t');for k=1:Nch, fprintf('%d\t',k);end; fprintf('\n')
fprintf('SNR:\t');for k=1:Nch, fprintf('%.2f\t',SNRdBb2bhat(k));end; fprintf(' [dB] \n')

%% Close simulation

for m=unique([1 fid]), fprintf(m,'\n\nEND of simulation\n\n'); end
if fid ~= 1, fclose(fid); end % close output file
