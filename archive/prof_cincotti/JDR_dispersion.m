clear
close all
uscita1=8;
uscita2=5;
uscita3=6;
M=8; % numero di porte del E/D
Numerocampioni=(M+1)*3000; % serve per i grafici 
FSR=200; %bit rate del codice in GHz 

% parametro di dispersione per L Km
L=10;
D=17*L;% ps/nm
lambda=1550; %wavelength [nm]
c=299792458; %speed of light [nm/ns]

tmax=10*(M+1)/FSR; %intervallo di tempo max 
tempo=linspace(0,tmax,Numerocampioni); % asse dei tempi
tempoconv=linspace(0,2*tmax,2*Numerocampioni-1); % asse dei tempi per graficare le convluzioni

f=(Numerocampioni-1)/(2*tmax)*linspace(-1,1,Numerocampioni); %frequenza misurata in GHz
coefficientiHAD=hadamard(M);
 %%%% EFFETTO DISPERSIONE
    % affinche le dimensioni tornino D deve essere misurato in ns/nm,
    % moltiplico per 10^-3
filter=exp(-1i*pi*D*10^-3*lambda^2*f.^2/c);% transfer function of the dispersive fiber

% full width half maximum of the Gaussian pulse from the laser
fwhm=0.003; % [ns]
deltat=fwhm/(2*sqrt(log(2))); % variance of the Gaussian pulse
TRA=1+(M-1)/2;



% codewords nel dominio della frequenza
spettroFT=zeros(M,Numerocampioni); 
spettroGM=zeros(M,Numerocampioni); 

for k=1:M %number of ports
            for j= 1:M
           akj=exp(-1i*2*pi*k*j/M)/M;                
          spettroFT(k,:)=spettroFT(k,:)+akj*exp(-1i*2*pi*f*j/FSR).*exp(-2*deltat^2*f.^2*pi^2);
          spettroGM(k,:)=spettroGM(k,:)+coefficientiHAD(k,j)/M*exp(-1i*2*pi*f*j/FSR).*exp(-2*deltat^2*f.^2*pi^2);
          % codewords nel dominio del tempo
            end
            signalFT(k,:)=ifft(ifftshift(spettroFT(k,:)));
          signalGM(k,:)=ifft(ifftshift(spettroGM(k,:)));
        spettroFT_disp(k,:)=spettroFT(k,:).*filter;
        spettroGM_disp(k,:)=spettroGM(k,:).*filter;
signalFT_disp(k,:)=ifft(ifftshift(spettroFT_disp(k,:)));
signalGM_disp(k,:)=ifft(ifftshift(spettroGM_disp(k,:)));
end
signalFT=signalFT/max(max(signalFT));
signalGM=signalGM/max(max(signalGM));
signalFT_disp=signalFT_disp/max(max(signalFT_disp));
signalGM_disp=signalGM_disp/max(max(signalGM_disp));

figure (1)
subplot(2,2,1), plot(tempo,real(signalFT(uscita1,:)),'-k','LineWidth', 1.5)
axis([0 0.1 0 1.1])
title 'FM codewords'
subplot(2,2,2),plot(tempo,real(signalGM(uscita1,:)),'-k','LineWidth', 1.5)
axis([0 0.1 -1.1 1.1])
title 'GM codewords'
subplot(2,2,3), plot(tempo,abs(signalFT_disp(uscita1,:)),'-r','LineWidth', 1.5)
axis([0 0.1 0 1.1])
xlabel('time [ns]')
title 'FM codewords after dispersive fiber'
subplot(2,2,4),plot(tempo,signalGM_disp(uscita1,:),'-r','LineWidth', 1.5)
axis([0 0.1 -1.1 1.1])
xlabel('time [ns]')
title 'GM codewords after dispersive fiber'


  
 figure (2)
subplot(1,2,1), plot(f,abs(spettroFT(uscita1,:)),'-k','LineWidth', 1.5)
axis([-FSR/2 FSR/2 0 1.1])
xlabel('frequency [GHz]')
ylabel('spectrum')
title 'FM codewords'
subplot(1,2,2),plot(f,abs(spettroGM(uscita1,:)),'-k','LineWidth', 1.5)
axis([-FSR/2 FSR/2 0 1.1])
title 'GM codewords'
xlabel('frequency [GHz]')
ylabel('spectrum')




if 0
% codewords nel dominio del tempo SERVE SOLO PER CONTROLLO
signalFT=zeros(M,Numerocampioni); 
signalGM=zeros(M,Numerocampioni); 
for k=1:M %number of ports
            for j= 1:M
           akj=exp(-1i*2*pi*k*j/M);                
           % per graficarli al centro, li ritardo di tmax/2
          signalFT(k,:)=signalFT(k,:)+akj*exp(-(tempo-j/FSR-tmax/2).^2/(2*deltat^2));
          signalGM(k,:)=signalGM(k,:)+coefficientiHAD(k,j)*exp(-(tempo-j/FSR-tmax/2).^2/(2*deltat^2));
          end
    end
    signalFT(k,:)=signalFT(k,:)/max(abs(signalFT(k,:)));
end

  % calcolo le correlazione per FT
 
  label1_FT=zeros(1,Numerocampioni);
  label2_FT=zeros(1,Numerocampioni);
  label3_FT=zeros(1,Numerocampioni);
for k  = 1:M 
   label1_FT(Numerocampioni*k/(FSR*tmax))=exp(-i*pi*2*k*(uscita1)/M); 
   label2_FT(Numerocampioni*k/(FSR*tmax))=exp(-i*pi*2*k*(uscita2)/M);
   label3_FT(Numerocampioni*k/(FSR*tmax))=exp(-i*pi*2*k*(uscita3)/M);
end
  
      

 autocorrelazione_FT=conv(signalFT(uscita1,:),label1_FT)/M;
 correlazione_FT=conv(signalFT(uscita1,:),label2_FT)/M;
 correlazione1_FT=conv(signalFT(uscita1,:),label3_FT)/M;

autocorrelazione_FT_disp=conv(signalFT_disp(uscita1,:),label1_FT)/M;
 correlazione_FT_disp=conv(signalFT_disp(uscita1,:),label2_FT)/M;
 correlazione1_FT_disp=conv(signalFT_disp(uscita1,:),label3_FT)/M;

CCP= max(abs(correlazione_FT));
ACP= max(abs(autocorrelazione_FT));


figure (4)
subplot(3,2,1),plot(tempoconv,abs(autocorrelazione_FT),'-k','LineWidth',1.2)
axis([0 0.3 -1.1 1.1])
 ylabel('Matched port' )
 title 'Correlation FT codewords'
subplot(3,2,3),plot(tempoconv,abs(correlazione_FT),'-k','LineWidth',1.2)
axis([0 0.3 -1.1 1.1])
 ylabel('Unmatched port ')
subplot(3,2,5),plot(tempoconv,abs(correlazione1_FT),'-k','LineWidth',1.2)
xlabel('time [ns]')
axis([0 0.3 -1.1 1.1])
 ylabel('Unmatched port')
subplot(3,2,2),plot(tempoconv,abs(autocorrelazione_FT_disp),'-r','LineWidth',1.2)
axis([0 0.3 -1.1 1.1])
 ylabel('Matched port' )
   title 'Correlation FT codewords after dispersion'
subplot(3,2,4),plot(tempoconv,abs(correlazione_FT_disp),'-r','LineWidth',1.2)
axis([0 0.3 -1.1 1.1])
 ylabel('Unmatched port ')
subplot(3,2,6),plot(tempoconv,abs(correlazione1_FT_disp),'-r','LineWidth',1.2)
axis([0 0.3 -1.1 1.1])
xlabel('time [ns]')
ylabel('Unmatched port')


% calcolo le correlazione per GM anche se poi devo usare gli switch (cioè
% metto a zero il segnale prima e dopo)
      
  label1_GM=zeros(1,Numerocampioni);
  label2_GM=zeros(1,Numerocampioni);
  label3_GM=zeros(1,Numerocampioni);
for k  = 1:M 
   label1_GM(Numerocampioni*k/(FSR*tmax))=coefficientiHAD(uscita1,k); 
   label2_GM(Numerocampioni*k/(FSR*tmax))=coefficientiHAD(uscita2,k); 
   label3_GM(Numerocampioni*k/(FSR*tmax))=coefficientiHAD(uscita3,k); 
end
 autocorrelazione_GM=conv(signalGM(uscita1,:),label1_GM)/M;
 correlazione_GM=conv(signalGM(uscita1,:),label2_GM)/M;
 correlazione1_GM=conv(signalGM(uscita1,:),label3_GM)/M;

autocorrelazione_GM_disp=conv(signalGM_disp(uscita1,:),label1_GM)/M;
 correlazione_GM_disp=conv(signalGM_disp(uscita1,:),label2_GM)/M;
 correlazione1_GM_disp=conv(signalGM_disp(uscita1,:),label3_GM)/M;
% adesso uso gli switch ovvero cancello il segnale
% 
minimo=(M)/FSR*(2*Numerocampioni-1)/(2*tmax);
massimo=(M+2)/FSR*(2*Numerocampioni-1)/(2*tmax);
autocorrelazione_GM(1:minimo)=0;
autocorrelazione_GM(massimo:2*Numerocampioni-1)=0;
correlazione_GM(1:minimo)=0;
correlazione_GM(massimo:2*Numerocampioni-1)=0;
correlazione1_GM(1:minimo)=0;
correlazione1_GM(massimo:2*Numerocampioni-1)=0;
autocorrelazione_GM_disp(1:minimo)=0;
autocorrelazione_GM_disp(massimo:2*Numerocampioni-1)=0;
correlazione_GM_disp(1:minimo)=0;
correlazione_GM_disp(massimo:2*Numerocampioni-1)=0;
correlazione1_GM_disp(1:minimo)=0;
correlazione1_GM_disp(massimo:2*Numerocampioni-1)=0;
figure (5)
subplot(3,2,1),plot(tempoconv,abs(autocorrelazione_GM),'-k','LineWidth',1)
axis([0 0.3 -1.1 1.1])
 ylabel('Matched port' )
  title 'Processing GM codewords'

subplot(3,2,3),plot(tempoconv,abs(correlazione_GM),'-k','LineWidth',1)
axis([0 0.3 -1.1 1.1])
 ylabel('Unmatched port ')
subplot(3,2,5),plot(tempoconv,abs(correlazione1_GM),'-k','LineWidth',1)
xlabel('time [ns]')
axis([0 0.3 -1.1 1.1])
 ylabel('Unmatched port')
subplot(3,2,2),plot(tempoconv,abs(autocorrelazione_GM_disp),'-r','LineWidth',1)
axis([0 0.3 -1.1 1.1])
  title 'Processing GM codewords after dispersion'

 ylabel('Matched port' )
subplot(3,2,4),plot(tempoconv,abs(correlazione_GM_disp),'-r','LineWidth',1)
axis([0 0.3 -1.1 1.1])
 ylabel('Unmatched port ')
subplot(3,2,6),plot(tempoconv,abs(correlazione1_GM_disp),'-r','LineWidth',1)
axis([0 0.3 -1.1 1.1])
xlabel('time [ns]')
ylabel('Unmatched port')



