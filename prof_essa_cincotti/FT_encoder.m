clear
close all

colormap("prism")

uscita1=2;
uscita2=3;
uscita3=4;
M=8; % numero di porte del E/D
alfa=0; %rapporto w/di dello slab diffraction deve essere <1 é un parametro dell' E/D che per il momento trascuriamo
%
Numerocampioni=(M+1)*3000; % serve per i grafici 
FSR=200; %bit rate del codice in GHz 

tmax=10*(M+1)/FSR; %intervallo di tempo max 
tempo=linspace(0,tmax,Numerocampioni); % asse dei tempi
tempoconv=linspace(0,2*tmax,2*Numerocampioni-1); % asse dei tempi per graficare le convluzioni

f=(Numerocampioni-1)/(2*tmax)*linspace(-1,1,Numerocampioni); %frequenza misurata in GHz

  spettroFT=zeros(M, Numerocampioni); % spettro dei codici FT
  spettroHAD=zeros(M, Numerocampioni); % spettro dei codici Hadamard
  coefficientiHAD=hadamard(M);
  %spettroslab=zeros(Nprimo, Numerocampioni); % spettro dei codici considerando l'effetto dello slab diffraction
  for n=1:M
      for k=1: M
          spettroFT(n,:)=spettroFT(n,:)+exp(-i*2*pi*n*k/M).*exp(-1i*2*pi*f*k/FSR)/M;
          %spettroslab(n,:)=spettroslab(n,:)+exp(-i*2*pi*TF*n*k/M).*exp(-i*2*pi*f*k/FSR).*exp(-(pi*alfa*(k-M/2)/M)^2)/M;  %tiene in considerazione lo slab diffraction;
        spettroHAD(n,:)=spettroHAD(n,:)+coefficientiHAD(n,k).*exp(-1i*2*pi*f*k/FSR)/M;
      end
  end
  

  colors = prism(M);  % Get M colors from the 'prism' colormap

% Use these colors in your plots
figure
hold on
for n = 1:M
    plot(f, abs(spettroFT(n, :)), 'Color', colors(n, :), "LineWidth", 1.5)
end
axis([-100 100 0 1.1])
xlabel('f [GHz]')
hold off


plot(f,abs(spettroFT), "LineWidth", 1.5)
axis([-100 100 0 1.1])
 xlabel('f [GHz]')


% rappresentazione dei FT codewords nel tempo
 
  label1=zeros(1,Numerocampioni);
  label2=zeros(1,Numerocampioni);
  label3=zeros(1,Numerocampioni);
  %label1slab=zeros(1,Numerocampioni);
  %label2slab=zeros(1,Numerocampioni);
  %label3slab=zeros(1,Numerocampioni);
for k  = 1:M 
   label1(Numerocampioni*k/(FSR*tmax))=exp(-i*pi*2*k*(uscita1)/M); 
   label2(Numerocampioni*k/(FSR*tmax))=exp(-i*pi*2*k*(uscita2)/M);
   label3(Numerocampioni*k/(FSR*tmax))=exp(-i*pi*2*k*(uscita3)/M);
   %label1slab(Numerocampioni*k/(FSR*tmax))=exp(-i*pi*2*k*(uscita1)/M).*exp(-(pi*alfa*(k-M/2)/M)^2); %tiene in considerazione lo slab diffraction;
   %label2slab(Numerocampioni*k/(FSR*tmax))=exp(-i*pi*2*k*(uscita2)/M).*exp(-(pi*alfa*(k-M/2)/M)^2);
   %label3slab(Numerocampioni*k/(FSR*tmax))=exp(-i*pi*2*k*(uscita3)/M).*exp(-(pi*alfa*(k-M/2)/M)^2);
   
    end
  
  
      

 autocorrelazione=conv(label1,label1)/M;
 correlazione=conv(label1,label2)/M;
 correlazione1=conv(label1,label3)/M;
%  autocorrelazioneslab=conv(label1slab,label1slab);%tiene in considerazione lo slab diffraction;
% correlazioneslab=conv(label1slab,label2slab);
% correlazione1slab=conv(label1slab,label3slab);

CCP= max(abs(correlazione))
ACP= max(abs(autocorrelazione));
%CCPslab= max(abs(correlazioneslab)); %tiene in considerazione lo slab diffraction;
%ACPslab= max(abs(autocorrelazioneslab));

figure
subplot(4,1,1),plot(tempo,abs(label1),'-k','LineWidth',1.2)
axis([0 tmax/10 0 1.1])
 ylabel('Intensity of FT codewords'), xlabel('t (ns)')
subplot(4,1,2),plot(tempoconv,abs(autocorrelazione),'-k','LineWidth',1.2)
axis([0 2*tmax/10 0 1.1])
 ylabel('Matched port' ), xlabel('t (ns)')
subplot(4,1,3),plot(tempoconv,abs(correlazione),'-k','LineWidth',1.2)
axis([0 2*tmax/10 0 1.1])
 ylabel('Unmatched port '), xlabel('t (ns)')
subplot(4,1,4),plot(tempoconv,abs(correlazione1),'-k','LineWidth',1.2)
axis([0 2*tmax/10 0 1.1])
 ylabel('Unmatched port'), xlabel('t (ns)')


 figure 
plot(f,abs(spettroHAD), "LineWidth", 1.5)
axis([-100 100 0 1.1])
 xlabel('f [GHz]')

figure 
plot(f,abs(spettroHAD(M,:)),'-k')
hold on
plot(f,abs(spettroFT(M,:)),'-r')
axis([-100 100 0 1.1])
ylabel('Comparison of FM and GM codeword spectra'), xlabel('f [GHz]')

 
