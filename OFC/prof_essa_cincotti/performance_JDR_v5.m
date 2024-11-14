clear all
close all
M=16;
nR=[0.01 0.1 1]';
n_H=zeros(3,1);
SNRdB = 0:0.1:20; % SNR in dB array of values [0, 1, 2, ..., 20]
SNR = 10 .^ (SNRdB / 10); % Convert SNR from dB to linear scale
for m=1:3
n_H=nR(m)./SNR; % rumore Hadamard

AUTO(m)=nR(m)*(2*M.^2+1)./(3*M);
n_F(m,:)=nR(m)./(SNR.*M); % rumore Fourier

   for k=1:M-1
CROSS(m,k)=nR(m)/(M.*(sin(pi*k./M)).^2);
   end

% Hadamard (Green Machine) theoretical BER calculation
SER_H(m,:)=exp(- (M * nR(m) + n_H));

for k=1:M-1
     SER_H(m,:)=SER_H(m,:)+(1 - exp(-n_H));
end
 SER_H(m,:)=SER_H(m,:)./M;

 % Fourier Machine theoretical BER calculation
SER_F(m,:)=exp(- (AUTO(m) + n_F(m,:)));
for k=1:M-1
     SER_F(m,:)=SER_F(m,:)+(1 - exp(-CROSS(m,k)-n_F(m,:)));
end


 SER_F(m,:)=SER_F(m,:)./M;
 BER_PSK(m,:)=exp(-nR(m)+n_H)/2+(1 - exp(-n_H))/2;
 SER_PSK(m,:)=1-(1-BER_PSK(m,:)).^log2(M);
end


figure (1)
semilogy(SNRdB,SER_H, ':')
hold on
semilogy(SNRdB,SER_F)
semilogy(SNRdB,SER_PSK,'--')
grid on
xlabel('SNR [dB]')
ylabel('SER')
legend ('Hadamard nR=0.01','Hadamard nR=0.1','Hadamard nR=1', 'Fourier nR=0.01','Fourier nR=0.1','Fourier nR=1','BPSK nR=0.01','BPSK nR=0.1','BPSK nR=1')

