%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Questo script genera e visualizza le parole di codice Hadamard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Numero di bit
n = 4;

% Numero di codewords di Hadamard
m = 2^n;

% Generazione della matrice Hadamard per n bit (m x m)
H = hadamard(m);

% Visualizzazione delle parole di codice
disp('Parole di codice Hadamard:');
disp(H);


% Dimostrazione empirica che le righe della matrice Hadamard sono ortogonali
% Con m parole si devono conforntare m!/2*(m-2)! coppie di righe
disp('Dimostrazione empirica che le righe della matrice Hadamard sono ortogonali:');
dot_res = zeros(1,m*(m-1)/2);
for i = 1:m
    for j = i+1:m
        %dot_res((i-1)*m+j) = dot(H(i,:),H(j,:));
    end
end
disp(dot_res);

% Visualizzazione grafica
figure;
imagesc(H); % Crea un'immagine della matrice
colormap(gray); % Scala di grigi
title('Visualizzazione delle parole di codice Hadamard');

%add labels, matrix lines and values for each bit
for i = 1:m
    for j = 1:m
        if H(i,j) == 1
            text(j,i,'1','HorizontalAlignment','center','VerticalAlignment','middle');
        else
            text(j,i,'-1','HorizontalAlignment','center','VerticalAlignment','middle', 'Color', 'white');
        end
    end
end

xlabel('Bit di codifica');
ylabel('codeword');


