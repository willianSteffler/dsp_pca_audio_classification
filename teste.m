clc;clear all; close all;

folder = 'audio';
samples=dir(folder);
totalFiles = numel(samples);
totalFiles =totalFiles - 2 %remove . e ..
menorAmostra = 1000000000;
qtdePrimeiroAndar = 104;
qtdeSegundoAndar = 42;

for i = 1 : totalFiles
    Filepath{i,1}=strcat(folder,'/',samples(i+2).name);
    [y,Fs] = audioread(Filepath{i,1});
    %size(y)
    if (size(y,1)<menorAmostra)
        menorAmostra = size(y,1);
        indice_menor_amostra = i;
    end
end

fprintf('% % % % % % % % % % % % % % % % % % % % % #\n1')
X = zeros(menorAmostra,totalFiles);
%reduz dimensão das amostras
for i = 1 : totalFiles
    Filepath{i,1};
    [y,Fs] = audioread(Filepath{i,1});
    
    dim = 1:1:length(y);                        %espaço original
    tgt = 1:length(y)/menorAmostra:length(y);   %espaço reduzido
    S = interp1(dim,y(:,1),tgt);                %sinal interpolado
    X(:,i) = S;
end

%size(X)

mu = mean(X,2);  %Média para cada ponto P de todas as amostras (2)
mu_1oAndar = mean(X(:,1:qtdePrimeiroAndar),2);
mu_2oAndar = mean(X(:,qtdePrimeiroAndar+1:qtdePrimeiroAndar+qtdeSegundoAndar),2);

%wt= V'*(K-mu) Padroes Centralizados (3)
Phi = X;
for i = 1:totalFiles,
    Phi(:,i) = Phi(:,i) - mu;
end

Cxx = Phi'*Phi/totalFiles;   %matriz de covariancia

%size(Cxx)

[U, Lambda] = eig(Cxx);
U=U(:,1:size(U,2));

V= Phi*U;
for i = 1 : size(V,2)
    V(:,i) = V(:,i) / norm(V(:,i));
end


S = X(:,140);               %Amostra de teste
Sm = S - mu;
S1 = S - mu_1oAndar;
S2 = S - mu_2oAndar;



% w = V'*(S-mu);
w = V'*Sm;
(Sm' * Sm) - (w' * w)

w = V'*S1;
% w = V'*(S-mu_1oAndar);
(S1' * S1) - (w' * w)

w = V'*S2;
% w = V'*(S-mu_2oAndar);
(S2' * S2) - (w' * w)