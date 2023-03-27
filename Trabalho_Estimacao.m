%% Trabalho de Estimação Linear
% Estimação LMS Equalizador
%
% Autor: Ecard e Eric
%
%

clear; close all; clc

% BlocoTotal = 1; %ensemble
n = 10000; % Amostras
deltaBloco = 0.01*n; % Tamanho do bloco
C = 1; % Potencia do ruido

%% Gerando Sinal
X = randn(n,1);
% S = double(X>(.5));
% S(S==0) = -1;
S = sign(tanh(X));
S_bloco = reshape(S, [deltaBloco n/deltaBloco]);

% Canal H
h = [randn(.4*deltaBloco,1);zeros(.6*deltaBloco,1)];
H = toeplitz(h);
H = tril(H);

B = zeros(deltaBloco, deltaBloco);

for indBloco = 1:(n/deltaBloco)
    
    v = C*randn(deltaBloco,1); % ruido
    y = H*S_bloco + v;
    
    Rss = xcorr(S_bloco(:,indBloco));
    Rss = toeplitz(Rss(deltaBloco:end));

    Rvv= xcorr(v);
    Rvv = toeplitz(Rvv(deltaBloco:end));
    
    % A mais simples (TESTE)
    % Porque o resultado é melhor? É mais custoso?
    w0 =  H'/(Rvv);
    w = inv(inv(Rss) + w0*H)*H'/(Rvv); %pag 82
    x_hat(:,indBloco) = w*y(:,indBloco);

    % Comparador
    S_hat(:,indBloco) = sign(tanh(x_hat(:,indBloco)));
    
    % Erro
    e_til_teste_amostra(:,indBloco) = (S_bloco(:,indBloco) ~= S_hat(:,indBloco));
    
    % Fazendo B = 0:
    aux = H'/(Rvv);
    aux2 = inv(Rss)+aux*H;
    P = inv(aux2);
    [L,D] = ldl(P);
    G = D*L'*H'/(Rvv);
    
    x_est(:,indBloco) = G*y(:,indBloco);
    % Comparador
    S_est(:,indBloco) = sign(tanh(x_est(:,indBloco)));
    
    % Erro
    e_til_amostra(:,indBloco) = (S_bloco(:,indBloco) ~= S_est(:,indBloco));
    
    % pag 99
    % Fazendo B diferente de zero
    % Atualização a cada bloco utilizando a estimação do bloco anterior
    
    B = inv(L)-eye(deltaBloco); %inv(L) - I
    if indBloco > 1
        aux(:,indBloco) = B*S_est2(:,indBloco-1);
        x_est2(:,indBloco) = x_est(:,indBloco) - aux(:,indBloco);
        % Comparador
        S_est2(:,indBloco) = sign(tanh(x_est2(:,indBloco)));
        S_est2(S_est2(:,indBloco)==0,indBloco) = -1;
    else
        S_est2(:,indBloco) = S_est(:,indBloco);
    end
    
    %     if indBloco > 1
    %         B(indBloco,:) = B_aux(indBloco,:);
    %         aux = B(indBloco,:).*S_est2(:,indBloco-1);
    %         x_est2(:,indBloco) = x_est(:,indBloco) - aux(:,indBloco);
    %         % Comparador
    %         S_est2(:,indBloco) = double(x_est2(:,indBloco)>0);
    %         S_est2(S_est2(:,indBloco)==0,indBloco) = -1;
    %     else
    %         S_est2(:,indBloco) = S_est(:,indBloco);
    %     end
    
    % Erro
    e_til2_amostra(:,indBloco) = (S_bloco(:,indBloco) ~= S_est2(:,indBloco));
end


%% Resultados

e_total_til_teste = sum(sum(e_til_teste_amostra))/(deltaBloco*indBloco);
e_total_til = sum(sum(e_til_amostra))/(deltaBloco*indBloco);
e_total_til2 = sum(sum(e_til2_amostra))/(deltaBloco*indBloco);

disp(['ERRO X Estimado Teste = ' num2str(e_total_til_teste*100, '%f') '%'])
disp(['ERRO X Estimado (B=0) = ' num2str(e_total_til*100, '%f') '%'])
disp(['ERRO X Estimado (B=/=0) = ' num2str(e_total_til2*100, '%f') '%'])