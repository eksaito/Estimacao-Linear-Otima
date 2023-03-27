%% Trabalho de Estimação Linear
% Estimação LMS Equalizador
%
% Autor: Ecard e Eric
%
% 0º loop -> delta
% 1º loop -> SNR
% 2º loop -> ensemble
%               -> randn H
% 3º loop -> estimação de blocos
%               -> randn X

clear;  clc

n_amostras = 16; % Amostras
n_H = 9; % Coef. canais

aux = true;
delta = 4:8;

for ind = 1:length(delta)+1 %delta
    for ind_SNR = 1:20 %1º loop -> SNR
        for n_loop = 1: 10000 %2º loop -> ensemble
            C = sqrt(ind_SNR);
            
            X = randn(n_amostras,1);
            S_bloco = sign(tanh(X));
            %S_bloco = S(1:n_amostras);
            
            h = normalize(randn(n_H,1));
            H = convmtx(h,n_amostras);
            
            if ind == 1
                H1 = H;
            else
                H1 = H(n_H-delta(ind-1):end-(n_H-delta(ind-1))+1,:);
            end
            %% Rxx e y
            v = randn(size(H1,1),1)/C; % ruido/SNR
            
            Rss = xcorr(S_bloco); % sinal
            Rss = toeplitz(Rss(length(S_bloco):end));
            
            Rvv= xcorr(v); % ruido
            Rvv = toeplitz(Rvv(size(H1,1):end));
            
            y = H1*S_bloco + v;
            %% Estimação de G e B
            
            aux = H1'/Rvv;
            aux2 = inv(Rss) + aux*H1;
            P = inv(aux2);
            [L,D] = ldl(P);
            
            G = D*L'*H1'/(Rvv);
            B = inv(L)-eye(n_amostras);
            %% Estimação
            for ind_bloco = 1:10
                
                % Sinal
                X = randn(n_amostras,1);
                S_bloco = sign(X);
                v = randn(size(H1,1),1)/C; %ruido
                y = H1*S_bloco + v;
                
                % G
                X_est = G*y;
                S_est = sign(X_est);
                e_til_G(:,ind_bloco) = (S_bloco ~= S_est);
                
                % B
                X_est2 = zeros(n_amostras,1);
                S_est2 = zeros(n_amostras,1);
                S_est2(1) = S_est(1);
                
                ind_amostras = 1; %B = 0
                X_est2(ind_amostras) = X_est(ind_amostras)-B(ind_amostras,:)*S_est2;
                S_est2(ind_amostras) = sign(X_est2(ind_amostras));
                
                for ind_amostras = 2:n_amostras
                    X_est2(ind_amostras) = X_est(ind_amostras)-B(ind_amostras,:)*S_est2;
                    S_est2(ind_amostras) = sign(X_est2(ind_amostras));
                end
                
                e_til_B(:,ind_bloco) = (S_bloco ~= S_est2);
                
            end
            % Prob. Erro ensemble
            prob_erro_G_en(n_loop) = sum(sum(e_til_G))/(n_amostras*ind_bloco);
            prob_erro_B_en(n_loop) = sum(sum(e_til_B))/(n_amostras*ind_bloco);
        end
        
        % Prob. Erro SNR
        prob_erro_G_SNR(ind,ind_SNR) = sum(prob_erro_G_en)/n_loop;
        prob_erro_B_SNR(ind,ind_SNR) = sum(prob_erro_B_en)/n_loop;
        
        %         disp(['Pot.Sinal = ' num2str(ind_SNR)])
        %         disp(['Prob. Erro (G) = ' num2str(prob_erro_G_SNR(ind_SNR)*100, '%f') '%'])
        %         disp(['Prob. Erro (B) = ' num2str(prob_erro_B_SNR(ind_SNR)*100, '%f') '%'])
    end
    %%
    
    
%     x = 1:ind_SNR;
%     %berTheory = (1/2)*erfc(sqrt(x));
%     figure(ind)
%     hold on
%     semilogy(x,prob_erro_G_SNR(ind,:),'r')
%     semilogy(x,prob_erro_B_SNR(ind,:),'b')
%     % semilogy(x,berTheory,'k')
%     legend('Estimação G','Estimação G+B')
%     xlabel('SNR (dB)')
%     ylabel('Prob. Erro (log)')
%     if ind==1
%         title('H')
%     else
%         title(['H - delta = ' num2str(delta(ind-1))])
%     end
%     grid
%     xlim([1 ind_SNR])
    
end
%%
% for ind = 1:length(delta)+1
%     x = 1:ind_SNR;
%     %berTheory = (1/2)*erfc(sqrt(x));
%     figure(ind)
%     hold on
%     semilogy(x,prob_erro_G_SNR(ind,:),'r')
%     semilogy(x,prob_erro_B_SNR(ind,:),'b')
%     % semilogy(x,berTheory,'k')
%     legend('Estimação G','Estimação G+B')
%     xlabel('SNR (dB)')
%     ylabel('Prob. Erro (log)')
%     if ind==1
%         title('H')
%     else
%         title(['H - delta = ' num2str(delta(ind-1))])
%     end
%     grid
%     xlim([1 ind_SNR])
% end
%%
x = 1:ind_SNR;
figure(ind+1)
semilogy(x,prob_erro_G_SNR)
legend('H', 'delta = 4', 'delta = 5', 'delta = 6', 'delta = 7', 'delta = 8')
xlim([1 ind_SNR])
xlabel('SNR (dB)')
ylabel('Prob. Erro (log)')
title(['Estimação G (' num2str(n_loop) ')'])

figure(ind+2)
semilogy(x,prob_erro_B_SNR)
legend('H', 'delta = 4', 'delta = 5', 'delta = 6', 'delta = 7', 'delta = 8')
xlim([1 ind_SNR])
xlabel('SNR (dB)')
ylabel('Prob. Erro (log)')
title(['Estimação G+B (' num2str(n_loop) ')'])