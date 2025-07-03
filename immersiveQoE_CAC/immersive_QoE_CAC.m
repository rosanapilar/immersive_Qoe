%%
% Running the program, get four values:
% 1. Uniform rendering power allocation strategy
% 2. Random rendering power allocation strategy
% 3. Optimial power Allocation (GT) - bound
% 4. Optimial Allocation base on our Predictions

clc;clear;
A1 = readmatrix("gd50.txt"); % The ground truth
A2 = readmatrix("predall50.txt"); % Prediction results
A3 = readmatrix("experiment50.txt");% Randomly generated sparse interactions
% A3 = readmatrix('experiment2.txt');% Randomly generated sparse interactions
A3 = A3+1; % Start with 1
T = readmatrix("parametros.csv");

FINAL=zeros(50,6);
Attention = -1.*ones(50,59);
%% Generate attention matrix

%Rosana - parte do RDk e EUk
Rd = [];  %vetorde taxa de downlink
Eu = [];  %vetor de taxa de BEP de uplink 
Lat= [];
cone = [];  %vetor do coeficiente de conexao (KPIs objetivos)
coneo = [];  %vetor do coeficiente de conexao (KPIs objetivos)


Rd = zeros(1,50);
Eu = zeros(1,50);
Lat = zeros(1,50);
Lato = zeros(1,50);
cone = zeros(1,50);
coneo = zeros(1,50);
render = zeros(1,50);


T = readmatrix("parametros.csv");

% Extrai os vetores das colunas da tabela
Rd = T(:,1)';    %
Eu = T(:,2)';    % coluna 2, transposta para vetor linha
Lat = T(:,3)';
cone = T(:,4)';




for u = 1:50
usernum = u;% change to try different users (1~30)

Atemp = [];
for k = 1: length(A3(usernum,:))
    if A3(usernum,k)>=0
        Atemp(k) = A3(usernum,k);
    end
end

uoal = [];
uoalpre = [];
cixu = [];

for k = 1:length(Atemp)
    uoal(k) = A1(usernum,Atemp(k)); % User attention to different objects (GT)
    uoalpre(k) = A2(usernum,Atemp(k)); % Predicted user attention for different objects
    Attention(usernum,k) = uoal(k);
end

numO = length(Atemp); % Total number of objects in one virtual tour

%% Initialize rendering power
% PthR = 15; % Minimum rendering power per object
% PkR = 1000; % The total rendering power of one virtual metaverse service assigned to user k

PthR = 15;
PkR = numO*20;

if PthR.*length(Atemp)>PkR
    disp('not availiable');
    finish
end

PnkR = zeros(1,length(uoal)); % Initialize the power assigned to each object



%% Optimial Allocation Predictions  %%ciente da atenção
PnkR = zeros(1,length(uoal));
uxing = sum(uoalpre)/PkR;
PnkR = uoalpre./uxing;
j = 1;
t1 = [];t2 = [];
while min(PnkR)<PthR 
    [a,b] = min(PnkR); 
    t1(j) = b;
    t2(j) = uoalpre(b); 
    uxing = (sum(uoalpre)-sum(t2))/(PkR - PthR*j); 
    PnkR = uoalpre./uxing; 
    for q = 1:j
    PnkR(t1(q)) = PthR;
    end
    j = j+1;
   %sum(PnkR);
end
FINAL(u,1) = sum(uoal.*log(PnkR./PthR));
%FINAL(u,2) = sum(uoal.*log(PnkR./PthR));
FINAL(u,2) = (Rd(usernum) * (1-Eu(usernum)) /(1-Lat(usernum)))* sum(uoal.*log(PnkR./PthR)); -- 24/10/24
%FINAL(u,2) = cone(usernum)*0.89* sum(uoal.*log(PnkR./PthR)); %-- 24/10/24
%FINAL(u,3) = coneo(usernum); %-- 24/10/24
FINAL(u,3) = Rd(u) * sum(uoal .* log(PnkR ./ PthR));
FINAL(u,4) = (1 - Eu(u)) * sum(uoal .* log(PnkR ./ PthR));
FINAL(u,5) = (1 / (1 - Lat(u))) * sum(uoal .* log(PnkR ./ PthR));
% sum(PnkR)


%% Optimial Allocation GT
% PnkR = zeros(1,length(uoal));
% uxing = sum(uoal)/PkR;
% PnkR = uoal./uxing;
% j = 1;
% t1 = [];t2 = [];
% while min(PnkR)<PthR 
% % When the condition that the minimum rendering power 
% % must be greater than PthR is not satisfied
%     [a,b] = min(PnkR); 
%     % a records the minimum renderning power,
%     % b records the corresponding position
%     t1(j) = b;
%     t2(j) = uoal(b);
%     uxing = (sum(uoal)-sum(t2))/(PkR - PthR*j); %Solve for the new u*
%     PnkR = uoal./uxing; 
%     for q = 1:j
%     PnkR(t1(q)) = PthR;
%     end
%     j = j+1;
% %     sum(PnkR)
% end
% FINAL(u,4) = sum(uoal.*log(PnkR./PthR));
% sum(PnkR)

end
%% Plot
figure
%FINAL=[FINAL;mean(FINAL)];
wzi = 14;
bar(FINAL);grid on;
bar(FINAL(:,1:2)); % plota as 2 primeiras coluna
%axis([0 32 0 75])
%axis([0 32 0 40])
%axis([0 32 0 2000])
xlabel('Usuários')
ylabel('QoE')
%legend('QoE Optimal (H.Du)','QoE-CAC', 'Qoe com Latencia Otimizada')
legend('QoE Optimal (H.Du)','QoE-CAC')
set(gca,'fontname','Times New Roman','FontSize',wzi);
 Diff1 = (FINAL(:,1)-FINAL(:,2))./FINAL(:,2);%
% Diff2 = (FINAL(:,4)-FINAL(:,3))./FINAL(:,3);% 
% 
% x = mean(Diff1')
% mean(Diff2')
% disp([Diff1'.*100])
% Diff2'.*100
% 

% 
% plotando separado

figure;
bar(FINAL(:,1:5)); % plota as 5 primeiras colunas
grid on;
xlabel('Usuários');
ylabel('QoE');
%legend('QoE Optimal (H.Du)','QoE-CAC','QoE-Rd', 'QoE-Eu', 'QoE-Lat');
legend('QoE Optimal (H.Du)','QoE-CAC','QoE-Rd');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);




%plotando parametros

figure;
i = 1:size(T, 1);

% Separando os dados

% Criando o gráfico

plot(i, Rd, '-o', 'Color', 'b',  'DisplayName', 'R_d(u)');
hold on;
plot(i, Eu, '-x', 'Color', 'r', 'DisplayName', 'E_u(u)');
plot(i, Lat, '-s', 'Color', 'g', 'DisplayName', 'Lat(u)');
% %plot(i, cone, '-d',  'DisplayName', 'cone(u)');


% Personalização do gráfico
xlabel('Usuários');
ylabel('Valor');
title('Valores dos parâmetros de Recursos de Rede');
legend('Location', 'best');
hold off;