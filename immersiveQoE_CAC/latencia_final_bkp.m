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

FINAL=zeros(50,3);
Attention = -1.*ones(50,59);
%% Generate attention matrix

%Rosana - parte do RDk e EUk
Rd = [];  %vetorde taxa de downlink
Eu = [];  %vetor de taxa de BEP de uplink 
Lat= [];
cone = [];  %vetor do coeficiente de conexao (KPIs objetivos)
coneo = [];  %vetor do coeficiente de conexao (KPIs objetivos)




Rd = zeros(1,30);
Eu = zeros(1,30);
Lat = zeros(1,30);
Lato = zeros(1,30);
cone = zeros(1,30);
coneo = zeros(1,30);
render = zeros(1,30);

for k = 1:50
     Rd(k) = randi([20, 42]);
     Lat(k) = 0.04 + (0.05 - 0.04) * rand; % Gera valores entre 40 ms e 50 ms
     Lato(k) = 0.02 + (0.04 - 0.02) * rand; % Preserva Lato(k) se necessário
  
    Eu(k)= 10^-8+ (10^-2 - 10^-8)*rand;
    
    %Eu(k)= rand(1)/10;  
end
for i = 1:50
 %normatizando os valores de Rd 
   Rd(i) =  (Rd(i)  - min (Rd)) / (max(Rd)  - min (Rd));  
   Eu(i) =  (Eu(i)  - min (Eu)) / (max(Eu)  - min (Eu)); 
   Lat(i) =  (Lat(i)  - min (Lat)) / (max(Lat)  - min (Lat));
   Lato(i) =  (Lato(i)  - min (Lato)) / (max(Lato)  - min (Lato));
    



  % cone(i)= Rd(i) * (1-Eu(i)) * Lat(i);  %% com latencia
   %cone(i)= Rd(i) * (1-Eu(i)) /(1 +Lat(i));  %% com latencia 24/1024
   %%% teste 11/12/24
   cone(i)= Rd(i) * (1-Eu(i)) /(1-Lat(i));  %% com latencia
   coneo(i)= (1-Lato(i));  %% com latencia

%%% teste somente com latencia
%    cone(i)= 1 /(1-Lat(i));  %% com latencia
%    coneo(i)= 1/(1-Lato(i));  %% com latencia
%  cone(i)= 1;
end
%%Rosana
fid = fopen('dados_Rd.csv', 'w');
fprintf(fid, 'Rd(i)\n');
for i = 1:50
    fprintf(fid, '%f\n', Rd(i));
end
fclose(fid);


% Escrevendo os dados linha por linha em um arquivo CSV
fid = fopen('parametros.csv', 'w');
fprintf(fid, 'Rd(i), Eu(i), Lat(i), cone(i)\n');
for i = 1:50
    fprintf(fid, '%f, %f, %f, %f\n', Rd(i), Eu(i),  Lat(i), cone(i));
end
fclose(fid);





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
%FINAL(u,2) = cone(usernum)* sum(uoal.*log(PnkR./PthR)); -- 24/10/24
FINAL(u,2) = cone(usernum)*0.89* sum(uoal.*log(PnkR./PthR)); %-- 24/10/24
%FINAL(u,3) = coneo(usernum); %-- 24/10/24

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

