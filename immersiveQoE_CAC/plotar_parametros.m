% Dados fornecidos
%Rd(i), Eu	(i), Lat(i)
% inputs: downlink, uplink, latencia, renderização
% output: qoe

dados = [
0.545455,	 0.012196,	 0.20180,	 0.007253;
0.855263,	 0.006703,	 0.21070,	 0.005425;
0.581726,	 0.042221,	 0.220120,	 0.011806;
0.590562,	 0.044881,	 0.121290,	 0.011632;
0.480899,	 0.017934,	 0.194180,	 0.010439;
0.481726,	 0.032209,	 0.227010,	 0.016164;
0.901281,	 0.013690,	 0.195220,	 0.007497;
0.827242,	 0.011033,	 0.246500,	 0.017973;
0.851922,	 0.033111,	 0.205560,	 0.012950;
0.481726,	 0.020300,	 0.167090,	 0.021721;
0.802562,	 0.009576,	 0.152980,	 0.033578;
0.753203,	 0.017311,	 0.166090,	 0.023380;
0.851922,	 0.040418,	 0.197300,	 0.026188;
0.950641,	 0.011005,	 0.245250,	 0.042412;
0.925961,	 0.032682,	 0.116040,	 0.021168;
0.654484,	 0.035410,	 0.173420,	 0.012688;
0.975320,	 0.014416,	 0.117980,	 0.041078;
0.901281,	 0.033760,	 0.146010,	 0.018375;
1.000000,	 0.027426,	 0.159840,	 0.021483;
0.555765,	 0.025862,	 0.194950,	 0.035633;
0.901281,	 0.000000,	 0.169700,	 0.009823;
0.580445,	 0.033587,	 0.193510,	 0.013731;
0.629804,	 0.032139,	 0.116030,	 0.021163;
0.728523,	 0.004559,	 0.184860,	 0.004914;
0.481726,	 0.007063,	 0.160260,	 0.042168;
0.728523,	 0.030585,	 0.115200,	 0.026688;
1.000000,	 0.038894,	 0.241440,	 0.028567;
1.000000,	 0.018423,	 0.180040,	 0.005994;
0.970166,	 0.018816,	 0.214250,	 0.042330;
1.000000,	 0.029594,	 0.213400,	 0.028356];

% Índices para o eixo X
i = 1:size(dados, 1);

% Separando os dados
Rd = dados(:,1);
Eu = dados(:,2);
Lat = dados(:,3);
cone = dados(:,4);

% Criando o gráfico
figure;
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