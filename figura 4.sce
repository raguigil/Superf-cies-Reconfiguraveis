
//Este script Scilab baseado em Matlab gera partes da Figura 4 no artigo:

//Emil Björnson, Henk Wymeersch, Bho Matthiesen, Petar Popovski, Luca
//Sanguinetti, e Elisabeth de Carvalho “Reconfigurable Intelligent
//Superfícies: uma perspectiva de processamento de sinais com aplicativos sem fio”,
//IEEE Signal Processing Magazine, a ser publicado em março de 2022.
//
//Download article: https://arxiv.org/pdf/2102.00742.pdf
//
//This is version 1.0 (Last edited: 2022-01-01)
//
clc;
clear;
function [y]=db2pow(ydb)
    y=10^(ydb/10)
endfunction
//Definir perdas de propagação por elemento
alpha = db2pow(-80); //Para o RIS
beta = db2pow(-60); //Do RIS
rho1 = db2pow(-80);//Opção 1 para caminho incontrolável
rho2 = db2pow(-110);//Opção 2 para caminho incontrolável
gamma = 1; //Perda em RIS

// //Largura de banda (Bandwidth) (Hz)
B = 1e6;

//Número do elemento RIS
N = 0:200;

//Definir SNR de transmissão (dependendo da largura de banda)
PBN0 = db2pow(100);

// Calcula SNRs usando (31) com configuração RIS ideal
SNR1 = PBN0*(sqrt(rho1)+N*sqrt(alpha*beta*gamma)).^2;
SNR2 = PBN0*(sqrt(rho2)+N*sqrt(alpha*beta*gamma)).^2;


//Número de realizações na simulação de Monte Carlo
numberOfRealizations = 10000;

//O sistema de coordenadas é selecionado para que o caminho incontrolável tenha zero
//fase, enquanto todos os outros caminhos têm fases uniformemente distribuídas
phaseShifts = rand(max(N),numberOfRealizations)*2*%pi;

//A fase de cada caminho é girada para ficar entre -pi/4 e pi/4, usando um
//RIS com 4 configurações
rotateShifts =pmodulo(phaseShifts,%pi/2)-%pi/4;

//Calcula os SNRs reais com um RIS com 4 configurações
SNR1b = PBN0*[rho1; mean(abs(sqrt(rho1)+cumsum(exp(1*%i*rotateShifts),1)*sqrt(alpha*beta*gamma)).^2,2)];
SNR2b = PBN0*[rho2; mean(abs(sqrt(rho2)+cumsum(exp(1*%i*rotateShifts),1)*sqrt(alpha*beta*gamma)).^2,2)];
// Plota os resultados da simulação
mtlb_hold on; mtlb_box on; mtlb_grid on;
plot(N,B*log2(1+SNR1)/1e6,'r-','LineWidth',3)
plot(N,B*log2(1+SNR1b)/1e6,'k--','LineWidth',3)
plot(N,B*log2(1+SNR2)/1e6,'b-.','LineWidth',3)
plot(N,B*log2(1+SNR2b)/1e6,'k:','LineWidth',3)
xlabel("$Número \ de \ elementos \ RIS \ [N]$");
ylabel("$Capacidade \ [Mbps]$");
legend("$\rho=-80 \ dB \ (ideal \ config)$","$\rho=-80 \ dB (4 \ configs)$","$\rho=-110 dB (ideal \ config)$","$\rho=-110 dB (4 \ configs)$",5);
