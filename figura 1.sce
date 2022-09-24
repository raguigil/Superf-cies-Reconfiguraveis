//Emil Björnson, Henk Wymeersch, Bho Matthiesen, Petar Popovski, Luca
//Sanguinetti, e Elisabeth de Carvalho “Reconfigurable IntelligentSurfaces: uma perspectiva de processamento de sinal com aplicativos sem fio”,
//IEEE Signal Processing, a ser publicada em março de 2021.
//IEEE Signal Processing, a ser publicada em março de 2022.
//
//Download do artigo: https://arxiv.org/pdf/2102.00742.pdf
//
//Esta é a versão 1.0 (Última edição: 2022-01-01

close all;
clear;
function vn = refcoefficient(omega,Cn)
// Calcula o coeficiente de reflexão usando o circuito de (3) em
//"Superfície refletora inteligente: modelo prático de mudança de fase e
//Beamforming Optimization" por Samith Abeywickrama, Rui Zhang, Chau Yuen.

L1 = 2.5e-9; //Indutância da camada inferior
L2 = 0.7e-9; //Indutância da camada superior
Rn = 1; //Resistência efetiva
Z0 = 377; //Impedância do espaço livre

//Calcula a impedância da superfície
Zn = 1*%i*omega*L1*(1*%i*omega*L2+1./(1*%i*omega*Cn)+Rn)./(1*%i*omega*L1+(1*%i*omega*L2+1./(1*%i*omega*Cn)+Rn));

 //Calcula o coeficiente de reflexão
vn = (Zn-Z0)./(Zn+Z0);

end

//Intervalo de valores de frequência
fRange = linspace(2,4,100);
//Valor de capacitância variável
Cn = [0.79 0.88 0.96 2.2]*1e-12;
//Preparar para armazenar resultados
vn = zeros(length(fRange),length(Cn));
for k = 1:length(fRange)
    //Calcular frequência angular
    omega = 2*%pi*fRange(k)*1e9;
    for m = 1:length(Cn)
         //Calcular o coeficiente de reflexão usando o circuito de (3) em
        // "Superfície Refletora Inteligente: Modelo Prático de Deslocamento de Fase e
        // Beamforming Optimization" por Samith Abeywickrama, Rui Zhang, Chau Yuen.
        vn(k,m) =  refcoefficient(omega,Cn(m));
    end
end
function [ydb]=pow2db(y)
    ydb=10*log10(y)
endfunction
function [y]=db2pow(ydb)
    y=10^(ydb/10)
endfunction
//set(groot,'defaultAxesTickLabelInterpreter','latex');
       //Plote os resultados da simulação
subplot(2,1,1)
mtlb_hold on; mtlb_box on; mtlb_grid on;
plot(fRange,atan(imag(vn(:,1)),real(vn(:,1))),'y','LineWidth',3);
plot(fRange,atan(imag(vn(:,2)),real(vn(:,2))),'g--','LineWidth',3);
plot(fRange,atan(imag(vn(:,3)),real(vn(:,3))),'b-.','LineWidth',3);
plot(fRange,atan(imag(vn(:,4)),real(vn(:,4))),'b:','LineWidth',3);
//set(gca,'fontsize',16);
xlabel('$Frequencia (f) [GHz]$');
ylabel('$Resposta \ de \ fase [rad]$');
legend("$0,79 \ pF$","$0,88pF$","$0,96 \ pF$","$2,2 \ pF$");
//ylim([-%pi %pi]);
//yticks(-pi:pi/2:pi);
//yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
subplot(2,1,2)
mtlb_hold on; mtlb_box on;mtlb_grid on;
plot(fRange,pow2db(abs(vn(:,1))),'y','LineWidth',3);
plot(fRange,pow2db(abs(vn(:,2))),'g--','LineWidth',3);
plot(fRange,pow2db(abs(vn(:,3))),'b-.','LineWidth',3);
plot(fRange,pow2db(abs(vn(:,4))),'b:','LineWidth',3);
//set(gca,'fontsize',16);
xlabel('$Frequencia [GHz]$');
ylabel('$Resposta \ de \ amplitude\ [dB]$');
legend("$0,79 pF$","$0,88pF$","$0,96pF$","$2,2pF$",4)
//ylim([-3 0]);
