clear all;
load('Cubo1216vpaper.mat');

Inicio = 1;
reduccion = 1;
ReducSpam = 200;
Mmax = round((nt-Inicio)/reduccion);

Datacruda=Ezt(Inicio:end);
Transmitida((Mmax))=0;
for I=1:round(Mmax);
Transmitida(I)=Datacruda(I*reduccion);
end

T = reduccion*dt;
Fs = 1/T;


L = Mmax;
%t = (0:L-1)*T;
%S = (0.7)*sin(2*pi*50*t);
S = Transmitida;
Y = fft(S,round(nt*ReducSpam));
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1)=2*P1(2:end-1);
f= Fs*(0:(L/2))/(L*ReducSpam);
%plot(f,P1,'green');
%subplot(4,1,4); plot(f,P1,'green');figure(gcf);
hold;

Datacruda2=Ezi0(Inicio:end);
Incidente(round(Mmax))=0;
for I=1:round(Mmax);
Incidente(I)=Datacruda2(I*reduccion);
end


%L = 1500;
%L = Mmax;
%t = (0:L-1)*T;
%S = (0.7)*sin(2*pi*50*t);
S2 = Incidente;
Y2 = fft(S2,nt*ReducSpam);
P22 = abs(Y2/L);
P12 = P22(1:L/2+1);
P12(2:end-1)=2*P12(2:end-1);
%f= 0.95*Fs*(0:(L/2))/(L*ReducSpam);
f= Fs*(0:(L/2))/(L*ReducSpam);
%subplot(4,1,4); plot(f,P12,'red');figure(gcf);

for I=1:numel(f);
%Eficiencia(I) = (10/9.4)*(1-((eta0-etaH2O)/(eta0+etaH2O))*((eta0-etaH2O)/(eta0+etaH2O)))*(P1(I)*P1(I))/(P12(I)*P12(I));
Eficiencia(I) = (1-((eta0-etaH2O)/(eta0+etaH2O))*((eta0-etaH2O)/(eta0+etaH2O)))*(P1(I)*P1(I))/(P12(I)*P12(I));

end

%KL = ((eta0-etaH2O)/(eta0+etaH2O))*((eta0-etaH2O)/(eta0+etaH2O));

plot(f,Eficiencia,'LineWidth',2,'Color','black');

clear all;

I=1;
N=1000;
Zaire = 377;
epsilonragua = 81;
Zagua = 377/(epsilonragua)^(1/2);
LC(N:2)=0;
c = 3e8;
f0 = 2450e6;
deltaf = 8e6;
v = c/(epsilonragua)^(1/2);
lambdav = v/f0;
%plot(LC(:,1),LC(:,2));
axis([100e6 4000e6 0 1])
hold on;
Z0transf = (Zagua*Zaire)^(1/2);
d1 = 0.25*lambdav;

I=1;

while (I < (N));
    

f=f0-(deltaf*(N-2*I));
w=6.283*f;
beta = w/v;

GammaStub1 = ((Zagua-Z0transf)/(Zagua+Z0transf));
Zin1 = Z0transf*(1+GammaStub1*exp(-1i*2*beta*d1))/(1-GammaStub1*exp(-1i*2*beta*d1));

ModGamma = abs((Zin1-Zaire)/(Zin1+Zaire));

ROE = (1+ModGamma)/(1-ModGamma);
Eficiencia = 1 - ModGamma*ModGamma;


TransfCuartodeOnda(I,1)=f;
TransfCuartodeOnda(I,2)=Eficiencia;

I=I+1;
end

plot(TransfCuartodeOnda(:,1),TransfCuartodeOnda(:,2),'--','Color','black');
legend('FDTD numerical model','Proposed analytical model');
