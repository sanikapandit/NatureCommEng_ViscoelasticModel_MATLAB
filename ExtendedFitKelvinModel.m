%Preliminary Math Modeling FBN Gels
%Author: Nina Moiseiwitsch
close all

%Plot Sample Data
figure(1)
inputcurve = Overnightrampholddata;
plot(inputcurve{:,1},inputcurve{:,4})
xlabel('Time (s)');
ylabel('Stress (N)');
%axis([0 300 0 0.05])

%Prompt user to enter experimentally determined values
% prompt = {'Enter peak stress (Fp)','Enter equilibrium stress (F as t --> infinity)','Enter F(0)', 'Enter time at peak (tp)', 'Enter u(0)'};
% dlgtitle = 'Input Experimental Values';
% experimental_inputs = inputdlg(prompt,dlgtitle)

%Define variables based on user input
%Fp = str2num(experimental_inputs{1});
Fp = 0.55
%Finf = str2num(experimental_inputs{2});
Finf = 0.45
%Fo = str2num(experimental_inputs{3});
Fo = 0
%tp = str2num(experimental_inputs{4});
tp = 30
%up = str2num(experimental_inputs{5});
up = 0.005

%define variables 
k1 = Finf/up;

T1 = 100;

F=zeros(1,100);

% define conditional f(x)
T2 = T1 + (1-exp(-tp/T1))^(-1) * (tp/(Finf))*(Fp-Finf-Fo*exp(-tp/T1));
tInit=0;
%tFinal=500;
%dt=tFinal/Nt;
%t=dt:dt:tFinal;
t=inputcurve{:,1}';
Nt=size(t,2);
FData=inputcurve{:,4}';
FData = FData - min(FData)

%How to add in additional variables is this section??
pInit(1)=k1;
pInit(2)=k1;
pInit(3)=k1;
pInit(4)=T1;
pInit(5)=T2;

options=optimset('Display','iter','TolFun',1e-9,'MaxFun',1e8,'TolX',1e-9,'MaxIter',5e8);
[pOpt,lsqVal,exitflag,output] = fminsearch(@ExtendedKelvinLSQ,pInit,options,Nt,t,tp,up,Fo,FData);

%How to add in additional variables is this section??
k1Opt=pOpt(1);
k2Opt=pOpt(2);
k3Opt=pOpt(3);
T1Opt=pOpt(4);
T2Opt=pOpt(5);

fprintf('k1Opt=%8.4f k2Opt=%8.4f k3Opt=%8.4f T1Opt=%8.4f T2Opt=%8.4f Cost=%8.5e\n',k1Opt,k2Opt,k3Opt,T1Opt,T2Opt,lsqVal)

% if T2Opt<=T1Opt
%     error('Tau constraint violated!')
% end

FOpt = ExtendedKelvinStress(pOpt,Nt,t,tp,up,Fo);

Residual = FOpt-FData';

figure(2)
plot(t,FData,'b')
%axis([0 200 0 0.05])
hold on
plot(t,FOpt,'r')

figure(3)
plot(t,Residual,'r')
