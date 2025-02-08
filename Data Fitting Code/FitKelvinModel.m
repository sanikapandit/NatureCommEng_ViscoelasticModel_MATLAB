%Preliminary Math Modeling FBN Gels
%Author: Nina Moiseiwitsch
close all

%Plot Sample Data
figure(1)
inputcurve = Overnightrampholddata;
plot(inputcurve{:,1},inputcurve{:,3})
xlabel('Time (s)');
ylabel('Stress (N)');
axis([0 300 0 0.7])
hold on

%Prompt user to enter experimentally determined values
% prompt = {'Enter peak stress (Fp)','Enter equilibrium stress (F as t --> infinity)','Enter F(0)', 'Enter time at peak (tp)', 'Enter u(0)'};
% dlgtitle = 'Input Experimental Values';
% experimental_inputs = inputdlg(prompt,dlgtitle)

%Define variables based on user input
%Fp = str2num(experimental_inputs{1});
Fp = 0.54028
%Finf = str2num(experimental_inputs{2});
Finf = 0.36
%Fo = str2num(experimental_inputs{3});
Fo = 0
%tp = str2num(experimental_inputs{4});
tp = 30
%up = str2num(experimental_inputs{5});
up = 0.0050054

%define variables 
k1 = Finf/up;

T1 = 10;

F=zeros(1,100);

% define conditional f(x)
T2 = T1 + (1-exp(-tp/T1))^(-1) * (tp/(Finf))*(Fp-Finf-Fo*exp(-tp/T1));
tInit=0;
%tFinal=500;
%dt=tFinal/Nt;
%t=dt:dt:tFinal;
t=inputcurve{:,1}';
Nt=size(t,2);
FData=inputcurve{:,3}';

pInit(1)=k1;
pInit(2)=T1;
pInit(3)=T2;

options=optimset('Display','iter','TolFun',1e-9,'MaxFun',1e8,'TolX',1e-9,'MaxIter',5e8);
[pOpt,lsqVal,exitflag,output] = fminsearch(@KelvinLSQ,pInit,options,Nt,t,tp,up,Fo,FData);

k1Opt=pOpt(1);
T1Opt=pOpt(2);
T2Opt=pOpt(3);

fprintf('k1Opt=%8.4f T1Opt=%8.4f T2Opt=%8.4f Cost=%8.5e\n',k1Opt,T1Opt,T2Opt,lsqVal)

if T2Opt<=T1Opt
    error('Tau constraint violated!')
end

FOpt = KelvinStress(pOpt,Nt,t,tp,up,Fo);

plot(t,FOpt,'r')
axis([0 t(end) 0 0.7])

residual = FOpt - FData';

figure(2)
plot(t,residual)
xlabel('Time (s)');
ylabel('Residual');




