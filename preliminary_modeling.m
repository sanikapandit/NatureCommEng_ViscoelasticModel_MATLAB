%Preliminary Math Modeling FBN Gels
%Author: Nina Moiseiwitsch
close all

%Plot Sample Data
figure(1)
inputcurve = SampleDataCurveJul21;
plot(inputcurve{:,"Time"},inputcurve{:,"Stress"})
xlabel('Time (s)');
ylabel('Stress (kN)');
axis([0 70 4.09 4.12])


%Prompt user to enter experimentally determined values
% prompt = {'Enter peak stress (Fp)','Enter equilibrium stress (F as t --> infinity)','Enter F(0)', 'Enter time at peak (tp)', 'Enter u(0)'};
% dlgtitle = 'Input Experimental Values';
% experimental_inputs = inputdlg(prompt,dlgtitle)

%Define variables based on user input
%Fp = str2num(experimental_inputs{1});
Fp = 4.112
%Finf = str2num(experimental_inputs{2});
Finf = 4.09*0.85;
%Fo = str2num(experimental_inputs{3});
Fo=4.0906
%tp = str2num(experimental_inputs{4});
tp = 50
%uo = str2num(experimental_inputs{5});
uo = 0.00799997

%define variables 
k1 = uo/Finf;

T1 = 50;

F=zeros(1,100);

% define conditional f(x)
T2 = T1 + (1-exp(-tp/T1))^(-1) * (tp/Finf)*(Fp-Finf-Fo*exp(-tp/T1));
if T2 <= T1
    error('Tau constraint violated!')
end

tInit=0;
tFinal=70.0;
Nt=140;
dt=tFinal/Nt;
t=dt:dt:tFinal;
for i=1:Nt
    if t(i)<tp
        F(1,i) = (k1*T2/T1)*uo*t(i)/tp +Fo*exp(-t(i)/T1) - ...
            (k1*(T2-T1)/(T1^2))*(uo/tp*(T1*t(i)-T1^2+(T1^2)*exp(-t(i)/T1)));
    else
        F(1,i) = (k1*T2/T1)*uo*t(i)/tp +k1*T2*(uo-uo*t(i)/tp)/T1 +Fo*exp(-t(i)/T1) - ...
            (k1*(T2-T1)/(T1^2))* (uo/tp*(T1*t(i)-T1^2+(T1^2)*exp(-t(i)/T1)) + ...
            (uo-uo*t(i)/tp)*T1 +uo/tp*T1*T1 - ...
            uo/tp*T1*T1*exp((tp-t(i))/T1));
    end
end

figure(2)
plot(t,F)
xlabel('time (s)')
ylabel('Stress (N)')

