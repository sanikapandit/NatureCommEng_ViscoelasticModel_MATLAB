function J = KelvinLSQ(p,Nt,t,tp,up,Fo,FData2)

k1 = p(1);
T1 = p(2);
T2 = p(3);

FData=FData2';

for i=1:Nt
    if t(i)<tp
        F(i,1) = (k1*T2/T1)*up*t(i)/tp +Fo*exp(-t(i)/T1) - (k1*(T2-T1)/(T1^2))*(up/tp*(T1*t(i)-T1^2+(T1^2)*exp(-t(i)/T1)));
    else
        F(i,1) = (k1*T2/T1)*up*t(i)/tp +k1*T2*(up-up*t(i)/tp)/T1 +Fo*exp(-t(i)/T1) - ...
            (k1*(T2-T1)/(T1^2))* (up/tp*(T1*t(i)-T1^2+(T1^2)*exp(-t(i)/T1)) + (up-up*t(i)/tp)*T1 +up/tp*T1*T1 - ...
            up/tp*T1*T1*exp((tp-t(i))/T1));
    end
end


residualVec = 1/sqrt(Nt)*(F - FData);

J = residualVec'*residualVec + (T1>T2);