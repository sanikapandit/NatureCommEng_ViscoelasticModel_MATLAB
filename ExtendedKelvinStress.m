function F = ExtendedKelvinStress(p,Nt,t,tp,up,Fo)

%lack of use of Fo? is this an issue?

k1 = p(1);
k2 = p(2);
k3 = p(3);
T1 = p(4);
T2 = p(5);

for i=1:Nt
    if t(i)<tp
        F(i,1) = up/tp*(k1 + k2 - k1*exp(-t(i)/T1) - k2*exp(-t(i)/T2) +k3*t(i));
        %F(i,1) = (k1*T2/T1)*up*t(i)/tp +Fo*exp(-t(i)/T1) - (k1*(T2-T1)/(T1^2))*(up/tp*(T1*t(i)-T1^2+(T1^2)*exp(-t(i)/T1)));
    else
        F(i,1) = k1*up/tp*(exp(-(t(i)-tp)/T1)-exp(-t(i)/T1)) + ...
            k2*up/tp*(exp(-(t(i)-tp)/T2)-exp(-t(i)/T2)) + k3*up;
        %F(i,1) = (k1*T2/T1)*up*t(i)/tp +k1*T2*(up-up*t(i)/tp)/T1 +Fo*exp(-t(i)/T1) - ...
            %(k1*(T2-T1)/(T1^2))* (up/tp*(T1*t(i)-T1^2+(T1^2)*exp(-t(i)/T1)) + (up-up*t(i)/tp)*T1 +up/tp*T1*T1 - ...
            % up/tp*T1*T1*exp((tp-t(i))/T1));
    end
end