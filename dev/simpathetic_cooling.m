f = @(n) n+1+2/n;
omega_ratio = 1;%sqrt(4/3);



T0_on_TF = 0.2;
eta = 6;

c = 4/5*pi^2*(omega_ratio)^3;
T1_on_TF = linspace(0,T0_on_TF);

dNb_on_Nf = -(T1_on_TF+c/3*T1_on_TF.^3-T0_on_TF-c/3*T0_on_TF.^3).*pi^2./f(eta);

figure(24)
hold on
plot(dNb_on_Nf,T1_on_TF)