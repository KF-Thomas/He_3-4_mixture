f = @(n) n+1+2/n;
omega_ratio = sqrt(4/3);%5;%



T0_on_TF = 0.2;
eta = 6;

c = 4/5*pi^2*(omega_ratio)^3;
T1_on_TF = linspace(0,T0_on_TF);

dNb_on_Nf = -(T1_on_TF+c/3*T1_on_TF.^3-T0_on_TF-c/3*T0_on_TF.^3).*pi^2./f(eta);

figure(24)
hold on
plot(dNb_on_Nf,T1_on_TF)

func = @(b,x) b(1).*x(:,2).*(x(:,1)+c/3.*x(:,1).^3-b(2))+b(3);

const_ratio = (zeta(3)*6)^(-1/3);
func = @(b,x) b(1).*x(:,2).*(x(:,1)+c/3.*x(:,1).^3)+b(1).*(const_ratio.*x(:,2).^(2/3).*b(2).^(1/3)+c/3.*const_ratio^3.*b(2))+b(2);

func2 = @(b,x) b(1).*(x).^b(2);
% func3 = @(b,x) b(1).*(x+c/3.*x.^3-b(2).^(1/3)*(0.448)-c/3*b(2)^3*(0.448)^3)+b(2);
func3 = @(b,x) (pi^2/3.*x-(pi^2/3*x./abs(b(1))).^b(2))./(1/b(2)-1);
T_ratio = T4.'./TF_c.';
T_unc  = T4_unc./TF_c;
T_ratio(abs(T_ratio)>10) = nan;
T_unc(abs(T_ratio)>10) = nan;
weights = 1./T_unc.^4;
fitobject=fitnlm([T_ratio,N3.'],N4./qe,func,[1,1e4],...
        'options',opt);
fitobject3=fitnlm(T_ratio.',N4./qe./N3,func3,[-0.644,3],'Weights',weights,...
        'options',opt)
fitobject3=fitnlm(T_ratio.',N4./qe./N3,func3,[-0.644,3],...
        'options',opt)
fitobject3=fitnlm(T_ratio.',N4./qe./N3,func3,[-0.644,3],'ErrorModel','combined')
fitobject3=fitnlm(T_ratio.',N4./qe./N3,func3,[-0.644,3],'ErrorModel','proportional')
func3 = @(b,x) real(log((pi^2/3.*x-(pi^2/3*x./abs(b(1))).^b(2))./(1/b(2)-1)));
fitobject3=fitnlm(T_ratio.',log(N4./qe./N3),func3,[-0.644,3])

fitobject2=fitnlm(N4./qe./N3,T_ratio.',func2,[1,1/3],...
        'options',opt);