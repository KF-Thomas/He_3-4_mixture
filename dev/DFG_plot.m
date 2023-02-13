%% comparison of momentum density of Boltzmann and Fermi-Dirac distrabution
% in harmonic trap usng the semi classical approximation

% constats
hebec_constants
m = const.mhe*3/4;
hb = const.hb;
kb = const.kb;

omega = [60 600 600].*2*pi; %trapping frequency
omegar = 600*2*pi;
lambda = 1/10;
a0 = sqrt(hb./(m.*omega));
omega_bar = geomean(omega);
N = 1e5;%atom number
T = 200e-9;%temperature

EF = hb.*omega_bar.*(6*N).^(1/3);
KF = (2*m*EF./hb^2).^0.5;%fermi wavevector
TF = EF/kb;
mu = find_mu(omega_bar,T,N,EF);

H = @(k,x,y,z) 1/(2*m).*(hb.^2.*k.^2)+m/2.*(omega(1).^2.*x.^2+omega(2).^2.*y.^2+omega(3).^2.*z.^2);
Hr = @(k,r) 1/(2*m).*(hb.^2.*k.^2)+m/2.*(omegar.^2.*r.^2);

p = @(k,x,y,z) 1/(2*pi)^3.*1./(exp((H(k,x,y,z)-mu)./(kb.*T))+1);
pr = @(k,r) 4*pi.*1./lambda.*1/(2*pi)^3.*r.^2./(exp((Hr(k,r)-mu)./(kb.*T))+1);
p_B = @(k,x,y,z) 1/(2*pi)^3.*1./(exp((H(k,x,y,z)-mu)./(kb.*T)));
N_B = @(k) (m./(2*pi.*kb.*T)).^(3/2).*(exp(-(H(k,0,0,0))./(kb.*T)));

% set up grid
x_vec = linspace(a0(1).*-30,a0(1).*30,1.5e2);
y_vec = linspace(a0(2).*-30,a0(2).*30,1.5e2);
z_vec = linspace(a0(3).*-30,a0(3).*30,1.5e2);
r_vec = linspace(0,a0(1).*10,1.5e3).';

% Fermi at T=0
n_F_0 = @(k) N/KF^3*8/pi^2*(1-(k./KF).^2).^(3/2);

% Fermi at finite T
n_F = @(k) int_loop(p,k,x_vec,y_vec,z_vec);

% Fermi at finite T (with reduced integral
n_F_r = @(k) int_loop_r(pr,k,r_vec);

%Boltzman
n_B =  @(k) int_loop(p_B,k,x_vec,y_vec,z_vec);

k_vec = linspace(-5*KF,5*KF);

stfig('distrabution comparisons');
clf
tic
plot(k_vec,n_F(k_vec).*KF.^3./N)
toc
hold on
tic
plot(k_vec,n_F_r(k_vec).*KF.^3./N)
toc
plot(k_vec,n_B(k_vec).*KF.^3./N)
plot(k_vec,n_F_0(k_vec).*KF.^3./N)
% plot(k_vec,N_B(k_vec).*KF.^3./N)
legend('nF','nFr','nB','nF0')
ylim([0 2])

stfig('compare statistics');
clf
plot(k_vec,p(k_vec,0,0,0))
hold on
plot(k_vec,p_B(k_vec,0,0,0))

v0 = @(b,m) sqrt((2.*kb.*b(1)./m)) ; % v

  fermi_dist = @(b,x) -real(polylog(1.5,-abs(b(2)).*exp(-x.^2./(v0(b,m).^2./omegar.^2))));
% y_vec = linspace()
%%
xi_vec = linspace(0.01,50,20);
T_vec = linspace(10,1e3,20).*1e-9;
T_fix = 100e-9;
xi_fix = 1;
for ii = 1:length(xi_vec)
  y = fermi_dist([T_vec(ii),xi_fix],y_vec);
  
figure(11)
hold on
plot(y_vec,y)
width(ii) = sqrt(trapz(y_vec,y_vec.^2.*y)./trapz(y_vec,y));
end
%%
function mu = find_mu(omega_bar,T,N,EF)
global const
kb = const.kb;
hb = const.hb;

mu_temp = linspace(-4.*EF,EF,5e3).';
N_E = @(E,mu_vec) E.^2./(2*(hb.*omega_bar).^3).*1./(exp((E-mu_vec)./(kb.*T))+1);

E = linspace(0,20.*EF,5e3);

N_temp = abs(trapz(E,N_E(E,mu_temp),2)-N);
[val,indx] = min(N_temp);
mu = mu_temp(indx);

end

function n_vec = int_loop(p,k,x_vec,y_vec,z_vec)
[X,Y,Z] = meshgrid(x_vec,y_vec,z_vec);
for ii = 1:length(k)
    n_vec(ii) = trapz(x_vec,trapz(y_vec,trapz(z_vec,p(k(ii),X,Y,Z))));
end
end

function n_vec = int_loop_r(p,k,r_vec)
for ii = 1:length(k)
    n_vec(ii) = trapz(r_vec,p(k(ii),r_vec));
end
end