%% comparison of momentum density of Boltzmann and Fermi-Dirac distrabution
% in harmonic trap usng the semi classical approximation

% constats
hebec_constants
m = const.mhe*3/4;
hb = const.hb;
kb = const.kb;
g = const.g0;
l = const.fall_distance;%0.847;

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
tic
mu = find_mu(omega_bar,T,N,EF);
toc

H = @(k,x,y,z) 1/(2*m).*(hb.^2.*k.^2)+m/2.*(omega(1).^2.*x.^2+omega(2).^2.*y.^2+omega(3).^2.*z.^2);
Hr = @(k,r) 1/(2*m).*(hb.^2.*k.^2)+m/2.*(omegar.^2.*r.^2);

pr = @(k,r) 4*pi.*1./lambda.*1/(2*pi)^3.*r.^2./(exp((Hr(k,r)-mu)./(kb.*T))+1);

% set up grid
r_vec = linspace(0,a0(1).*10,1.5e3).';
r_vec_2 = linspace(0,30e-3,2.5e1);
r_vec_3 = linspace(0e-3,60e-3,4.5e1);
r_vec_4 = linspace(-60e-3,60e-3,5e1);

% Fermi at finite T (with reduced integral
vx = @(x,t) x./t;
vy = @(y,t) y./t;
vr = @(r,t) r./t;
vz = @(t) (0.5.*g.*t.^2-l)./t;
t_vec = linspace(0.4,0.46,2e2);

n_F_r = @(k) int_loop_r(pr,k,r_vec);

n_tof = @(t,r) r.*n_F_r(m.*sqrt(vz(t).^2+vr(r,t).^2)./hb).*(0.5*g*t.^2+l)./t.^4;
n_tof_2 = @(t,r2,r) r2.*pr(m.*sqrt(vz(t).^2+vr(r2,t).^2)./hb,r).*(0.5*g*t.^2+l)./t.^4;
n_tof_r = @(r,t) abs(r).*n_F_r(m.*sqrt(vz(t).^2+vr(r,t).^2)./hb).*(0.5*g*t.^2+l)./t.^4;

n_tof_xy = @(r,t) n_F_r(m.*sqrt(vz(t).^2+vr(r,t).^2)./hb).*(0.5*g*t.^2+l)./t.^4;


fermi_dist = @(t) 2.*pi.*int_loop_r(n_tof,t,r_vec_2);
fermi_dist_2 = @(t) 2.*pi.*int_loop_r_2(n_tof_2,t,r_vec_2,r_vec);
fermi_dist_r = @(r) 2.*pi.*int_loop_r(n_tof_r,r,t_vec);
fermi_dist_xy = @(r) int_loop_r(n_tof_xy,r,t_vec);

k_vec = linspace(-5*KF,5*KF);
tic
f_vec = fermi_dist(t_vec);
f_vec = 1./trapz(t_vec,f_vec).*f_vec;
toc
tic
f_vec = fermi_dist_2(t_vec);
f_vec = 1./trapz(t_vec,f_vec).*f_vec;
toc
mb_vec = 1./trapz(t_vec,thermal_dist([T,0,0.001],t_vec,m,1)).*thermal_dist([T,0,0.001],t_vec,m,1);
mb_vec_r = 1./trapz(r_vec_3,r_vec_3.*thermal_dist([T,0,0.001,0],r_vec_3,m,2)).*r_vec_3.*thermal_dist([T,0,0.001,0],r_vec_3,m,2);

fr_vec = fermi_dist_r(r_vec_3);
fr_vec = 1./trapz(r_vec_3,fr_vec).*fr_vec;

fxy_vec = fermi_dist_xy(r_vec_4);
fxy_vec = 1./trapz(r_vec_4,fxy_vec).*fxy_vec;


stfig('time of flight distribtuion');
clf
tic
plot(t_vec,f_vec)
toc
hold on
xlabel('time')
ylabel('flux')
plot(t_vec,mb_vec)

stfig('time of flight (radial) distribtuion');
clf
plot(r_vec_3,fr_vec)
hold on
plot(r_vec_3,mb_vec_r)
xlabel('radial')
ylabel('flux')

stfig('time of flight (x-y) distribtuion');
clf
plot(r_vec_4,fxy_vec)
hold on
% plot(r_vec_4,mb_vec_r)
xlabel('x-y')
ylabel('flux')

% [X,Y] = meshgrid(r_vec_4,r_vec_4);
% R = sqrt(X.^2+Y.^2);
% for ii = 1:length(r_vec_4)
% R_diff(:,:,ii) = abs(R-r_vec_4(ii));
% end
% size(min(R_diff,3))
% (r_vec_4,fxy_vec)

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


function n_vec = int_loop_r(p,k,r_vec)
for ii = 1:length(k)
    n_vec(ii) = trapz(r_vec,p(k(ii),r_vec));
end
end

function n_vec = int_loop_r_2(p,k,r_vec_2,r_vec)
for ii = 1:length(k)
    n_vec(ii) = trapz(r_vec_2,trapz(r_vec,p(k(ii),r_vec_2,r_vec)));
end
end

function thermal = thermal_dist(b,x,m,ax)
global const
l=const.fall_distance;%0.847;
g=const.g0;
k= const.kb;
t0 = sqrt(2*l/g);

A = @(b,m) (m./(2*pi.*k.*b(1))).^(3/2); % A
v0 = @(b,m) sqrt((2.*k.*b(1)./m)) ; % v

part_1 = @(b,t)  (0.5*g.*(t-b(2)).^2 +l)./((t-b(2)).^2);
part_2 = @(b,t,m)  exp(-((0.5*g.*(t-b(2)).^2 -l).^2)./(v0(b,m).^2.*(t-b(2)).^2));
part_3 = @(b,x,t,m)  exp(-((x-b(4)).^2)./(v0(b,m).^2.*(t-b(2)).^2));

if ax == 1 % time axis
    thermal = b(3).*A(b,m).*v0(b,m).^2.*pi.*part_1(b,x).*part_2(b,x,m);
else % spatial axes
    tau = linspace(0.1.*t0,2.*t0,1e3).';
    thermal = b(3)*A(b,m).*v0(b,m).*sqrt(pi).*trapz(tau,(0.5.*g.*tau.^2+l)./tau.^3.*part_2(b,tau,m).*part_3(b,x,tau,m));%for a spatial dimension (either x or y)
end
end