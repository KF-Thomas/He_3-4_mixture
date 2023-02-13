clear all
% constats
const.c = 299792458; %speed of light (m/s)
const.h = 6.626070040*10^-34; %Planck constant (J s)
const.hb = 1.054*10^-34; %reduced Planck constant (J S)
const.kb = 1.3806488*10^-23; %Boltzmann constant (*m2 kg s-2 K-1*)
const.mu0 = 1.2566370614*10^-6;  % vacuum permeability [Tm/A]
const.epsilon0 = 8.858*10^-12;%(*electric permittivity of free space*)
%elemental
const.mub =9.274009994*10^-24; %Bohr magneton*(J/T)
const.electron = 1.60217657*10^-19;%(*charge of electron*)
const.me = 9.10938291*10^-31;%(*mass of electron*)
const.grav=6.67430*10^-11;  %Newtonian constant of gravitation %https://physics.nist.gov/cgi-bin/cuu/Value?bg

%Helium
% const.ahe_scat=15*10^-9;
const.ahe_scat=7.512000000000000e-09; %m^2 Moal et al https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.96.023203
const.b_freq=2.802*1e6*1e4; %hz/T
const.mhe = 1.66*10^-27*4.002;%(*helium mass*)
const.interaction_energy = 4*pi*const.hb^2*const.ahe_scat/const.mhe; % interaction strength []

%environmental
const.g0 =9.7960207-78e-8;
% (-35 17.0', 149 6.8) alt=560.78m
% https://d28rz98at9flks.cloudfront.net/15167/Rep_261.pdf 
%with correction of 78e-8 from http://www.publish.csiro.au/ex/ASEG2007ab149 applied
%9.796 cbr approx grav at (-35.283103 , 149.112634) from http://ddfe.curtin.edu.au/gravitymodels/GGMplus/data/ga/S40E145.ga.png

%customary
const.a0 = 5.29177210903e-11;%bohr radius (m) (0.00000000080e-11) %https://physics.nist.gov/cgi-bin/cuu/Value?Abohrrada0

%experiment
const.fall_distance = 0.8587; 
%%
mb = const.mhe;
mf = const.mhe*3/4;
mr = mb*mf/(mb+mf);
hb = const.hb;
h = const.h;
kb = const.kb;
omegab = 600*2*pi;
omegaf = 600*sqrt(4/3)*2*pi;
lambdab = 1/10;
lambdaf = 1/10;

abf = 28.8e-9;%0e-9;%

g = 4*pi*hb^2*const.ahe_scat/mb;
f = 2*pi*hb^2*abf/mr;

omega = [60 600 600].*2*pi; %trapping frequency
a0 = sqrt(hb./(mb.*omega));
omega_bar = geomean(omega);
a0_bar = sqrt(hb./(mb.*omega_bar));

Nb = 550e3;%atom number
Nf = 20e3;%9e3atom number
T = 180e-9;%temperature

%%
EF = hb.*omega_bar.*(6*Nf).^(1/3)*sqrt(4/3);
KF = (2*mf*EF./hb^2).^0.5;%fermi wavevector
TF = EF/kb;
Rf = (2*EF/(mf*omegaf^2))^0.5;


mub = hb.*omega_bar./2.*(15*const.a0*Nb/a0_bar).^(2/5);
muf = EF;%mub;%find_mu(omega_bar,T,Nf,EF)/2;

%% functions
Vb_ext = @(r,z) mb.*omegab.^2.*(r.^2+lambdab.^2.*z.^2)./2;
Vf_ext = @(r,z) mf.*omegaf.^2.*(r.^2+lambdaf.^2.*z.^2)./2;
H = @(p,m,V,mu) p.^2./(2.*m)+V-mu;
fn = @(p,m,V,mu,d) max(1/(2*pi*hb)^3.*1./(exp(H(p,m,V,mu)./(kb.*T))+d),0);

%% set up grid


%% initial guesses
options = optimset('PlotFcns',@optimplotfval);
x0 = [mub,muf];%[8.936824706016943e-31     5.489373962454873e-30]
x = fminsearch(@fun,x0,options);
%%
mub = x(1);
muf = x(2);



p_vec_1 = linspace(0,5*KF);
p_vec_2 = linspace(sqrt(mub*2*mb).*0.05,sqrt(mub*2*mb)*5);
p_vec = reshape(p_vec_1,[1,1,length(p_vec_1)]);
p_vec_b = reshape(p_vec_2,[1,1,length(p_vec_2)]);

r_vec = linspace(0,Rf.*1.5,1.0e2).';
z_vec = linspace(-Rf.*15,Rf.*15,1.5e2).';

[R,Z] = meshgrid(r_vec,z_vec);

nc = (mub-Vb_ext(R,Z))./g;
nc(nc<0) = 0;
%trapz(r_vec,trapz(z_vec,R.*nc)).*2*pi
nf = Nf.*lambdaf./(Rf.^3).*8./pi^2.*(1-(R.^2+lambdaf.^2.*Z.^2)./Rf.^2).^(3/2);
nf(real(nf)<0|abs(imag(nf))>0)=0;
nf_int = nf;
nnc = nf.*0;

%% evaluate
for ii = 1:20
    Vb_eff = Vb_ext(R,Z) + 2.*g.*nc + 2.*g.*nnc + f.*nf;
    Vf_eff = Vf_ext(R,Z) + f.*nc + f.*nnc;
nc = (mub -Vb_ext(R,Z) -2*g*nnc-f*nf)./g;
nc(nc<0) = 0;
nnc = trapz(p_vec_2,4.*pi.*p_vec_b.^2.*fn(p_vec_b,mb,Vb_eff,mub,-1),3);
nf = trapz(p_vec_2,4.*pi.*p_vec_b.^2.*fn(p_vec_b,mf,Vf_eff,muf,1),3);
nf(real(nf)<0|abs(imag(nf))>0)=0;

end

Vb_eff = Vb_ext(R,Z) + 2.*g.*nc + 2.*g.*nnc + f.*nf;
    Vf_eff = Vf_ext(R,Z) + f.*nc + f.*nnc;

%%
nf_k = squeeze(trapz(r_vec,r_vec.'.*trapz(z_vec,2.*pi.*fn(p_vec_b,mf,Vf_eff,muf,1))));
nb_k = squeeze(trapz(r_vec,r_vec.'.*trapz(z_vec,2.*pi.*fn(p_vec_b,mb,Vb_eff,mub,-1))));
nf_k_norm = squeeze(trapz(r_vec,r_vec.'.*trapz(z_vec,2.*pi.*fn(p_vec_b,mf,Vf_ext(R,Z),muf,1))));
nf_k_therm = squeeze(trapz(r_vec,r_vec.'.*trapz(z_vec,2.*pi.*fn(p_vec_b,mf,Vf_ext(R,Z),0,0))));
fn_T = @(p,m,V,mu,d,T) max(1/(2*pi*hb)^3.*1./(exp(H(p,m,V,mu)./(kb.*T))+d),0);
nf_k_T = squeeze(trapz(r_vec,r_vec.'.*trapz(z_vec,2.*pi.*fn_T(p_vec_b,mf,Vf_ext(R,Z),muf,1,195e-9))));

%% plot
figure(2)
clf
surf(R,Z,nf)
xlabel('radial')
ylabel('z')

figure(11)
clf
plot(z_vec.*1e6,nc(:,1))
hold on
plot(z_vec.*1e6,nf(:,1))
plot(z_vec.*1e6,nnc(:,1))
title('r=0')
xlabel('z (mum)')
set(gca,'FontSize',18)

figure(4)
clf
plot(z_vec,nnc(:,1))
xlabel('z')
set(gca,'FontSize',18)

%%
figure(3)
clf
plot(squeeze(p_vec_b),nf_k./trapz(squeeze(p_vec_b),nf_k))
hold on
plot(squeeze(p_vec_b),nb_k./trapz(squeeze(p_vec_b),nb_k))
plot(squeeze(p_vec_b),nf_k_norm./trapz(squeeze(p_vec_b),nf_k_norm))
nf_k_zero = Nf./KF.^3.*8/pi^2*(1-(squeeze(p_vec_b)./hb).^2./KF.^2).^(1.5);
nf_k_zero(real(nf_k_zero)<0|abs(imag(nf_k_zero))>0)=0;
plot(squeeze(p_vec_b),nf_k_zero./trapz(squeeze(p_vec_b),nf_k_zero))
plot(squeeze(p_vec_b),nf_k_T./trapz(squeeze(p_vec_b),nf_k_T))

n_k_theraml = exp(-squeeze(p_vec_b).^2./(2.*mf.*kb.*141.6e-9));
plot(squeeze(p_vec_b),n_k_theraml./trapz(squeeze(p_vec_b),n_k_theraml))

xlabel('momentum')
legend('fermi momentum','no boson','zero temp no boson','lower temp no bosons')

%%
figure
nb_therm_pol = polylog(3/2,exp((-mub-p_vec_2.^2/(2*mb))./(kb*T)));
plot(p_vec_2,nb_therm_pol./trapz(p_vec_2,nb_therm_pol))
hold on
nT = @(b) squeeze(trapz(r_vec,r_vec.'.*trapz(z_vec,2.*pi.*fn_T(p_vec_b,mb,Vb_ext(R,Z),b(2),-1,b(1)))));
% nT = @(b) squeeze(trapz(r_vec,r_vec.'.*trapz(z_vec,2.*pi.*fn_T(p_vec_b,mb,Vb_ext(R,Z),0,0,b))));
modelfun = @(b,x) nT(b)./trapz(px,nT(b));%exp(-x.^2./(2.*b.*mf.*kb)).*1./(sqrt(2*pi*b*mf*kb));%
% plot(px,modelfun([T,-mub],px))
%%
px = squeeze(p_vec_b);%p_vec_b;%
num_pts = 30;
thr = linspace(0.0,5,num_pts);
for jj = 1:num_pts
mask = px>=thr(jj)*1e-28;
ny = nf_k./trapz(p_vec_2(mask),nf_k(mask));
nT = @(b,x) exp(-x.^2./(2.*b.*mf.*kb));%squeeze(trapz(r_vec,r_vec.'.*trapz(z_vec,2.*pi.*fn_T(p_vec_b(mask),mf,Vf_ext(R,Z),muf,1,b))));
modelfun = @(b,x) nT(b,x)./trapz(x,nT(b,x));%.*1./(sqrt(2*pi*b*mf*kb));%

nTf = @(b,x) squeeze(trapz(r_vec,r_vec.'.*trapz(z_vec,2.*pi.*fn_T(p_vec_b(mask),mf,Vf_ext(R,Z),muf*0.5,1,b))));
% nTf = @(b,x) squeeze(trapz(r_vec,r_vec.'.*trapz(z_vec,2.*pi.*fn_T(p_vec_b(mask),mf,Vf_ext(R,Z),b,1,T))));
modelfun_f = @(b,x) nTf(b,x)./trapz(x,nTf(b,x));

beta0 = 200e-9;
% beta0mu = muf;

mdl = fitnlm(px(mask),ny(mask),modelfun,beta0)
mdl_f = fitnlm(px(mask),ny(mask),modelfun_f,beta0)
vec_he3_T(jj) = mdl.Coefficients.Estimate;
vec_he3(jj) = mdl_f.Coefficients.Estimate;
end
%%
px = squeeze(p_vec_b);%p_vec_b;%
thr = linspace(0.0,5,30);
for jj = 1:30 
mask = px>thr(jj)*1e-28;
ny = nb_k./trapz(p_vec_2(mask),nb_k(mask));
nT = @(b,x) exp(-x.^2./(2.*b.*mb.*kb));%real(polylog(3/2,exp((-p_vec_2(mask).^2/(2*mb))./(kb*max(1e-9,b(1)))))).';%squeeze(trapz(r_vec,r_vec.'.*trapz(z_vec,2.*pi.*fn_T(p_vec_b,mb,Vb_ext(R,Z),-abs(mub/10000),-1,b(1)))));
% nT = @(b) squeeze(trapz(r_vec,r_vec.'.*trapz(z_vec,2.*pi.*fn_T(p_vec_b,mb,Vb_ext(R,Z),0,0,b))));
modelfun = @(b,x) max(0,nT(b,x)./trapz(x,nT(b,x)));%exp(-x.^2./(2.*b.*mf.*kb)).*1./(sqrt(2*pi*b*mf*kb));%

nTb = @(b,x) real(polylog(3/2,exp((-p_vec_2(mask).^2/(2*mb))./(kb*max(1e-9,b(1)))))).';%squeeze(trapz(r_vec,r_vec.'.*trapz(z_vec,2.*pi.*fn_T(p_vec_b,mb,Vb_ext(R,Z),-abs(mub/10000),-1,b(1)))));
% nT = @(b) squeeze(trapz(r_vec,r_vec.'.*trapz(z_vec,2.*pi.*fn_T(p_vec_b,mb,Vb_ext(R,Z),0,0,b))));
modelfun_b = @(b,x) max(0,nTb(b,x)./trapz(x,nTb(b,x)));%exp(-x.^2./(2.*b.*mf.*kb)).*1./(sqrt(2*pi*b*mf*kb));%


beta0 = [55.5e-9];

mdl = fitnlm(px(mask),ny(mask),modelfun,beta0)
mdl_b = fitnlm(px(mask),ny(mask),modelfun_b,beta0)
vec_T(jj) = mdl.Coefficients.Estimate;
vec(jj) = mdl_b.Coefficients.Estimate;
end
%%
stfig('threshold vs fit temperature');
clf
plot(thr./sqrt(mb*kb*T).*1e-28,vec)
hold on
plot(thr./sqrt(mb*kb*T).*1e-28,vec_T)
plot(thr./sqrt(mf*kb*T).*1e-28,vec_he3_T)
plot(thr./sqrt(mf*kb*T).*1e-28,vec_he3)
legend('full BE','thermal bose','thermal fermi','full FD')
%%
function N_err = fun(x)
mub = x(1);
muf = x(2);
hebec_constants

mb = const.mhe;
mf = const.mhe*3/4;
mr = mb*mf/(mb+mf);
hb = const.hb;
kb = const.kb;
omegab = 600*2*pi;
omegaf = 600*sqrt(4/3)*2*pi;
lambdab = 1/10;
lambdaf = 1/10;

abf = 28.8e-9;%0e-9;%

g = 4*pi*hb^2*const.ahe_scat/mb;
f = 2*pi*hb^2*abf/mr;

omega = [60 600 600].*2*pi; %trapping frequency
a0 = sqrt(hb./(mb.*omega));
omega_bar = geomean(omega);
a0_bar = sqrt(hb./(mb.*omega_bar));

Nb = 1e4;%atom number
Nf = 1e4;%atom number
T = 200e-9;%temperature

%%
EF = hb.*omega_bar.*(6*Nf).^(1/3)*sqrt(3/4);
KF = (2*mf*EF./hb^2).^0.5;%fermi wavevector
TF = EF/kb;
Rf = (2*EF/(mf*omegaf^2))^0.5;


%% functions
Vb_ext = @(r,z) mb.*omegab.^2.*(r.^2+lambdab.^2.*z.^2)./2;
Vf_ext = @(r,z) mf.*omegaf.^2.*(r.^2+lambdaf.^2.*z.^2)./2;
H = @(p,m,V,mu) p.^2./(2.*m)+V-mu;
fn = @(p,m,V,mu,d) max(1/(2*pi*hb)^3.*1./(exp(H(p,m,V,mu)./(kb.*T))+d),0);

%% set up grid
p_vec_1 = linspace(0,5*KF);
p_vec_2 = linspace(sqrt(mub*2*mb).*0.05,sqrt(mub*2*mb)*5);
p_vec = reshape(p_vec_1,[1,1,length(p_vec_1)]);
p_vec_b = reshape(p_vec_2,[1,1,length(p_vec_2)]);

r_vec = linspace(0,Rf.*3,1.0e2).';
z_vec = linspace(-Rf.*15,Rf.*15,1.5e2).';

[R,Z] = meshgrid(r_vec,z_vec);

nc = (mub-Vb_ext(R,Z))./g;
nc(nc<0) = 0;
%trapz(r_vec,trapz(z_vec,R.*nc)).*2*pi
nf = Nf.*lambdaf./(Rf.^3).*8./pi^2.*(1-(R.^2+lambdaf.^2.*Z.^2)./Rf.^2).^(3/2);
nf(real(nf)<0|abs(imag(nf))>0)=0;
nf_int = nf;
nnc = nf.*0;



%% evaluate
for ii = 1:20
    Vb_eff = Vb_ext(R,Z) + 2.*g.*nc + 2.*g.*nnc + f.*nf;
Vf_eff = Vf_ext(R,Z) + f.*nc + f.*nnc;
nc = (mub -Vb_ext(R,Z) -2*g*nnc-f*nf)./g;
nc(nc<0) = 0;
nnc = trapz(p_vec_2,4.*pi.*p_vec_b.^2.*fn(p_vec_b,mb,Vb_eff,mub,-1),3);
nf = trapz(p_vec_2,4.*pi.*p_vec_b.^2.*fn(p_vec_b,mf,Vf_eff,muf,1),3);
nf(real(nf)<0|abs(imag(nf))>0)=0;

end

Nb_c = trapz(r_vec,trapz(z_vec,R.*nc)).*2*pi+trapz(r_vec,trapz(z_vec,R.*nnc)).*2*pi;
Nf_c = trapz(r_vec,trapz(z_vec,R.*nf)).*2*pi;

N_err = abs(Nb_c-Nb)./Nb+abs(Nf_c-Nf)./Nf;
end

function mu = find_mu(omega_bar,T,N,EF)

mu_temp = linspace(-4.*EF,EF,5e3).';
hb = 6.626070040*10^-34;
kb = 1.3806488*10^-23;
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

function [y errors] = polylog(n,z) 
%% polylog - Computes the n-based polylogarithm of z: Li_n(z)
% Approximate closed form expressions for the Polylogarithm aka de 
% Jonquiere's function are used. Computes reasonably faster than direct
% calculation given by SUM_{k=1 to Inf}[z^k / k^n] = z + z^2/2^n + ...
%
% Usage:   [y errors] = PolyLog(n,z)
%
% Input:   z < 1   : real/complex number or array
%          n > -4  : base of polylogarithm 
%
% Output: y       ... value of polylogarithm
%         errors  ... number of errors 
%
% Approximation should be correct up to at least 5 digits for |z| > 0.55
% and on the order of 10 digits for |z| <= 0.55!
%
% Please Note: z vector input is possible but not recommended as precision
% might drop for big ranged z inputs (unresolved Matlab issue unknown to 
% the author). 
%
% following V. Bhagat, et al., On the evaluation of generalized
% Bose胞instein and Fermi縫irac integrals, Computer Physics Communications,
% Vol. 155, p.7, 2003
%
% v3 20120616
% -------------------------------------------------------------------------
% Copyright (c) 2012, Maximilian Kuhnert
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:  
% 
%     Redistributions of source code must retain the above copyright
%     notice, this list of conditions and the following disclaimer. 
%     Redistributions in binary form must reproduce the above copyright
%     notice, this list of conditions and the following disclaimer in the
%     documentation and/or other materials provided with the distribution.  
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.          
% -------------------------------------------------------------------------
if nargin~=2
    errors=1;
    error('[Error in: polylog function] Inappropriate number of input arguments!')
end
if (isreal(z) && sum(z(:)>=1)>0) % check that real z is not bigger than 1 
    errors=1;
    error('[Error in: polylog function] |z| > 1 is not allowed')
elseif isreal(z)~=1 && sum(abs(z(:))>1)>0 % check that imaginary z is defined on unit circle
    errors=1;
    error('[Error in: polylog function] |z| > 1 is not allowed')
elseif n<=-4 % check that n is not too largly negative (see paper)
    errors=1;
    error('[Error in: polylog function] n < -4 might be inaccurate')
end
% display more digits in Matlab terminal:
%format long
alpha = -log(z); % see page 12
% if |z| > 0.55 use Eq. (27) else use Eq. (21):
if abs(z) > 0.55
    preterm = gamma(1-n)./alpha.^(1-n);
    nominator = b(0) + ...
        - alpha.*( b(1) - 4*b(0)*b(4)/7/b(3) ) + ...
        + alpha.^2.*( b(2)/2 + b(0)*b(4)/7/b(2) - 4*b(1)*b(4)/7/b(3) ) + ...
        - alpha.^3.*( b(3)/6 - 2*b(0)*b(4)/105/b(1) + b(1)*b(4)/7/b(2) - 2*b(2)*b(4)/7/b(3) );
    denominator = 1 + alpha.*4*b(4)/7/b(3) +...
        + alpha.^2.*b(4)/7/b(2) +...
        + alpha.^3.*2*b(4)/105/b(1) +...
        + alpha.^4.*b(4)/840/b(0);
    y = preterm + nominator ./ denominator;
else
    nominator = 6435*9^n.*S(n,z,8) - 27456*8^n*z.*S(n,z,7) + ...
        + 48048*7^n*z.^2.*S(n,z,6) - 44352*6^n*z.^3.*S(n,z,5) + ...
        + 23100*5^n*z.^4.*S(n,z,4) - 6720*4^n.*z.^5.*S(n,z,3) + ...
        + 1008*3^n*z.^6.*S(n,z,2) - 64*2^n*z.^7.*S(n,z,1);
    denominator = 6435*9^n - 27456*8^n*z + ...
        + 48048*7^n*z.^2 - 44352*6^n*z.^3 + ...
        + 23100*5^n*z.^4 - 6720*4^n*z.^5 + ...
        + 1008*3^n*z.^6 - 64*2^n*z.^7 + ...
        + z.^8;
    y = nominator ./ denominator;
end
% define b:
    function out = b(i)
        out = zeta(n-i);
    end
% define S as partial sums of Eq. 12:
    function out = S(n,z,j)
        out =0;
        for i=1:j
            out = out + z.^i./i^n;
        end
    end
    function [out] = zeta(x)
        %% Zeta Function  
        % Eq. 18
        % following V. Bhagat, et al., On the evaluation of generalized
        % Bose胞instein and Fermi縫irac integrals, Computer Physics Communications,
        % Vol. 155, p.7, 2003
        %
        % Usage: [out] = zeta(x)
        % with argument x and summation from 1 to j
        %
        % MK 20120615
        prefactor = 2^(x-1) / ( 2^(x-1)-1 );
        numerator = 1 + 36*2^x*eta(x,2) + 315*3^x*eta(x,3) + 1120*4^x*eta(x,4) +...
            + 1890*5^x*eta(x,5) + 1512*6^x*eta(x,6) + 462*7^x*eta(x,7);
        denominator = 1 + 36*2^x + 315*3^x + 1120*4^x + 1890*5^x + 1512*6^x +...
            + 462*7^x;
        out = prefactor * numerator / denominator;
        function [out] = eta(x,j)
            %% Eta Function  
            % Eq. 17 (partial sums)
            % following V. Bhagat, et al., On the evaluation of generalized
            % Bose胞instein and Fermi縫irac integrals, Computer Physics Communications,
            % Vol. 155, p.7, 2003
            %
            % Usage: [out] = eta(x,j)
            % with argument x and summation from 1 to j
            %
            % MK 20120615
            
            out=0;
            for k=1:j
                out = out + (-1)^(k+1) ./ k.^x;
            end
        end
    end
end
