close all
clear all
set(groot,'DefaultTextInterpreter','latex')

%measured detected number and temperature
N = 20e3;
T = 200e-6;
qe = 0.08;

% Temperature (Kelvin)
T_vec=logspace(log10(1),log10(230)).*1e-6;
% Fall Distance
l0=0.848;
%Detector Radius
R=4/100;


%time vector
t = 100e-3:1e-6:600e-3;

%find the detected fraction for the plots
[N4_frac, N3_frac, c4, c3, c4T, c3T] = detected_fraction(t,T_vec,l0,R);

[N4_frac_measure, N3_frac_measure] = detected_fraction(t,T,l0,R);

N./[N4_frac_measure, N3_frac_measure]./qe

%Phase space density

%n0 in a cigar trap n0=N/(2 pi k T)^3 *(m wz wr^2)
%lambda = sqrt(2 pi /(kT)) *hbar
% n0 Lambda^3 = (2 pi/(k T))^1.5 * (N hbar^3 m^2.5 wz wr^2)

%%
figure(12)
clf
loglog(T_vec,N4_frac,'k-','linewidth',2)
hold on
loglog(T_vec,N3_frac,'r--','linewidth',2)
legend('He$^4$','He$^3$','interpreter','latex')
ylabel('Detected Fraction')
xlabel('Temperature (K)')
set(gcf,'color','white')
set(gca,'FontSize',19)

text(1.5e-6,0.029,['Detector Radius = ',num2str(R.*1e3),' mm'],'FontSize',19)
text(1.5e-6,0.02,['Fall Distance = ',num2str(l0.*1e3),' mm'],'FontSize',19)
grid on
%%
figure(909)
clf
T_indx = 25;
plot(t,c4(:,T_indx),'r',t,c3(:,T_indx),'b--',t,c4T(:,T_indx),'k-.',t,c3T(:,T_indx),'go','linewidth',2)
set(gcf,'color','white')
set(gca,'FontSize',19)
legend('He$^4$ detected','He$^3$ detected','He$^4$','He$^3$','interpreter','latex')
xlabel('time (s)')
ylabel('Flux (counts/s)')
text(0.2,20,['T = ',num2str(T_vec(T_indx).*1e6), ' $\mu$K'],'FontSize',19)
grid on

%% General Function to find detected fraction
function [N4_frac, N3_frac, c4, c3, c4T, c3T] = detected_fraction(t,T,l0,R)

% Constants
kb = 1.38e-23;
g = 9.8;
m4 = 4/6.022e26;
m3 = 3/6.022e26;

%%He4
A4 = (m4./(2*pi*kb.*T)).^1.5.*ones(size(t)).';
v4 = sqrt(2.*kb.*T./m4).*ones(size(t)).';

%%He3
A3 = (m3./(2*pi*kb.*T)).^1.5.*ones(size(t)).';
v3 = sqrt(2*kb.*T./m3).*ones(size(t)).';

%work out counts
f1 = ((.5*g*t.^2+l0)./t.^2).'.*ones(size(T));
f2 = ((.5*g*t.^2-l0).^2).'.*ones(size(T));
f34 = (v4.*t.').^2;
f33 = (v3.*t.').^2;

c4=A4.*pi.*v4.^2.*f1.*(1-exp(-R^2./f34)).*exp(-f2./f34);
c3=A3.*pi.*v3.^2.*f1.*(1-exp(-R^2./f33)).*exp(-f2./f33);
c4T=A4.*pi.*v4.^2.*f1.*(1-exp(-(1000*R)^2./f34)).*exp(-f2./f34);
c3T=A3.*pi.*v3.^2.*f1.*(1-exp(-(1000*R)^2./f33)).*exp(-f2./f33);
% 

N4_frac=sum(c4)./sum(c4T);
N3_frac=sum(c3)./sum(c3T);
end