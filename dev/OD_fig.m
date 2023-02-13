load('OD_v_detuning.mat')
dshift=0;
stfig('Optical density vs detuning');
clf
indx = x_all>125.5&x_all<128;%1:45;%[2:9,10,12];%[2:9,10,11];%1:11;%1:18;%1:15;%

amp_func = @(b,x) b(1)./(1+4.*(2.*pi.*x(:,1)-b(2)).^2./b(3).^2);
amp_guesses = [2,126.5*2*pi,5.2];
fitobjecta=fitnlm(x_all(indx).',y_all(indx).',amp_func,amp_guesses)
xt=linspace(min(x_all(indx))*0.95,max(x_all(indx))*1.05);
amp_fit=predict(fitobjecta,xt.');
plot(xt-dshift,amp_fit,'LineWidth',2)
hold on
scatter(x_all(indx)-dshift,y_all(indx),'ko')
box on
xlabel('$\Delta$ ($2\pi$MHz)')
ylabel('$D_0$')
set(gca,'FontSize',17)
xlim([min(x_all(indx))*0.95,max(x_all(indx))*1.05]-dshift)
%%
load('OD_v_inten.mat')
stfig('amp vs power');
clf
indx2 = [1,3:11];
power = (ref_max);%ref_mean;%
% int_y = (amp).*(2*pi.*cell2mat(x).*cell2mat(y));
amp_y = amp;
damp_y = sqrt(damp.^2+offset.^2+doffset.^2);
amp_y(1) = 2.6;
damp_y(1) = 0.25;
weights = 1./damp_y.^2;
% errorbar(power,int_y,int_y.*sqrt((damp./amp).^2+(cell2mat(dx)./cell2mat(x)).^2+(cell2mat(dy)./cell2mat(y)).^2),'kx')

amp_func = @(b,x) b(1)./(1+(x).^b(3).*b(2));
amp_guesses = [0.85,11,2];
fitobjecta=fitnlm(power(indx2),amp_y(indx2).',amp_func,amp_guesses,'Weights',weights(indx2))
xt=linspace(0,max(power).*1.1);

amp_fit=predict(fitobjecta,xt.');

plot(xt./max(power),amp_fit,'linewidth',2)
hold on
errorbar(power(indx2)./max(power),amp_y(indx2),damp_y(indx2),'kx')
box on
xlim([0,1.1])
xlabel('$\sqrt{I_0}$ (arb. unit)')
ylabel('$D_0$')
set(gca,'FontSize',17)
%%
stfig('OD example');
load('OD_exampe.mat')
xy_factor=1e6;

        font_name='cmr10';
        linewidth=1.5;
        font_size=12;
pcolor((xvals_sub+1.47904828029349e-03)*xy_factor,(yvals_sub-7.67236372643944e-04)*xy_factor,OD_meas)
caxis([0 1.5])
hcb=colorbar;
colorTitleHandle = get(hcb,'Title');
titleString = '$D(x,y)$';
set(colorTitleHandle ,'String',titleString,'Interpreter','latex');
ifh=gca;
    shading interp
    sh.EdgeColor='k';
   box on
    xlabel('x ($\mathrm{\mu}$m)')
    ylabel('y ($\mathrm{\mu}$m)')
    zlabel('Optical density')
    shading interp
    set(gca,'FontSize',17)
    xlim([-400 400])
    ylim([-400 400])