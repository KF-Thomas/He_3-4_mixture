opts.data_root = 'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\';
%  opts.data_root = 'Z:\EXPERIMENT-DATA\2020_Momentum_Bells\';
% opts.data_root = 'C:\Users\BEC Machine\Documents\DATA_BACKUP\';

% data_dir = '20220406_TF_vs_evap';
data_dir = '20220729_he3_4_vs_hold_time';%'20220819_he3_he4_data_4';%
% data_dir = '20220706_large_he3_clouds_detuning';
% data_dir = '20220328_TF_dep_trial_5_fracNwizard_detuning';
% data_folder = {'20220324_he3_and_4_seperated_clouds'};
% Find all the data folders in the data directory

f_winfreak = 10767.5;%in MHz
% expected detuning 33.574 GHz
delt = 20.346-20.4;

full_data_dir = fullfile(opts.data_root, data_dir);

he4_time = [0.4392,0.46;0.4387,0.46;0.4378,0.46;0.437,0.46;0.435,0.46;0.432,0.46];%[0.435,0.46];%[1.05,1.09];
he3_time = [0.41,0.4392;0.41,0.4387;0.41,0.4378;0.41,0.437;0.41,0.435;0.41,0.432];%[0.41,0.435];%[1.01,1.05];

do_bimod = 0;

hebec_constants

bins=502;
spatial_blur=2;
tmin = 0.36;%0.9;%
tmax = 0.46;%1.0;%
ymin = -30e-3;
ymax = 30e-3;
xmin = -30e-3;
xmax =  30e-3;
time_cen = 0;
XEdges=linspace(xmin,xmax,bins);
YEdges=linspace(ymin,ymax,bins);
TEdges=linspace(tmin,tmax,bins);
bin_area=((ymax-ymin)/bins)*((xmax-xmin)/bins);

%% Plot of flux fig

for ii = [1,2]
load(fullfile(full_data_dir,['clean_data_',num2str(ii),'.mat']))
plt_indx = 2;
ax = ii;
shots_indx = 1:47;%1:23;%
initials = [1e-7 0.65 1e3 -0.012];
% variable(plt_indx)
% he4_fit_plt = he4_fits(plt_indx,:);%
% he3_fit_plt = he3_fits(plt_indx,:);%flux_he4_all{jj}
new_he3_flux = mean(cell2mat(flux_he3_all{plt_indx}(shots_indx)),2).';%new_he3_flux_all(plt_indx,:).';
new_he4_flux = mean(cell2mat(flux_he4_all{plt_indx}(shots_indx)),2).';%new_he4_flux_all(plt_indx,:).';
[val,ind] = max(new_he4_flux);
mask_he4 = ~(abs(bin_centres_he4-bin_centres_he4(ind))<0.0023);%logical(he4_masks(plt_indx,:)).';
if ax~= 1
threshold_plt = 8e-3;%-0.022;%0.0022;%0.004;
    threshold_plt_3 = -1;
else
    threshold_plt = -1;
    threshold_plt_3 = -1;
end
[he4_fit_avg,fit_flux,mask_he4,~,he4_fit_unc, he4_fit_std] = fit_model(bin_centres_he4,new_he4_flux.',initials,threshold_plt,2000,'he4','gauss',ax);
[he3_fit_avg,fit_flux_3,mask_he3,~,he3_fit_unc, he3_fit_std] = fit_model(bin_centres_he3,new_he3_flux.',initials,threshold_plt_3,2000,'he3','gauss',ax);


stfig('full flux');
if ii == 1
clf
end
hold on
box on

scale_factor = 1;
flux_scale = 1e3;
if ax ~= 1

    scale_factor = 1e3;
    subplot(3,1,3)

    % plot(bin_centres_he4(mask_he4).*scale_factor,he4_flux_avg_all(plt_indx,mask_he4),'kx',...
    %             bin_centres_he4(~mask_he4).*scale_factor,he4_flux_avg_all(plt_indx,~mask_he4),'rx',...
    %             bin_centres_he4.*scale_factor,new_he4_flux,'linewidth',2)
    grayColor = [.55 .55 .55];
    temp_flux = new_he4_flux;
    temp_flux(:,~mask_he4) = nan;
    plot(bin_centres_he4.*scale_factor,temp_flux./flux_scale,'-','linewidth',4,'Color', grayColor)
    hold on
        plot(bin_centres_he4(~mask_he4).*scale_factor,new_he4_flux(:,~mask_he4)./flux_scale,'r-','linewidth',4)
    grayColor = [.55 .55 .55];
    plot(bin_centres_he4.*scale_factor,fit_flux./flux_scale,'k--','linewidth',2.2)
text(-28,500,['t$\in$[',num2str(he4_time(1,1)+20.346-20.4),' s,',num2str(he4_time(1,2)+20.346-20.4),' s]'],'FontSize',22)
text(-28,650,'(c)','FontSize',28)
% caxis([0 1e7])
% colormap(gca,cmap)
lgd = legend('Included $^4$He$^*$','Removed $^4$He$^*$','Thermal Fit','location','southwest');
    lgd.FontSize = 25;


else
    subplot(3,1,1)
    plot((bin_centres_he4(mask_he4)+delt).*scale_factor,new_he4_flux(:,mask_he4)./flux_scale,'k-','linewidth',4);%,...
%         (bin_centres_he4(~mask_he4)+delt).*scale_factor,new_he4_flux(:,~mask_he4)./flux_scale,'rx')
    hold on
end
if ax == 2
    xlabel('$x$-axis position (mm)')
    ylabel('Flux (counts/mm)')
elseif ax == 3
    xlabel('$y$-axis position (mm)')
    ylabel('Flux (counts/mm)')
else
    xlabel('Arrival time (s)')
    ylabel('Flux (kHz)')
end


set(gca,'FontSize',19)
if ax == 1
    xlim([he3_time(1,1),he4_time(1,2)]+delt)
else
    xlim([-30 30])
end

if ax~=1
    stfig('full flux');
axes('Position',[.7 .16 .18 .14])
box on
counts_txy=cell2mat(he4_txy_all{plt_indx}.');
[counts,centers]=hist3(counts_txy(:,[3,2]),'edges',{YEdges,XEdges});
counts=counts/bin_area;

if  ~spatial_blur==0
    counts=imgaussfilt(counts,spatial_blur);
end
imagesc(centers{2}.*1e3,centers{1}.*1e3,counts.^1)
colorbar
xlim([-30,20])
ylim([-25 20])
set(gca,'Ydir','normal')
set(gcf,'Color',[1 1 1]);
% title('$x$-$y$ plane','FontSize',13)
ylabel('$y$ (mm)','FontSize',13)
xlabel('$x$ (mm)','FontSize',13)
set(gca,'FontSize',14)
end

if do_bimod
    new_he4_flux_bimodal = new_he4_flux_bimodal_all(plt_indx,:).';
    plot(bin_centres_he4.*scale_factor,new_he4_flux_bimodal,'r--','linewidth',2)
end
if ax ~= 1
    subplot(3,1,2)
    box on
    hold on
    yy2 = smooth(bin_centres_he3.*scale_factor,new_he3_flux(:,:)./flux_scale,0.05,'rloess');
    plot(bin_centres_he3.*scale_factor,yy2,'-','linewidth',4,'Color', grayColor)
%     plot(bin_centres_he3(mask_he3).*scale_factor,new_he3_flux(:,mask_he3)./flux_scale,'kx',...
%     bin_centres_he3(~mask_he3).*scale_factor,new_he3_flux(:,~mask_he3)./flux_scale,'bx')
plot(bin_centres_he3.*scale_factor,fit_flux_3./flux_scale,'k--','linewidth',2.2)
text(-28,25,['t$\in$[',num2str(he3_time(1,1)+20.346-20.4),' s,',num2str(he3_time(1,2)+20.346-20.4),' s]'],'FontSize',22)
text(-28,35,'(b)','FontSize',28)
lgd = legend('$^3$He$^*$','DFG Fit','location','southwest');
    lgd.FontSize = 25;
else
    subplot(3,1,1)
    plot((bin_centres_he3(mask_he3)+delt).*scale_factor,new_he3_flux(:,mask_he3)./flux_scale,'b-','linewidth',4);%,...
%     (bin_centres_he3(~mask_he3)+delt).*scale_factor,new_he3_flux(:,~mask_he3)./flux_scale,'bx')
    lgd = legend('$^4$He$^*$','$^3$He$^*$');
    lgd.FontSize = 25;
    text(0.358,800,'(a)','FontSize',28)
end
% plot(bin_centres_he3.*scale_factor,he3_flux_avg_all(plt_indx,:).','bx')
% plot(bin_centres_he3.*scale_factor,new_he3_flux,'bx')


if ax == 2
    xlabel('$x$-axis position (mm)')
    ylabel('Flux (counts/mm)')
elseif ax == 3
    xlabel('$y$-axis position (mm)')
    ylabel('Flux (counts/mm)')
else
    xlabel('Arrival time (s)')
    ylabel('Flux (kHz)')
end

set(gca,'FontSize',19)
if ax == 1
%     xlim([he3_time(1),he4_time(2)])
else
    xlim([-30 30])
end

if ax~=1
    stfig('full flux');
axes('Position',[.7 .46 .18 .14])
box on
counts_txy=cell2mat(he3_txy_all{plt_indx}.');
[counts,centers]=hist3(counts_txy(:,[3,2]),'edges',{YEdges,XEdges});
counts=counts/bin_area;

if  ~spatial_blur==0
    counts=imgaussfilt(counts,spatial_blur);
end
imagesc(centers{2}.*1e3,centers{1}.*1e3,counts.^1)
colorbar
xlim([-30,20])
ylim([-25 20])
set(gca,'Ydir','normal')
set(gcf,'Color',[1 1 1]);
% title('$x$-$y$ plane','FontSize',13)
ylabel('$y$ (mm)','FontSize',13)
xlabel('$x$ (mm)','FontSize',13)
set(gca,'FontSize',14)
end

end
fig = gcf;
set(fig,'Position',[1430 310 1224 1338])
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'C:\Users\BEC Machine\cloudstor\PROJECTS\He_3-4_mixture\figs\tof_example_fig_v4','-dsvg')

%% functions
% full_model_space = thermal_space(b,x,m)
% full_model_time = thermal_time(b,t,m)
function thermal = thermal_dist(b,x,m,ax)
global const
l=const.fall_distance;%0.847;
g=const.g0;
k= const.kb;
if m == 5.008237293822000e-27
    t0 = 0.37;
else
    t0 = 0.388;
end
% t0 = sqrt(2*l/g);

A = @(b,m) (m./(2*pi.*k.*b(1))).^(3/2); % A
v0 = @(b,m) sqrt((2.*k.*b(1)./m)) ; % v

part_1 = @(b,t)  (0.5*g.*(t-b(2)).^2 +l)./((t-b(2)).^2);
part_2 = @(b,t,m)  exp(-((0.5*g.*(t-b(2)).^2 -l).^2)./(v0(b,m).^2.*(t-b(2)).^2));
part_3 = @(b,x,t,m)  exp(-((x-b(4)).^2)./(v0(b,m).^2.*(t-b(2)).^2));

if ax == 1 % time axis
    thermal = b(3).*A(b,m).*v0(b,m).^2.*pi.*part_1(b,x).*part_2(b,x,m);
else % spatial axes
%     tau = linspace(0.1.*t0,2.*t0,1e3);
%     thermal = b(3)*A(b,m).*v0(b,m).*sqrt(pi).*trapz(tau,(0.5.*g.*tau.^2+l)./tau.^3.*part_2(b,tau,m).*part_3(b,x,tau,m),2);%for a spatial dimension (either x or y)
    thermal = b(3)*A(b,m).*v0(b,m).*pi.* exp(-((x-b(4)).^2)./(v0(b,m).^2.*t0.^2));
end
end
%
function fermi_dirac = fermi_dirac_dist(b,x,m,ax)
global const
% m = const.mhe*3/4;
hb = const.hb;
kb = const.kb;
l=const.fall_distance;%0.847;
g=const.g0;

omega = [60 600 600].*2*pi; %trapping frequency
omegar = 600*2*pi;
lambda = 1/10;
a0 = sqrt(hb./(m.*omega));
omega_bar = geomean(omega);
N = 1e5;%b(4);%atom number
T = b(1);%temperature

EF = hb.*omega_bar.*(6*N).^(1/3);
KF = (2*m*EF./hb^2).^0.5;%fermi wavevector
TF = EF/kb;
mu = b(4);%find_mu(omega_bar,T,N,EF);

H = @(k,x,y,z) 1/(2*m).*(hb.^2.*k.^2)+m/2.*(omega(1).^2.*x.^2+omega(2).^2.*y.^2+omega(3).^2.*z.^2);
Hr = @(k,r) 1/(2*m).*(hb.^2.*k.^2)+m/2.*(omegar.^2.*r.^2);

pr = @(k,r) 4*pi.*1./lambda.*1/(2*pi)^3.*r.^2./(exp((Hr(k,r)-mu)./(kb.*T))+1);

% set up grid
r_vec = linspace(0,a0(1).*10,1.5e3).';
r_vec_2 = linspace(0,30e-3,2.5e1);%30e-3


% Fermi at finite T (with reduced integral
vx = @(x,t) x./t;
vy = @(y,t) y./t;
vr = @(r,t) r./t;
vz = @(t) (0.5.*g.*t.^2-l)./t;
% t_vec = linspace(0.4,0.46,2e2);

n_F_r = @(k) int_loop_r(pr,k,r_vec);

n_tof = @(t,r) r.*n_F_r(m.*sqrt(vz(t).^2+vr(r,t).^2)./hb).*(0.5*g*t.^2+l)./t.^4;

fermi_dist = @(t) 2.*pi.*int_loop_r(n_tof,t,r_vec_2);

fermi_dirac = fermi_dist(x-b(2)).';
fermi_dirac = b(3).*fermi_dirac./trapz(x,fermi_dirac);
% plot(t_vec,tshermal_dist([T,0,0.001],t_vec,m,1))

end

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


function thomas_fermi = thomas_fermi_dist(b,t,ax)
global const
l=const.fall_distance;%0.847;
g=const.g0;
t0 = sqrt(2*l/g);

if ax == 1
    thomas_fermi = min(ones(size(t)).*1e10,b(5).*real(max(zeros(size(t)),(1-((t-b(2)-t0)./b(4)).^2)).^(3/2)));
else
    thomas_fermi = min(ones(size(t)).*1e10,b(6).*real(max(zeros(size(t)),(1-((t-b(4))./b(5)).^2)).^(3/2)));
end
end

%% fitting function
function [he_fit_avg,new_he_flux,mask,he_fit,he_fit_unc_avg,he_fit_std] = fit_model(bin_centre,he_flux_avg,initials,threshold,N_he,ele,model,ax)
global const
g=const.g0;
m4 = const.mhe;
m3 = 5.008237293822000e-27;
if strcmp(ele,'he4')
    t0_cen = 0.65;
    mass = m4;
else
    t0_cen = 0.627;
    mass = m3;
end

warning('off','all')

thermal = @(b,x,m) min(ones(size(x)).*1e10,thermal_dist(b,x,m,ax));
thomas_fermi = @(b,t) thomas_fermi_dist(b,t,ax);

full_model = @(b,x) min(ones(size(x)).*1e10,thermal(b,x,mass));
modelfun_he4=@(b,t) thomas_fermi(b,t) + full_model(b,t);
modelfun_he3=@(b,x) min(ones(size(x)).*1e10,fermi_dirac_dist(b,x,mass));

t_cen = trapz(bin_centre,bin_centre.*he_flux_avg)./trapz(bin_centre,he_flux_avg);
fo = statset('TolFun',10^-8,...
    'TolX',10^-10,...
    'MaxIter',10^10,...
    'UseParallel',1);

if N_he > 300%2000
    for threshold_c = threshold
        mask = abs(t_cen-bin_centre)>threshold_c;
        try
            [he_fit_c,~,J,CoV] = nlinfit(bin_centre(mask),he_flux_avg(mask),full_model,initials);
        catch
            try
                if ax == 1
                    initials = [5e-7 t0_cen max(he_flux_avg)./1e2];
                else
                    initials = [5e-7 t0_cen max(he_flux_avg)./1e2,-0.012];
                end
                [he_fit_c,~,J,CoV] = nlinfit(bin_centre,he_flux_avg,full_model,initials);
            catch
                he_fit_c = nan.*initials;
                CoV = he_fit_c;
            end
        end
        if threshold_c == threshold(1)
            he_fit = he_fit_c;
            he_fit_unc = diag(CoV.^(0.5)).';
        else
            he_fit = [he_fit;he_fit_c];
            he_fit_unc = [he_fit_unc;diag(CoV.^(0.5)).'];
        end
        he_fit_avg = nanmean(he_fit,1);
        if size(he_fit,1) == 1
            he_fit_std = he_fit.*0;
        else
            he_fit_std = nanstd(he_fit,1)./sqrt(size(he_fit,1));
        end
        he_fit_unc_avg = nanmean(he_fit_unc,1);
    end
elseif N_he > 200
    mask = ones(size(he_flux_avg));
    [he_fit,~,J,CoV] = nlinfit(bin_centre,he_flux_avg,full_model,initials);
    he_fit_avg = he_fit;
    he_fit_unc_avg = diag(CoV.^(0.5)).';
    he_fit_std = he_fit_unc_avg;
else
    mask = ones(size(he_flux_avg));
    he_fit = initials.*nan;
    he_fit_avg = he_fit;
    he_fit_unc_avg = he_fit;
    he_fit_std = he_fit_unc_avg;
end

new_he_flux = full_model(he_fit_avg,bin_centre);

if strcmp(model,'bimodal')
    if strcmp(ele,'he4')
        mask = ones(size(he_flux_avg));
        if sum(isnan(he_fit_avg))>1
            clear he_fit_avg
            T = 1e-6;
            sig_guess = 5e-3;
            mu_guess = t_cen;
            amp_guess = max(he_flux_avg)./1e2;
            if ax == 1
                initials_bimod = [T,mu_guess,amp_guess/2,sig_guess./4,max(he_flux_avg)/2];
            else
                initials_bimod = [T,0.65,amp_guess/2,mu_guess,sig_guess,max(he_flux_avg)/2];
            end
            fitobject=fitnlm(bin_centre,he_flux_avg,modelfun_he4,...
                initials_bimod,'Options',fo);
            he_fit_avg = fitobject.Coefficients.Estimate;
            new_he_flux = modelfun_he4(fitobject.Coefficients.Estimate,bin_centre);
        else
            sig_guess = sqrt(he_fit_avg(1)*c/m4)./10;
            initials_bimod = [sig_guess,max(he_flux_avg)];
            mdl_func = @(b,t) modelfun_he4( [he_fit_avg,b],t);
            fitobject=fitnlm(bin_centre,he_flux_avg,mdl_func,...
                initials_bimod,'Options',fo);
            he_fit_avg = [he_fit_avg,fitobject.Coefficients.Estimate.'];
            new_he_flux = modelfun_he4(he_fit_avg,bin_centre);
        end
    elseif strcmp(ele,'he3') % TO DO: add proper DFG function
        if sum(isnan(he_fit))>1
            T = 1e-6;
            sig_guess = 5e-3;
            mu_guess = t_cen;
            amp_guess = max(he_flux_avg)./1e2;
            initials_bimod = [T,mu_guess,amp_guess/2,max(he_flux_avg)/2];
            fitobject=fitnlm(bin_centre,he_flux_avg,modelfun_he3,...
                initials_bimod,'Options',fo);
            he_fit_avg = fitobject.Coefficients.Estimate;
            new_he_flux = modelfun_he3(fitobject.Coefficients.Estimate,bin_centre);
        else
            T = he_fit_avg(1);
            %                     sig_guess = sqrt(he_fit(1)*c/m3);
            mu_guess = he_fit_avg(2);
            amp_guess = he_fit_avg(3);
            %             initials_bimod = [T,mu_guess,amp_guess,sig_guess./4,max(he_flux_avg)/10];
            initials_bimod = [200e-9,0.0,4e2,3e-29];%[T,mu_guess,amp_guess/2];
%             fitobject=fitnlm(bin_centre.',he_flux_avg,modelfun_he3,...
%                 initials_bimod,'Options',fo);
            modelfun_he3_temp = @(b,x) modelfun_he3([200e-9,b(1),b(2),3e-29],x);
            fitobject=fitnlm(bin_centre.',he_flux_avg,modelfun_he3_temp,...
                [0.0,4e2],'Options',fo);
            he_fit_avg = fitobject.Coefficients.Estimate;
            new_he_flux = modelfun_he3(fitobject.Coefficients.Estimate,bin_centre);
        end
    end
end

warning('on','all')
end
