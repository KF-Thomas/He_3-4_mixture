% beam profile difference


%fit probe beam
% user param

%im_raw=imread('./data/beam_profile/20180724T211918.841.png');
%im_raw=imread('C:\Users\BEC Machine\cloudstor\PROJECTS\Tune_out_v2_git\data\beam_profile\output_of_fiber_6.png');

% testfiledir='C:\Users\BEC Machine\cloudstor\PROJECTS\Tune_out_v2_git\data\beam_profile\input_slide_rail_mon_20220228\';
% filedir_1='Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\beam_images\on';
% filedir_2='Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\beam_images\off';
fildedir_dark = '';
filedir_comp = 'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\beam_images\no_atoms_ref';
filedir = 'Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\beam_images\evap_refs';

files = dir(filedir);
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags); % A structure with extra info.
% Get only the folder names into a cell array.
subFolderNames = {subFolders(3:end).name}; % Start at 3 to skip . and ..
subFolderNames = {filedir_comp,subFolderNames{:}};

% dist=[100,105,110,55,60,65,70,75,80,85,90,95];
% dist=[100,120,160,170,180,185,190,192,193,198,50];

% dist= [6.9,6.92,6.93,6.95,6.98,6.99,7.0,7.01,7.018,7.02,7.03,7.05,7.08,7.1];

dist = [6.1,6.25,6.3,6.4,6.5,6.6,6.7];
I_lim = [1,4,6,9,11,9,17];

%need Ilight, I atoms, and Idark
hebec_constants
%physical parameters
%we use the 23S1 to 23P2 transition
Gamma=1.0216e+07;%Hz transition linewidth
det=0;%detuning
dipole_moment = 3;
omega0=2*pi*2.767321865671393e+14;
%general case for Isat and sigma0
Isat=const.c*const.epsilon0*Gamma^2*const.hb^2/(4*dipole_moment);
sigma0=const.hb*omega0*Gamma/(2*Isat);
%sigma+ light and two level system
Isat=omega0^3*Gamma*const.hb/(12*pi*const.c^2);
lambda0=2*pi*const.c/omega0;
sigma0=3*lambda0^2/(2*pi);

I0=0.*Isat;

sigma=sigma0/(1+4.*(det/Gamma).^2+I0/Isat);

OD_sat = 1.58;%10;%0.6;%


x=[];
y=[];
% dist = [310,320,330,350,400,450,500,550];
% dist = [310,320,330,380,400,420,450,500,580];
% dist = 1:10;
% dist = [330, 340, 350, 400, 450, 500, 550, 560, 570];
% dist = [0,100,150,200,250,300,350,400,450,500,550,580,50];
wavelen = 1083.3e-9;
pix_size=20e-6;%3e-6;
% close('all')
% files = {'C:\Users\BEC Machine\cloudstor\PROJECTS\Tune_out_v2_git\data\beam_profile\input_slide_rail_0.png'};

for ii = 1:length(subFolderNames)
    if ii == 1
filedir_c = subFolderNames{ii};
    else
filedir_c = fullfile(filedir,subFolderNames{ii});
    end

files_total = dir(fullfile(filedir_c, '*.png'));
files_total = {files_total.name};
num_shots = size(files_total,2);

for j=1:2
    files = files_total((0+j):2:(num_shots-2+j));
    
    for i=1:length(files)
        im_raw=imread(fullfile(filedir_c,files{i}));
        %     im_raw=imread(fullfile(files{i}));


        %%
        beam_image=double(im_raw)/255;
        % beam_image=sum(beam_image,3);
        if i == 1
            total_image{ii,j} = beam_image;
        else
            total_image{ii,j} = total_image{ii,j}+beam_image;
        end



        %     % find the mean of the image excluding this subframe
        %     im_background=sum(beam_image(:))-sum(subframe(:));
        %     im_background=im_background/(numel(beam_image)-numel(subframe));
        %     subframe=subframe-im_background;
        %     subframe=subframe./max(subframe(:));
    end
    total_image{ii,j} = total_image{ii,j}./length(files);
    %     total_image{j}=total_image{j}./max(max(total_image{j}));
    sub_size=size(total_image{ii,j});
    sub_size_m=pix_size*sub_size;
    xvals_sub=linspace(-sub_size(2)/2,sub_size(2)/2,sub_size(2))*pix_size;
    yvals_sub=linspace(-sub_size(1)/2,sub_size(1)/2,sub_size(1))*pix_size;


    %%

    x_lim=[-1,1]*30e-6;
    y_lim=x_lim;
    xy_factor=1e6;

    font_name='cmr10';
    linewidth=1.5;
    font_size=12;

    stfig(['camera dir',num2str(j)]);
    clf
    set(gca, 'FontName', font_name)
    set(gca, 'FontSize', font_size)



    sh=surf(xvals_sub*xy_factor,yvals_sub*xy_factor,total_image{ii,j},...
        'FaceAlpha',0.4);
    ifh=gca;
    shading interp
    sh.EdgeColor='k';
    if j ==1
        colormap(gca,plasma())
    end
    %     caxis([0 0.8])
    %pbaspect([1,1,1])
    % xlim(x_lim*xy_factor)
    % ylim(y_lim*xy_factor)
    ifh.View= [7.8   14.04];
    box on
    xlabel('y ($\mathrm{\mu}$m)')
    ylabel('z ($\mathrm{\mu}$m)')
    zlabel('Intensity (arb. u.)')
    shading interp

    stfig(['integrated profile dir:',num2str(j)]);
    if ii==1
    clf
    end
    subplot(2,1,1)
    hold on
    plot(xvals_sub*xy_factor,trapz(yvals_sub,total_image{ii,j}));
    xlabel('y ($\mathrm{\mu}$m)')
    ylabel('Integrated Intensity (arb. u.)')
    subplot(2,1,2)
    hold on
    plot(yvals_sub*xy_factor,trapz(xvals_sub,total_image{ii,j},2));
    xlabel('z ($\mathrm{\mu}$m)')
    ylabel('Integrated Intensity (arb. u.)')
    if j == 1
        ref_max(ii) = max(trapz(xvals_sub,total_image{ii,j},2));
    end

    %%
% 
%     fig_name=['beam_profile_img_',files{i}];
%     %fig_name='pal_mean_v_acc_dyn_static';
%     fig_dir='C:\Users\BEC Machine\cloudstor\PROJECTS\Tune_out_v2_git\data\beam_profile\figs';
    %     export_fig(fullfile(fig_dir,strcat(fig_name,'.svg')))
    %     export_fig(fullfile(fig_dir,strcat(fig_name,'.png')))



    %%


    %colors_main=[[98,136,63];[95,109,187];[180,72,117]]./255;
    colors_main=[[40,136,40];[95,109,187];[218, 129, 57]]./255;

    hsv=colorspace('RGB->HSV',colors_main(:,:));
    hsv(:,2)=hsv(:,2);
    colors_shaded=colorspace('HSV->RGB',hsv);

end
%%
if ii>1
stfig('optical density 1');
clf
set(gca, 'FontName', font_name)
set(gca, 'FontSize', font_size)

I_light = total_image{ii,2};%total_image{ii,2};
I_atoms = total_image{ii,1};

OD_meas = log(I_light./I_atoms);%log((I_light-I_dark)./(I_atoms-I_dark));measured optical density
%basic thresholding
% A(abs(A)>3)=nan;
% power htesholding
I_min = 10;%I_lim(ii); %[1,4,4,6,9,11,9,17]
OD_meas(((I_atoms<I_min)|(I_light<I_min))) = nan;
% spike and averaging and interpolation

OD_mod = log((1-exp(-OD_sat))./(exp(-OD_meas)-exp(-OD_sat)));
OD_actual = OD_mod+(1-exp(-OD_mod)).*I0./Isat;
OD_actual(isnan(OD_actual)) = 0;


maximum = max(max(OD_meas));
[x_indx,y_indx]=find(OD_meas==maximum);

[xmesh,ymesh]=meshgrid(xvals_sub,yvals_sub);
    pos_val = zeros(size(xmesh,1),size(ymesh,2),3);
    pos_val(:,:,1) = xmesh;
    pos_val(:,:,2) = ymesh;
    pos_val(:,:,3) = OD_meas; 
    pos_val_rows=reshape(pos_val,[],3,1);
    % stfig('conversion check')
 opt = statset('TolFun',1e-10,'TolX',1e-10,...
        'MaxIter',1e4,... %1e4
        'UseParallel',1);
 r_lim = 750e-6;
  r_mask = (sqrt((pos_val(:,:,1)-xvals_sub(y_indx)).^2+(pos_val(:,:,2)-yvals_sub(x_indx)).^2))<r_lim;
    r_mask = logical(r_mask);
OD_meas(~r_mask) = nan;

sh=surf(xvals_sub*xy_factor,yvals_sub*xy_factor,OD_meas,...
    'FaceAlpha',0.75);
ifh=gca;
shading interp
sh.EdgeColor='k';
%     caxis([0 0.8])
%pbaspect([1,1,1])
% xlim(x_lim*xy_factor)
% ylim(y_lim*xy_factor)
ifh.View= [7.8   14.04];
box on
xlabel('y ($\mathrm{\mu}$m)')
ylabel('z ($\mathrm{\mu}$m)')
zlabel('Optical density')
shading interp
    
    cof_names={'sigma x','sigma y', 'cen x','cen y', 'theta', 'offset', 'amp'};
    
    modelfun = @(b,x) gauss_2d_rot_function(x,[b(1),b(2)],[b(3),b(4)],b(5),b(6),b(7),'sigma');
    
    fit_params_guess = [std(pos_val_rows(:,1),max(pos_val_rows(:,3),0))/5,...
                        std(pos_val_rows(:,2),max(pos_val_rows(:,3),0))/5,...
                        xvals_sub(y_indx),yvals_sub(x_indx),0,-0.2,max(OD_meas(:))]; %Inital guess parameters
    fitobject=fitnlm(pos_val_rows(:,1:2),pos_val_rows(:,3),modelfun,fit_params_guess,...
        'options',opt,...
        'CoefficientNames',cof_names)
    fitparam=fitobject.Coefficients;
    x{ii} = fitparam{1,1};
    dx{ii} = fitparam{1,2};
    y{ii} = fitparam{2,1};
    dy{ii} = fitparam{2,2};
    amp(ii) = fitparam{7,1};
    damp(ii) = fitparam{7,2};
    offset(ii) = fitparam{6,1};
    doffset(ii) = fitparam{6,2};
    hold on
    %
    %subplot(1,4,3)

    fit_pos_val_rows=pos_val_rows;
    %fit_pos_val_rows(:,3)=modelfun(fit_params_guess,pos_val_rows(:,1:2));
    fit_pos_val_rows(:,3)=predict(fitobject,pos_val_rows(:,1:2));
    fit_pos_vals=reshape(fit_pos_val_rows,size(pos_val));
    %imagesc(fit_pos_vals(:,:,3))
    sh=surf(fit_pos_vals(:,:,1)*xy_factor,fit_pos_vals(:,:,2)*xy_factor,fit_pos_vals(:,:,3),...
        'FaceAlpha',0.8,'FaceColor',[100,200,53]./255,'edgecolor','none');

stfig('OD slice');
hold on
OD_avg = mean(OD_meas(:,abs(xvals_sub*xy_factor-550)<50),2);
plot(yvals_sub*xy_factor,OD_avg)
xlabel('z ($\mathrm{\mu}$m)')
ylabel('optical density')




stfig('difference 4');
clf
set(gca, 'FontName', font_name)
set(gca, 'FontSize', font_size)


sh=surf(xvals_sub*xy_factor,yvals_sub*xy_factor,total_image{1}-total_image{2},...
    'FaceAlpha',0.75);
ifh=gca;
shading interp
sh.EdgeColor='k';
%     caxis([0 0.8])
%pbaspect([1,1,1])
% xlim(x_lim*xy_factor)
% ylim(y_lim*xy_factor)
ifh.View= [7.8   14.04];
box on
xlabel('y ($\mathrm{\mu}$m)')
ylabel('z ($\mathrm{\mu}$m)')
zlabel('Intensity (arb. u.)')
shading interp
stfig('integrated profile diff');
% clf

E = @(b,x) 1./b(2).*besselj(1,b(1).*(x-b(3)))./(x-b(3))./b(1);
fit_params_guess = [1,0.06,-200];
fitobject=fitnlm(xvals_sub*xy_factor,trapz(yvals_sub,total_image{1}-total_image{2}),E,fit_params_guess);
fitparam_E=fitobject.Coefficients;

subplot(2,1,1)
hold on
plot(xvals_sub*xy_factor,trapz(yvals_sub,total_image{1}-total_image{2}));
xlabel('y ($\mathrm{\mu}$m)')
ylabel('Integrated Intensity (arb. u.)')
subplot(2,1,2)
hold on
dif_z = trapz(xvals_sub,total_image{1}-total_image{2},2);
plot(yvals_sub*xy_factor,dif_z);
xlabel('z ($\mathrm{\mu}$m)')
ylabel('Integrated Intensity (arb. u.)')
[v1,id1]=max(dif_z(yvals_sub>0));
[v2,id2]=max(dif_z(yvals_sub<0));

(yvals_sub(id1+128)-yvals_sub(id2))*xy_factor;

[X,Y]= meshgrid(xvals_sub,yvals_sub);
X_w = X.*OD_meas;
Y_w = Y.*OD_meas;
X_OD = nansum(OD_meas-fitparam{6,1},2)./size(yvals_sub,1);
Y_OD = nansum(OD_meas-fitparam{6,1})./size(xvals_sub,2);
stfig('integrated OD profile');
subplot(2,1,1)
hold on
plot(xvals_sub*xy_factor,Y_OD);
xlabel('y ($\mathrm{\mu}$m)')
ylabel('Integrated Intensity OD')
subplot(2,1,2)
hold on
plot(yvals_sub*xy_factor,X_OD);
xlabel('z ($\mathrm{\mu}$m)')
ylabel('Integrated Intensity OD')

stfig('integrated OD')
hold on
plot(yvals_sub,trapz(xvals_sub,OD_actual,2))
OD_meas(isnan(OD_meas)) = 0;
num(ii) = 1/sigma.*trapz(yvals_sub,trapz(xvals_sub,OD_actual-fitparam{6,1},2));
OD_int_adj(ii) = trapz(yvals_sub,trapz(xvals_sub,OD_actual-fitparam{6,1},2));
OD_int(ii) = trapz(yvals_sub,trapz(xvals_sub,OD_actual,2));
fit_int(ii) = fitparam{7,1}*(2*pi*fitparam{1,1}*fitparam{2,1});
fit_int_unc(ii) = fit_int(ii).*sqrt((fitparam{7,2}/fitparam{7,1}).^2+(fitparam{1,2}/fitparam{1,1}).^2+(fitparam{2,2}/fitparam{2,1}).^2);
int_amp(ii) = max(trapz(xvals_sub,OD_actual,2));
end
end
%%
stfig('spot size vs distance');
clf
errorbar(dist,cell2mat(x),cell2mat(dx),'kx')
hold on
errorbar(dist,cell2mat(y),cell2mat(dy),'bx')

modelfun = @(beta,dist) (beta(1)*sqrt(1+((dist.*1e-3-beta(2))./(pi*beta(1).^2/(wavelen))).^2));
initials = [15e-5 100e-3];
new_coeff = nlinfit(dist.*1e-3,2.*abs(cell2mat(x)),modelfun,initials);
dist_f = linspace(-180,max(dist),1e4).';
new_x = modelfun(new_coeff,dist_f.*1e-3);
hold on
plot(dist_f,new_x,'LineWidth',1.8)

% end

%%
stfig('amp vs detuning');
clf
b =[6.953052447686723e+01,8.104204261854445e+00];
freq = (dist.*b(2)+b(1));
scatter(freq,amp)
amp_func = @(b,x) b(1)./(1+4.*((x(:,1)-b(2))./b(3)).^2+b(4));
amp_guesses = [0.8,126.,10.2,0];
fitobjecta=fitnlm(freq.',amp.',amp_func,amp_guesses)
xt=linspace(min(freq),max(freq));
amp_fit=predict(fitobjecta,xt.');
hold on
plot(xt,amp_fit)
%%
stfig('amp vs power');
clf
power = dist;
errorbar(power,amp-offset,damp,'kx')
amp_func = @(b,x) b(1)./(1+(x-b(3)).*b(2));
amp_guesses = [0.85,11,5];
fitobjecta=fitnlm(power(2:end),amp(2:end).',amp_func,amp_guesses)
xt=linspace(min(power),max(power));
amp_fit=predict(fitobjecta,xt.');
hold on
plot(xt,amp_fit)