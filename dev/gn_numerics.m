% 3d gaussian
lc = 300e-6;
res = 20e-9;

g1 = @(t) exp(-t.^2./(2.*lc^2));
g2 = @(t) 1-abs(g1(t)).^2;
g3 = @(t1,t2) 1-g1(t1).^2-g1(t2).^2-g1(t1-t2).^2+2.*g1(t1).*g1(t2)*g1(t1-t2);

t_max = 3;
num_pts = 1e3;

t1 = rand(1,num_pts).*t_max.*lc;
t2 = rand(1,num_pts).*t_max.*lc;

g1_data = [g1(t1);g1(t2)];
g2_data = [g2(t1);g2(t2)];
g3_data = g3(t1,t2.');

% resolution

t1_res = t1+2.*randn(1,num_pts).*res;
t2_res = t2+2.*randn(1,num_pts).*res;

% saturation

%% binning

figure(121)
clf
subplot(1,2,1)
scatter(t1_res,g1_data(1,:))
subplot(1,2,2)
scatter(t1_res,g2_data(1,:))

figure(122)
clf
scatter3(t1_res,t2_res,g3_data)