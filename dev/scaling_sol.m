% syms u(t) v(t)
% ode1 = diff(u,2) == 1/(u^3*v);
% ode2 = diff(v,2) == 1/(u^2*v^2);
% odes = [ode1; ode2]
% cond1 = u(0) == 1;
% cond2 = v(0) == 1;
% conds = [cond1; cond2];
% [uSol(t),vSol(t)] = dsolve(odes,conds)

[t,y] = ode45(@myODE,[0 500000],[1; 1;0;0]);

figure(1)
clf
plot(t,y(:,3),'-x',t,y(:,4),'-o')
hold on
% plot(t,sqrt(1+t.^2))
xlabel('Time t');
ylabel('Solution y');
legend('y1','y2')

figure(2)
clf
plot(t,y(:,1).^2.*y(:,2)./t.^3,'-o')
xlim([100,20000])
hold on
% plot(t,sqrt(1+t.^2))
xlabel('Time t');
ylabel('Solution y');
legend('y1','y2')

figure(3)
clf
plot(t,y(:,1)./y(:,2).*1/(5))
xlim([10,500000])
hold on
% plot(t,sqrt(1+t.^2))
xlabel('Time t');
ylabel('Solution y');
legend('y1','y2')

function dy = myODE(t,y)
eta = 1/5;
  dy(1,1) = y(3,1);
  dy(2,1) = y(4,1);
  dy(3,1) = 1/(y(1,1)^3*y(2,1));
  dy(4,1) = eta^2/(y(1,1)^2*y(2,1)^2);
end