clc,clear,close all
% Parameters
m = 1; % kg
c = 0.01; % Ns/m
k = 0.2; % N/m
alpha = 1; % N/m^3

% Time
dt = 1e-1;
t = 0:dt:1000;

options = odeset;
options.RelTol = 1e-8;
options.AbsTol = 1e-9;

v0 = 0.1;
[~,y1] = ode45(@(t,y) sys(t,y,m,c,k,alpha),t,[0;v0],options);
v0 = 1;
[~,y2] = ode45(@(t,y) sys(t,y,m,c,k,alpha),t,[0;v0],options);

% Compute Wavelet Transforms
[freq1,mods1] = MoDAL.WaveletTransform(t,y1(:,1),0,0.25);
[freq2,mods2] = MoDAL.WaveletTransform(t,y2(:,1),0,0.25);

figure
tiledlayout('flow');
nexttile;
plot(t,y1(:,1),'k')
xlabel('Time [\cdot]')
ylabel('Displacement [\cdot]')
title('v_0 = 0.1')
nexttile;
imagesc(t,freq1,mods1')
set(gca,'ydir','normal')
ylabel('Frequency [\cdot]') 
xlabel('Time [\cdot]')

figure
tiledlayout('flow');
nexttile;
plot(t,y2(:,1),'k')
xlabel('Time [\cdot]')
ylabel('Displacement [\cdot]')
title('v_0 = 1')
nexttile;
imagesc(t,freq2,mods2')
set(gca,'ydir','normal')
ylabel('Frequency [\cdot]') 
xlabel('Time [\cdot]')

function dy = sys(t,y,m,c,k,alpha)
dy(1,1) = y(2);
dy(2,1) = -1/m*(c*y(2)+k*y(1)+alpha*y(1).^3);
end

