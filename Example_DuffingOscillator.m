


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
v0 = 0.5;
[~,y2] = ode45(@(t,y) sys(t,y,m,c,k,alpha),t,[0;v0],options);

% Compute Wavelet Transforms
[freq1,mods1] = WaveletTransform(t,y1(:,1),0,0.25,"numFreq",200);
[freq2,mods2] = WaveletTransform(t,y2(:,1),0,0.25,"numFreq",200);

figure
subplot(2,1,1)
plot(t,y1(:,1),'k')
xlabel('Time [\cdot]'); ylabel('Displacement [\cdot]')
title('v0 = 0.1'); set(gca,'fontsize',14)

subplot(2,1,2)
imagesc(t,freq1,mods1'); set(gca,'ydir','normal')
xlabel('Time [\cdot]'); ylabel('Frequency [\cdot]')
title('v0 = 0.1'); set(gca,'fontsize',14)
colormap(flipud(gray))

figure
subplot(2,1,1)
plot(t,y2(:,1),'k')
ylabel('Displacement [\cdot]'); xlabel('Time [\cdot]')
title('v0 = 0.5'); set(gca,'fontsize',14)

subplot(2,1,2)
imagesc(t,freq2,mods2'); set(gca,'ydir','normal')
xlabel('Time [\cdot]'); ylabel('Frequency [\cdot]')
title('v0 = 0.5'); set(gca,'fontsize',14)
colormap(flipud(gray))

function dy = sys(t,y,m,c,k,alpha)
dy(1,1) = y(2);
dy(2,1) = -1/m*(c*y(2)+k*y(1)+alpha*y(1).^3);
end



