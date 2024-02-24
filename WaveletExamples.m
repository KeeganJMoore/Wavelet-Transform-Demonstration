

clc,clear,close all

t = -10:1e-3:10;

gaussianWavelet = -t.*exp(-t.^2./2);

mexicanHatWavelet = (2/(sqrt(3)*pi^0.25))*exp(-t.^2/2).*(1-t.^2);

haarWavelet = 0*t;
haarWavelet(sum(t <= 0):sum(t<1/2)) = 1;
haarWavelet(sum(t <= 1/2):sum(t<1)) = -1;

omega0 = 4; % Central frequency of the Morlet wavelet in rad/s
morletWavelet = 1/(pi^(1/4))*exp(1i*omega0*t).*exp(-t.^2/2);


figure
subplot(2,2,1)
plot(t,gaussianWavelet,'k')
ylim([-1 1])
xlabel('Time [\cdot]'); ylabel('Amplitude [\cdot]')
title('Gaussian Wavelet'); set(gca,'fontsize',14)

subplot(2,2,2)
plot(t,mexicanHatWavelet,'k')
ylim([-1 1])
xlabel('Time [\cdot]'); ylabel('Amplitude [\cdot]')
title('Mexican Hat Wavelet'); set(gca,'fontsize',14)

subplot(2,2,3)
plot(t,haarWavelet,'k')
xlim([-2 2]); ylim([-1.5 1.5])
xlabel('Time [\cdot]'); ylabel('Amplitude [\cdot]')
title('Haar Wavelet'); set(gca,'fontsize',14)

subplot(2,2,4)
plot(t,real(morletWavelet),'k',t,imag(morletWavelet),'r')
ylim([-1 1])
xlabel('Time [\cdot]'); ylabel('Amplitude [\cdot]')
title('Morlet Wavelet'); set(gca,'fontsize',14)
set(gcf,"Renderer","painters")

%% Plot Morlet in Complex Plane
figure
plot3(t, real(morletWavelet), imag(morletWavelet),'k', 'LineWidth',0.5)
hold on
plot3(t, real(morletWavelet), zeros(size(t))-1.5,'r')
plot3(t, zeros(size(t))-2, imag(morletWavelet),'b')
box on
axis([-6  6    -2  2    -1.5  1.5])
view([-125  30])
xlabel('Time [s]', 'Rotation',-30)
ylabel('Real Axis', 'Rotation',10)
zlabel('Imag Axis')
legend('Morlet wavelet, \morletWavelet(\tau)',...
    'Re(\morletWavelet(\tau))','Im(\morletWavelet(\tau))')
set(gca,'fontsize',14)


%% Plot Dilation and Translation

clc,clear,close all

t = linspace(-10,10,200);
a = 1; tau = [0 3 6]';
in1 = (t-tau)/a;
mexicanHatWavelet1 = (2/(sqrt(3)*pi^0.25))*exp(-in1.^2/2).*(1-in1.^2);

a = [1 2 4]'; tau = 0;
in2 = (t-tau)./a;
mexicanHatWavelet2 = (2/(sqrt(3)*pi^0.25))*exp(-in2.^2/2).*(1-in2.^2);


figure
S1 = subplot(2,2,1);
plot(t,mexicanHatWavelet1)
xlabel('Time [\cdot]'); ylabel('Amplitude [\cdot]')
title('Translation'); set(gca,'fontsize',14)
S1.Position(2) = 0.6;
S1.Position(3) = S1.Position(3)*2;
S1.Position(4) = 0.9*S1.Position(4);

S2 = subplot(2,2,3);
plot(t,mexicanHatWavelet2)
xlabel('Time [\cdot]'); ylabel('Amplitude [\cdot]')
title('Dilation (Stretch/Compress)'); set(gca,'fontsize',14)
S2.Position
S2.Position(3) = S2.Position(3)*2;
S2.Position(4) = 0.9*S2.Position(4);


%% 
