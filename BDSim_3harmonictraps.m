  %% Simulate 1, 2 or 3 Particles that diffuse in gaussian potentials, due to brownian Motion. 
  %%The trap of the traps can move lateraly
%close all;
clear; clc;
%%

%% Grids
SamplingRate= 1e5;
time=5;                                                                                                                         ;
NumParticles= 2;
ztrappos=[0,-800e-9,800e-9];
ztrappos=ztrappos(1:NumParticles);                          
BoxLength= 500e-9+sum(abs(ztrappos));
Numpoints=SamplingRate*time;
xpos= zeros(Numpoints,NumParticles);
ypos= zeros(Numpoints,NumParticles);
zpos= zeros(Numpoints,NumParticles);                                           





%%

%% Parameters and constants
dt= 1/SamplingRate;                                     % Time step [s]
t= (0:dt:time-dt)';
kB= 1.3806e-23;                                         % Boltzmann constant [J/K]
T= 300;                                                 % Absolute temperature [K]
NMradius= [535/2,535/2,535/2];                                          % particle radius [nm]
NMradius=NMradius(1:NumParticles);
eta= 1.85e-14*exp(4209/T+0.0452*T-3.376e-5*T^2);        % Coefficient of viscosity of water at room temp [Pa*s]
gamma= 6*pi*eta*NMradius*1e-9;                          % Viscous drag force on the particle [N*s/m]
D= kB*T./gamma;                                         % Diffusion constant [m^2/s]
E=[15,15,15]*kB*T;
E=E(1:NumParticles); 


% Traps
% position


xtrap=zeros(Numpoints,NumParticles);
ytrap=zeros(Numpoints,NumParticles);
ztrap=zeros(Numpoints,NumParticles)+ztrappos;


%add any trajectory of the traps, that you like:

%period=1; %time of trap to osscilate
%A=0.2e-6; %amplitude of oscillation

%circular movement
%xtrap(:,1)=xtrap(:,1)+A*sin(2*pi*t/period);
%ytrap(:,1)=ytrap(:,1)+A/2*cos(2*pi*t/period);

%oscillatory motion
%ytrap(:,1)=ytrap(:,1)+A*cos(2*pi*t/period);

%complicated motion --> sum of 5 sine waves
A=[1 1 1 0.5 0.2]*1e-7;
f=[10 15,20,100,1000];
phi=2*pi*rand(5,1);

disp=sum((cos(2*pi*(t*f)+phi').*A)');

%ytrap(:,2)=ytrap(:,2)+disp';
%figure;
%plot(t(1:20000),disp(1:20000));
%set(gca, 'XTick',xt, 'XTickLabel',xt*1e3);
%yt = get(gca, 'YTick');set(gca, 'YTick',yt, 'YTickLabel',yt*1e6);
%xlabel('Delay (ms)');ylabel('Displacement (µm)');
%spring konstant
%stiff coverslip
k_stick_coverslip=1;
kx=[2e-6, 4e-6*k_stick_coverslip,2e-6]; 
ky=[2e-6, 20e-6*k_stick_coverslip, 4e-6];
kz=[0.5e-6, 20e-6*k_stick_coverslip, 8e-6];

kx=kx(1:NumParticles);
ky=ky(1:NumParticles);
kz=kz(1:NumParticles);

%Gaussian shape:
sigmax=sqrt(E./kx);
sigmay=sqrt(E./ky);
sigmaz=sqrt(E./kz);
 



%% Some example Values from PRL Alex 2005 that are not necessarily used:
%%Diameter [µm] n (refractive) phaseshift kx ky kz [pN/µm]
%%0.22k 1.57 0.31 1.46 2.36 0.27
%%0.69 1.57 0.98 24.2 26.5 4.29
%%0.85 1.57 1.20 29.4 30.4 5.03
%%1.03 1.57 1.46 28.2 25.2 6.20
%%1.66 1.57 2.35 11.0 10.0 3.85
%%0.64 1.43 0.38 9.45 12.6 1.58
%%1.00 1.43 0.59 10.6 10.0 2.91

%% Thermal displacement 
dthermX = sqrt(2*D*dt).*normrnd(0,1,[Numpoints,NumParticles]);
dthermY = sqrt(2*D*dt).*normrnd(0,1,[Numpoints,NumParticles]);
dthermZ = sqrt(2*D*dt).*normrnd(0,1,[Numpoints,NumParticles]);




%% Dynamics

xpos(1,:)=xtrap(1,:);ypos(1,:)=ytrap(1,:);zpos(1,:)=ztrap(1,:);     % initial position 
figure('Name','Particle trajectory');
%%set(gcf,'Position',[10 10 1000 1000]);hold on;



%Dynamics

%Force


for i=1:Numpoints
    xpos(i+1,:)=xpos(i,:)+ dthermX(i,:) - kx.*(xpos(i,:)-xtrap(i,:))*dt./gamma;
    ypos(i+1,:)=ypos(i,:)+ dthermY(i,:) - ky.*(ypos(i,:)-ytrap(i,:))*dt./gamma;
    zpos(i+1,:)=zpos(i,:)+ dthermZ(i,:)- kz.*(zpos(i,:)-ztrap(i,:))*dt./gamma;
    t(i+1)=t(i)+dt;
end

   
%Plotting
Npic=20;    %number of frames of movie

lim1=max([max(abs(xpos)) max(abs(ypos)) max(abs(zpos))])*3; %set limits
for j=1:Npic
    start=int32((j-1)*Numpoints/Npic);
    stop=int32(j*Numpoints/Npic);
    %plotting



    sp1 = subplot(2,2,1);  
    title('X-Y trajectory');
    histogram2(xpos(1:stop,:),ypos(1:stop,:),[100 100],'facecol','flat');
    grid on;axis equal;view(2);colormap jet;
    xlim([-lim1/2 lim1/2]);ylim([-lim1/2 lim1/2]);
    xt = get(gca, 'XTick');set(gca, 'XTick',xt, 'XTickLabel',xt*1e9);
    yt = get(gca, 'YTick');set(gca, 'YTick',yt, 'YTickLabel',yt*1e9);
    xlabel('X (nm)');ylabel('Y (nm)');
    
    
    sp1 = subplot(2,2,2);  
    title('X-Y trajectory');
    histogram2(xpos(1:stop,:),zpos(1:stop,:),[100 100],'facecol','flat');grid on;axis equal;view(2);colormap jet;
    xlim([-lim1/2 lim1/2]);ylim([-lim1/2 lim1/2]);
    xt = get(gca, 'XTick');set(gca, 'XTick',xt, 'XTickLabel',xt*1e9);
    yt = get(gca, 'YTick');set(gca, 'YTick',yt, 'YTickLabel',yt*1e9);
    xlabel('X (nm)');ylabel('Z (nm)');

    sp1 = subplot(2,2,3);  
    title('X-Y trajectory');
    histogram2(ypos(1:stop,:),zpos(1:stop,:),[100 100],'facecol','flat');grid on;axis equal;view(2);colormap jet;
    xlim([-lim1/2 lim1/2]);ylim([-lim1/2 lim1/2]);
    xt = get(gca, 'XTick');set(gca, 'XTick',xt, 'XTickLabel',xt*1e9);
    yt = get(gca, 'YTick');set(gca, 'YTick',yt, 'YTickLabel',yt*1e9);
    xlabel('Y (nm)');ylabel('Z (nm)');%}
    pause(0.1);
    
end




hold off;

%% format data for totalsimulation analysis
% 
%
if NumParticles == 1
    wave1.data=[xpos'; ypos'; zpos'];

    
else
    wave1.data=[sum(xpos'); sum(ypos'); sum(zpos')];
end

wave1.NumParticles=NumParticles;
wave1.sampleRate=SamplingRate;
wave1.w_c_theo=[kx(1:NumParticles) ky(1:NumParticles) kz(1:NumParticles)]./([gamma(1:NumParticles) gamma(1:NumParticles) gamma(1:NumParticles)]);
wave1.kappa_theo=[kx(1:NumParticles) ky(1:NumParticles) kz(1:NumParticles)];
wave1.gamma_theo=gamma(1:NumParticles);
wave1.filename='Simulation of Particles.MicPyVoltage';