%mios analysis
%@author Subhrokoli Ghosh
%close all;clc;
clearvars  -except wave1;
%warning off;
tic
%% Loading of data file [load_MicPyVoltage.m should be on the same folder]
if (exist('wave1')) ~= 1
    wave1 = load_MicPyVoltage(); 
end
data = wave1.data';
sample_rate = wave1.sampleRate;
name = string(erase(wave1.filename,".MicPyVoltage"));
date= '16-03-2022';
%% Default values %%%
% everything is calculated in units of µm 
handles.max_time = length(data)/sample_rate;
handles.temp = 300;      
handles.sample_rate = sample_rate;          
handles.time_step= 1/sample_rate;
handles.wavepoints = length(data);
handles.pot_bins = 60;
handles.radius = 0.535/2;     %Adjust radius [µm]
handles.rad = handles.radius;%*1e-9;
handles.viscosity = 1.85e-14*exp(4209/handles.temp+0.0452*handles.temp -3.376e-5*handles.temp^2)*1e-6; % Kg*s/µm 
handles.kbt = (1.38066e-23)*(handles.temp)*1e12;
handles.visc_drag = 6*pi*handles.viscosity*handles.rad;
handles.kappa = []; 
handles.g = [];
handles.title = ['X','Y','Z'];
handles.colors = ['r', 'g', 'b'];
TimeDevider=1;
Decadebin = 0;

%% Data initialization
n=round(length(data)/TimeDevider);newD = zeros(n,3);                       % For scalled data 
fminCutoff = 65e3;fmaxCutoff = 100e3;CutPowerRatio = 50;     % For noise removal qaroud 75kHZ
fgalvo1=[0.2e4,0.2e4,1.0e4];fgalvo2=[2e4,2e4,1.8e4];        % For galvo peak removal
beta= [-0.12,-0.17,0.5];
fp=zeros(3,handles.pot_bins-1);                         % For force from the potential
window_size=handles.wavepoints/10;                                        % For calculating moving average
AveData = zeros(n,3);StdData = zeros(n,3);
LinFitfactor =0.1; ExpFitfactor=2; PlotFactor = 5;
ratio_pot = zeros(3,1);ratio_acfL = zeros(3,1);ratio_psd = zeros(3,1);
tdata= transpose(linspace(0,handles.max_time,handles.wavepoints)); %old definiton by subhro
tdata=tdata(1:n);

%% Dynamics
%% parameters for PSD fitting
fmin=[1,1,10];
fmax = [inf,inf,inf];
power=17;   %adjust this, when taking data samples of different sizes!! See binningPSD


for i=1:3
    newD(:,i) = (data(:,i) - mean(data(:,i)))*1e6;        % Mean removal and rescaling in µm
    %newD(:,i) = allfun.RmvNoise1(newD(:,i),sample_rate,fminCutoff,fmaxCutoff,CutPowerRatio); % Noise removal
    if i==4
        newD(:,i) = allfun.lowfreqfilt(newD(:,i),1000,10,n);
        
    end
    ACFIT(i) = allfun.DoACFit(newD(:,i),sample_rate,LinFitfactor,ExpFitfactor,handles.kbt);      % Autocorrelation and fitting
    POT(i) = allfun.Potential(newD(:,i),ACFIT(i).kappa_exp,handles.pot_bins);      % Potential and fitting with conversion factor
    
    %% PSD
    PSD(i) =allfun.DoPSD(newD(:,i),sample_rate);
    PSDmean(i)=allfun.logbin(PSD(i).freq,PSD(i).psdata,300);
    %PSD1(i) = allfun.GalvoNrmv(newD(:,i),sample_rate,fgalvo1(i),fgalvo2(i),power,beta(i));
    w_cAC=ACFIT(i).kappa_exp/(2*pi*handles.visc_drag);
    ip = [ 0.1,(w_cAC*2*pi)^2, 1e-12];
    PSDFIT(i) = allfun.lorFit(PSDmean(i).FREQ,PSDmean(i).PSD,(PSDmean(i).nbin).^2, ip,[fmin(i), fmax(i)]);
    %PSDFIT1(i) = allfun.psdfit(PSD1(i).frequency,PSD1(i).binPSD,fmin(i),fmax(i),ip,handles.visc_drag);
    fprintf('Stiffness from PSD data %s is : %2.2f pN/µm  \n',handles.title(i), PSDFIT(i).CornerFreq*handles.visc_drag*1e6);
    handles.kappa(i) = ACFIT(i).kappa_exp;

    % Convertion of data units
    handles.g(i) = POT(i).ConvFactor;
    %newD(:,i) = newD(:,i)./handles.g(i);        % Data conversion from Volt to Distance unit
    handles.kappa(i) = ACFIT(i).kappa_exp;
    %% Force from potential
    f_dummy= -diff(POT(i).POT);
    disp_dummy = diff(POT(i).fitdata_xc);
    fp(i,:)= f_dummy./disp_dummy; 
end
%%
for i=1:3     
    %% Moving mean
    AveData(:,i) = movmean(newD(:,i),window_size);
    StdData(:,i) =movstd(newD(:,i),window_size); 
    
    %% Ratio of stiffnessess
    ratio_pot(i)= round(ACFIT(i).kappa_exp/ACFIT(3).kappa_exp);
    ratio_acf(i)=round(ACFIT(i).kappa_linear/ACFIT(3).kappa_linear);
    ratio_psd1(i)= round(PSDFIT(i).CornerFreq/PSDFIT(3).CornerFreq);
    %ratio_psd2(i)= round(PSDFIT1(i).k_psd/PSDFIT1(3).k_psd);
    %% Cross-correlation
    
    j=mod(i,3)+1; k=mod(i+1,3)+1;

    CC(i) = DoCC(newD(:,j),newD(:,k),sample_rate,power); 
end

%% Data resize for plot
dview=1e4;
dgap= handles.wavepoints/dview;
for i=1:3
    dataview(:,i)=newD(1:dgap:end,i);
end
tview=tdata(1:dgap:end);

%% infotable

%info.kappa=wave1.kappa_theo*1e6;
%info.Particles=wave1.NumParticles;
%info.sampleRate=wave1.sampleRate;
%info.w_c_1=wave1.w_c_theo;
%info.gamma_theo=wave1.gamma_theo;

%disp(info)
%Info=struct2table(info);
%disp(Info)

%% ------------------------Plotting---------------------------------
%-------------------------------------------------------------------

     
%% Figure-1
lim1= -max(newD,[],'all')*1.2;lim2= max(newD,[],'all')*1.2;
figure('Name','Total analysis Simulation'); 
set(gcf,'Position',[0 0 1600 900]);
hold on;

for i=1:3
    %Timeseries of data
    subplot(3,4,1+4*(i-1));                        
    allfun.plot_timeseries(tview,tdata,dataview,AveData,lim1,lim2,i)

    %trajectory
    subplot(3,4,2+4*(i-1))
    allfun.plot_2Dhist(newD,lim1,lim2,i);

    %PSD Data
    subplot(3,4,3+4*(i-1))
    allfun.plot_PSD(PSDmean,PSDFIT,handles.kbt,i);

end


sp9 = subplot(3,4,4);hold on;                 % Potential
scatter(POT(1).fitdata_xc,POT(1).POT./handles.kbt,handles.colors(1),'linewidth',1);grid on;
plot(POT(1).fitdata_xc,POT(1).fitdata_POT./handles.kbt,handles.colors(1),'linewidth',1);
scatter(POT(2).fitdata_xc,POT(2).POT./handles.kbt,handles.colors(2),'linewidth',1);
plot(POT(2).fitdata_xc,POT(2).fitdata_POT./handles.kbt,handles.colors(2),'linewidth',1);
scatter(POT(3).fitdata_xc,POT(3).POT./handles.kbt,handles.colors(3),'linewidth',1);
plot(POT(3).fitdata_xc,POT(3).fitdata_POT./handles.kbt,handles.colors(3),'linewidth',1);
%xticks(-5e-7:5e-8:5e-7);
xt = get(gca, 'XTick');set(gca, 'XTick',xt, 'XTickLabel',xt*1e9);
ytickformat('%2.1f');
title('Potential');
xlabel('Displacement (nm)');ylabel('Potential (k_BT)');
str_x= sprintf('k_x: %2.2f pN/µm \nk_y: %2.2f pN/µm \nk_z:% 2.2f pN/µm',...
    ACFIT(1).kappa_exp*1e6, ACFIT(2).kappa_exp*1e6, ACFIT(3).kappa_exp*1e6);
dim = [.25 .1 .15 .15];
h=annotation('textbox',dim,'String',str_x,'Position',sp9.Position,'Hori','center','Vert','middle','FitBoxToText','on');hold off;
h.FontWeight = 'bold';

sp12= subplot(3,4,8);hold on;      % Force
flim2 = ceil(max(max(fp))*1e12)/1e12;
flim1 = -flim2;%floor(min(min(fp))*1e12)/1e12;
scatter(POT(1).fitdata_xc(2:end-1),fp(1,1:end-1),handles.colors(1));
scatter(POT(2).fitdata_xc(2:end-1),fp(2,1:end-1),handles.colors(2));
scatter(POT(3).fitdata_xc(2:end-1),fp(3,1:end-1),handles.colors(3));hold off;
%xticks(-5e-7:5e-8:5e-7);
yticks(flim1:round((flim2 - flim1)*1e12)/10e12 : flim2);
xt = get(gca, 'XTick');set(gca, 'XTick',xt, 'XTickLabel',xt*1e9);
yt = get(gca, 'YTick');set(gca, 'YTick',yt, 'YTickLabel',yt*1e12);
xlabel('Displacement (nm)');ylabel('Force (pN)');grid on;axis tight;grid on;
legend(handles.title(1),handles.title(2),handles.title(3));
title('Force from potential');

sp15= subplot(3,4,12);      % Cross-correlation FFTs
%ccxlmin=min([1/ACFIT(1).tau_linear,1/ACFIT(2).tau_linear,1/ACFIT(3).tau_linear]);
plot(CC(1).frequency,1e0.*CC(1).binPSD,'r-','LineWidth',0.1);set(gca,'xscale','log');set(gca,'yscale','log');hold on;
plot(CC(2).frequency,1e-2.*CC(2).binPSD,'g-','LineWidth',0.1);set(gca,'xscale','log');set(gca,'yscale','log');
plot(CC(3).frequency,1e-4.*CC(3).binPSD,'b-','LineWidth',0.1);set(gca,'xscale','log');set(gca,'yscale','log');hold off;
xlabel('Frequency (Hz)');ylabel('PSD-CC');grid on;
xlim([min([1/ACFIT(1).tau_linear,1/ACFIT(2).tau_linear,1/ACFIT(3).tau_linear]) CC(3).frequency(end)]);
legend(handles.title(1),handles.title(2),handles.title(3));
title('Cross-correlation FFTs');

str = sprintf('Dataname = '+name+'; Date='+date+'; \nRadius = %2.2f nm; Sampling rate = %2.1f kHz; Temperature = %2.1f K; \nKappa ratios--> PSD_1 -> k_x:k_y:k_z = %2d:%2d:%2d;  Exp-ACF -> k_x:k_y:k_z = %2d:%2d:%2d; Lin-ACF --> k_x:k_y:k_z = %2d:%2d:%2d',...
    handles.radius, sample_rate/1e3, handles.temp, ratio_psd1(1),ratio_psd1(2),ratio_psd1(3), ratio_pot(1),ratio_pot(2),ratio_pot(3),ratio_acf(1),ratio_acf(2),ratio_acf(3));
if isfield(wave1,'w_c_theo')
    str2= "; w_c_x:w_c_y:w_c_z:=" + sprintf(" %2.2f; " ,wave1.w_c_theo) + "Hz";
    str = str +str2;
end
sgtitle(str,'FontWeight','bold','FontSize',12);

%% Correlationsplot
figure('Name','Total analysis Simulation Correlation');
set(gcf,'Position',[0 0 1400 600]);

for i = 1:3
    mask(:,i)= ACFIT(i).lags < ACFIT(i).ACFTime/sample_rate*ExpFitfactor;
end

t=[0.001,0.001,0.005];
for i = 1:3
    
    allfun.doubleExpAC(ACFIT(i).lags,ACFIT(i).ACdata,t(i),i,handles.kbt,handles.visc_drag)
    %allfun.plot_AC(ACFIT,mask,i);
    subplot(3,2,2*i);
    allfun.plot_cc(CC,mask,i)

end

%% Special Auto-y
modus1='moving';
modus2='still';
color='g';
figure;
sp8 = subplot(1,3,1);hold on;                     % Semilog
mask2= ACFIT(2).lags  < 1e-3;
plot(ACFIT(2).lags(mask2),ACFIT(2).ACdata(mask2),color,'linewidth',0.5);
plot(ACFIT(2).lags(mask2),ACFIT(2).lags(mask2).*ACFIT(2).LinFitCoff(2)+ACFIT(2).LinFitCoff(1),color,'LineStyle','--');
%plot(ACFIT(3).lags(mask2), ACFIT(2).ExpFitData(mask2),'b--','linewidth',1.5);
hold off;grid on;axis tight;set(gca,'yscale','log');
%legend(modus1,modus2);title('Autocorrelation Y');
dim = [.71 .07 .15 .15];
str_x= sprintf('κ_{yl}: %2.2f pN/µm , still ',...
     ACFIT(2).kappa_linear*1e6);
h=annotation('textbox',dim,'String',str_x,'Position',sp8.Position,'Hori','left','Vert','bottom','FitBoxToText','on');hold off;


xt = get(gca, 'XTick');set(gca, 'XTick',xt, 'XTickLabel',xt*1e3);
xlabel('Delay (ms)');ylabel('Autocorrelation');




sp8 = subplot(1,3,2);hold on;                     % Semilog
mask2= ACFIT(2).lags  < 10e-3;
plot(ACFIT(2).lags(mask2),ACFIT(2).ACdata(mask2),color,'linewidth',0.5);
%plot(ACFIT(3).lags(mask2), ACFIT(2).ExpFitData(mask2),'b--','linewidth',1.5);
hold off;grid on;axis tight;set(gca,'yscale','log');
%legend(modus1,modus2);title('Autocorrelation Y');

xt = get(gca, 'XTick');set(gca, 'XTick',xt, 'XTickLabel',xt*1e3);
xlabel('Delay (ms)');ylabel('Autocorrelation');
      

sp8 = subplot(1,3,3);hold on;                     % linear
mask2= ACFIT(2).lags <  1000e-3;
plot(ACFIT(2).lags(mask2),ACFIT(2).ACdata(mask2),color,'linewidth',0.5);
%plot(ACFIT(3).lags(mask2), ACFIT(2).ExpFitData(mask2),'b--','linewidth',1.5);
hold off;grid on;axis tight;
%legend(modus1,modus2);title('Autocorrelation Y');

xlabel('Delay (s)');ylabel('Autocorrelation');

