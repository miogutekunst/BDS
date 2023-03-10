classdef allfun
    %% @author Subhrokoli Ghosh
    % This is a copilation of all the functions used for Brownian dynamics
    % simulation and data analysis
    % @p : parameter
    
    methods(Static)      
 %%  Function RmvNoise
        function y = RmvNoise(Signal,samplerate,fmin,fmax,CutPowerRatio)
 % Descrition: This function is used for removing noise from PSDs in
 % frequency domain within a frequency band defined by @p fmin  and @p fmax
 % this code is version 1. This can be used to remove laser noise at high
 % frequencies
 % @p Signal: timeseries signal
 % @p samplerate: sampling rate of the signal
 % @p fmin: lower limit of the frequency band to be changed
 % @p fmax: upper limit of the frequency band to be changed
 % @p CutPowerRatio: (input noise amplitude/output noise amplitude)^2
            n = length(Signal);
            dt=1/samplerate;
            tt = dt*(n-1);
            fhat = fft(Signal,n);
            psd1 = fhat.*conj(fhat)/n;

            cutP = psd1(tt*fmin:tt*fmax);
            cutfft = fhat(tt*fmin:tt*fmax);

            ind1 = cutP<sqrt(CutPowerRatio)*cutP(end);
            cutfft = ind1.*cutfft;
            fhat(tt*fmin:tt*fmax)=cutfft;
            y = real(ifft(fhat));
            y=y-mean(y);
        end
  
 %%  Function RmvNoise1
        function y = RmvNoise1(Signal,samplerate,fmin,fmax,CutPowerRatio)
 % Descrition: This function is used for removing noise from PSDs in
 % frequency domain within a frequency band defined by @p fmin  and @p fmax
 % this code is version 1. This can be used to remove laser noise at high
 % frequencies
 % @p Signal: timeseries signal
 % @p samplerate: sampling rate of the signal
 % @p fmin: lower limit of the frequency band to be changed
 % @p fmax: upper limit of the frequency band to be changed
 % @p CutPowerRatio: (input noise amplitude/output noise amplitude)^2
            n = length(Signal);
            dt=1/samplerate;
            tt = dt*(n-1);
            fhat = fft(Signal,n);
            psd1 = fhat.*conj(fhat)/n;

            cutP = psd1(tt*fmin:tt*fmax);
            cutfft = fhat(tt*fmin:tt*fmax);
            aveP = (mean(cutP(1:10)) + mean(cutP(end-10:end)))/2;

            ind1 = cutP<sqrt(CutPowerRatio)*aveP;
            ind2 = cutP>=sqrt(CutPowerRatio)*aveP;
            cutfft = ind1.*cutfft + ind2.*aveP.*cutfft;
            fhat(tt*fmin:tt*fmax)=cutfft;
            y = real(ifft(fhat));
            y=y-mean(y);
        end
  
 %% Function GalvoNrmv
        function y = GalvoNrmv(Signal,samplerate,fmin,fmax,power,beta)
 % Descrition: This function is used for removing noise from PSDs in
 % frequency domain within a frequency band defined by @p fmin  and @p fmax
 % this code is version 1. This can be used to remove noise from Galvo
 % mirrors in setup at B3 Lab at BNP IMTEK
 % @p Signal: timeseries signal
 % @p samplerate: sampling rate of the signal
 % @p fmin: lower limit of the frequency band to be changed
 % @p fmax: upper limit of the frequency band to be changed
 % @p CutPowerRatio: (input noise amplitude/output noise amplitude)^2
            y = struct();
            n=length(Signal);
            s1 = normrnd(0,1/6,[n,1]);
            Decadebin=0;
            psdbin = binningPSD(Signal,power,samplerate);
            psdbin1 = binningPSD(s1,power,samplerate);
            index1 = find(psdbin.frequency>fmin,1,'first');
            index2 = find(psdbin.frequency<fmax,1,'last');
            cutfreq = psdbin.frequency(index1:index2);
            cutpsd = psdbin.binPSD(index1:index2);
            cutpsd1 = psdbin1.binPSD(index1:index2);
            newp = cutpsd1./cutfreq.^(2+beta);
            y.AveRatio = (cutpsd(1)+cutpsd(end))/(newp(1)+newp(end));
            psdbin.binPSD(index1:index2) = newp.*y.AveRatio;
            y.frequency = psdbin.frequency;
            y.binPSD =psdbin.binPSD;
        end
        
 %% Function DoACFitSim
        function OP = DoACFitSim(data,sample_rate,ExpFitfactor,kBT)
 % Description: This function computes autocorrelation function (ACF) of a given
 % input scaled data. It also performs linear and exponential fitting of
 % the data. For linear fitting, it taks the first 30 points of the
 % ACF,lags data. For exponential fitting (2P and 3P), the data lenght is chosen from
 % the @p ExpFitfactor.
 %@p ExpFitfactor: # times autocorrelation length of the data
 % Lin fit fun: a*x+b
 % Exp fit fun (2P): a*exp(b*x)
 % Exp fit fun (3P): a+b*exp(c*x)
            OP = struct();
            n=length(data);
            
            [ACdata,lags] = autocorr(data,n-1);
            OP.lags = lags./sample_rate;
            OP.ACdata = ACdata.*var(data);
            OP.ACFTime = find(OP.ACdata < var(data)*exp(-1),1);

            %Linear fitting
            xldata = OP.lags(1:30);
            yldata = OP.ACdata(1:30);
            OP.LinFitCoff = polyfit(xldata,yldata,1);
            OP.LinFitData = polyval(OP.LinFitCoff,OP.lags);
            OP.kappa_linear = kBT/OP.LinFitCoff(2);
            OP.gamma_linear = -kBT/OP.LinFitCoff(1);
            OP.tau_linear = OP.gamma_linear/OP.kappa_linear;


            %Exponential fitting
            dlength = ceil(ExpFitfactor*OP.ACFTime);
            xedata = OP.lags(1:dlength);
            yedata = OP.ACdata(1:dlength);
            ExpFitFn2P= @(c,x) c(1)*exp(c(2).*x);
            options = optimset('MaxIter',5000,'MaxFunEvals',5000);
            ip2P = [OP.LinFitCoff(2), -1/OP.tau_linear];
            OP.ExpFitCoff = fminsearch(@(c) norm(yedata - ExpFitFn2P(c,xedata)), ip2P, options);
            OP.ExpFitData = ExpFitFn2P(OP.ExpFitCoff,OP.lags);
            OP.tau_exp2P = -1/OP.ExpFitCoff(2);
            OP.kappa_exp2P = kBT/OP.ExpFitCoff(1);
            OP.gamma_exp2P = OP.tau_exp2P*OP.kappa_exp2P;

            ExpFitFn3P= @(c,x) c(1) + c(2)*exp(c(3).*x);
            ip3P = [var(data), OP.LinFitCoff(2), -1/OP.tau_linear];
            OP.ExpFitCoff3P = fminsearch(@(c) norm(yedata - ExpFitFn3P(c,xedata)), ip3P, options);
            OP.ExpFitData3P = ExpFitFn3P(OP.ExpFitCoff3P,OP.lags);
            OP.tau_exp3P = -1/OP.ExpFitCoff3P(3);
            OP.kappa_exp3P = kBT/OP.ExpFitCoff3P(2);
            OP.gamma_exp3P = OP.tau_exp3P*OP.kappa_exp3P;
        end
  
  %% Function DoACFit
        function OP = DoACFit(data,sample_rate,LinFitfactor,ExpFitfactor,kBT)
 % Description: This function computes autocorrelation function (ACF) of a given
 % input unscaled data. It also performs linear and exponential fitting of
 % the data. For linear fitting, it takes the first 30 points of the
 % ACF,lags data. For exponential fitting (2P and 3P), the data lenght is chosen from
 % the @p ExpFitfactor.
 %@p ExpFitfactor: # times autocorrelation length of the data
 % Lin fit fun: a*x+b
 % Exp fit fun (2P): a*exp(b*x)
 % Exp fit fun (3P): a+b*exp(c*x)
            OP = struct();
            n=length(data);
            [OP.NormACdata,OP.lags] = autocorr(data,n-1);
            OP.ACdata=OP.NormACdata*var(data);
            OP.ACFTime = find(OP.ACdata < exp(-1)*var(data),1); 
            
            %if OP.ACFTime > 20e-3*sample_rate % upper limit for ACTime used for fitting
            %    OP.ACFTime = 6e-3*sample_rate;
            %    'ACT at limit';
            %end

            %if max(OP.ACdata(ceil(sample_rate):end)) > 6e-3*sample_rate
            %    OP.ACFTime = 6e-3;
            %end

            
                
            OP.lags = OP.lags./sample_rate;

            %Linear fitting
            dlength = ceil(LinFitfactor*OP.ACFTime);
            
            xdata = OP.lags(1:dlength);
            ydata = OP.ACdata(1:dlength);
            LinFitFn= @(b,x) b(1) + b(2).*x;
            OP.LinFitCoff = fminsearch(@(b) norm(ydata - LinFitFn(b,xdata)), [var(data),-var(data)/OP.ACFTime*sample_rate]);
            tau_lin = find(LinFitFn(OP.LinFitCoff,OP.lags)< LinFitFn(OP.LinFitCoff,xdata(1))*exp(-1),1);
            OP.tau_linear = OP.lags(OP.ACFTime);
            %OP.kappa_linear = gamma/OP.tau_linear;
            OP.kappa_linear = kBT/OP.LinFitCoff(1);
            OP.gamma_linear = -kBT/OP.LinFitCoff(2);
            OP.wc_linear=-OP.LinFitCoff(2)/OP.LinFitCoff(1)/2/pi;
            OP.LinFitData = LinFitFn(OP.LinFitCoff,xdata);

            %Exponential fitting
            dlength1 = ceil(ExpFitfactor*OP.ACFTime);
            xedata = OP.lags(1:dlength1);
            yedata = OP.ACdata(1:dlength1);
            ExpFitFn= @(c,x) c(1)*exp(c(2).*x);
            OP.ExpFitCoff = fminsearch(@(c) norm(yedata - ExpFitFn(c,xedata)), [var(data),-var(data)/OP.ACFTime*sample_rate]);
            
            OP.kappa_exp    = kBT/OP.ExpFitCoff(1);
            OP.wc_exp       =-OP.ExpFitCoff(2);
            OP.tau_exp       =1/OP.wc_exp;
            OP.gamma_exp    =OP.kappa_exp/OP.wc_exp;
            OP.ExpFitData   = ExpFitFn(OP.ExpFitCoff,OP.lags);
            
        end
    
        %% Function DoACFit2
        function OP = DoACFit2(data,sample_rate,LinFitfactor,ExpFitfactor,kBT)
 % Description: This function computes autocorrelation function (ACF) of a given
 % input unscaled data. It also performs linear and exponential fitting of
 % the data. For linear fitting, it takes the first 30 points of the
 % ACF,lags data. For exponential fitting (2P and 3P), the data lenght is chosen from
 % the @p ExpFitfactor.
 %@p ExpFitfactor: # times autocorrelation length of the data
 % Lin fit fun: a*x+b
 % Exp fit fun (2P): a*exp(b*x)
 % Exp fit fun (3P): a+b*exp(c*x)
            OP = struct();
            n=length(data);
            [OP.NormACdata,OP.lags] = autocorr(data,n-1);
            OP.ACdata=OP.NormACdata*var(data);
            OP.ACFTime = find(OP.ACdata < exp(-1)*var(data),1);
            
            if OP.ACFTime > 20e-3*sample_rate % upper limit for ACTime used for fitting
                OP.ACFTime = 6e-3*sample_rate;
                'ACT at limit';
            end

            if max(OP.ACdata(ceil(sample_rate):end)) > 6e-3*sample_rate
                OP.ACFTime = 6e-3;
            end

            
                
            OP.lags = OP.lags./sample_rate;
            
        
            %Exponential fitting
            dlength1 = ceil(ExpFitfactor*OP.ACFTime);
            xedata = OP.lags(1:dlength1);
            yedata = OP.ACdata(1:dlength1);
            ExpFitFn= @(c,x) c(1)*exp(c(2).*x) + c(3)*exp(c(4).*x) ;
            OP.ExpFitCoff = fminsearch(@(c) norm(yedata - ExpFitFn(c,xedata)), [var(data),-var(data)/OP.ACFTime*sample_rate,var(data)/20,-var(data)/OP.ACFTime*sample_rate*20]);
            
            OP.ExpFitData   = ExpFitFn(OP.ExpFitCoff,OP.lags);
            
            OP.kappa_exp    =kBT/OP.ExpFitCoff(1);
            OP.wc_exp       =-OP.ExpFitCoff(2);
            OP.gamma_exp    =OP.kappa_exp/OP.wc_exp;
            OP.tau_exp      =1/OP.wc_exp;

            OP.kappa_exp2    = kBT/OP.ExpFitCoff(3);
            OP.wc_exp2       =-OP.ExpFitCoff(4);
            OP.gamma_exp2    =OP.kappa_exp2/OP.wc_exp2;
            OP.tau_exp2       =1/OP.wc_exp2;
        end
  %% Function DoACFitDual
        function OP = DoACFitDual(data,sample_rate,ExpFitfactor,KappaRatio)
 % Description: This function computes autocorrelation function (ACF) of a given
 % input scaled data from asymetric confining potential. @p KappaRatio
 % defines the ratio of stifnessess of left-half and right half of the
 % potential. It also performs exponential fitting of the ACF.  
 % The data lenght is chosen from the @p ExpFitfactor.
 % @p ExpFitfactor: # times autocorrelation length of the data
 % Fit fun: a*exp(-b*x)+(a/c)*exp(-d*x)
            OP = struct();
            n=length(data);
            [OP.ACdata,OP.lags] = autocorr(data,n-1);
            OP.ACFTime = find(OP.ACdata < exp(-1),1);
            OP.lags = OP.lags./sample_rate;

            % %Linear fitting
            % dlength = ceil(LinFitfactor*OP.ACFTime);
            % xdata = OP.lags(1:dlength);
            % ydata = OP.ACdata(1:dlength);
            % LinFitFn= @(b,x) b(1) + b(2).*x;
            % OP.LinFitCoff = fminsearch(@(b) norm(ydata - LinFitFn(b,xdata)), ones(2,1));
            % OP.tau_linear = OP.lags(OP.ACFTime);
            % OP.kappa_linear = gamma/OP.tau_linear;
            % OP.LinFitData = LinFitFn(OP.LinFitCoff,OP.lags);

            %dual Exponential fitting
            dlength1 = ceil(ExpFitfactor*OP.ACFTime);
            xedata = OP.lags(1:dlength1);
            yedata = OP.ACdata(1:dlength1);
            FitFn = fittype( 'a*exp(-b*x)+(a/c)*exp(-d*x)', 'independent', 'x', 'dependent', 'y');
            opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
            opts.Lower = [0 0 KappaRatio 0];
            opts.Upper = [1 Inf KappaRatio Inf];
            opts.StartPoint = [1 1e3 KappaRatio 1e3];
            OP.FitCoff = fit(xedata,yedata,FitFn, opts);



            %ExpFitFn2= @(c,x) c(1)*exp(c(2).*x) + c(3)*exp(c(4)*x);
            % 
            % OP.ExpFitCoff = fminsearch(@(c) norm(yedata - ExpFitFn2(c,xedata)),[1,-1e-3,1,-5e-2]);
            % OP.tau1_exp = -1/OP.ExpFitCoff(2);
            % OP.tau2_exp = -1/OP.ExpFitCoff(4);
            % OP.kappa1_exp = OP.tau1_exp*gamma;
            % OP.kappa2_exp = OP.tau2_exp*gamma;
            % OP.ExpFitData = ExpFitFn2(OP.ExpFitCoff,OP.lags);
        end
  
  %% Function DoCC
        function op=DoCC(x,y,sample_rate,Power)
% Description: This fucntion computes the crosscorrelation of two given
% timeseries @p x and y. It also computes Power spectral density of the
% crosscorrelation data and output is binned into 2^power points in
% frequency domain.
% @p x,y : input time series whose crosscorrelation has to be performed
% @p sample_rate: sampling rate of both data. Note: They should be same.
% @p Power: length of binned data: 2^Power
            n=length(x);
            [xc,lags] = xcorr(x,y,n);
            lags=lags(n+1:end-1);
            %cXY = binningPSD(xc,power,sample_rate);
            [freq,psdata]=DoPSD(xc,sample_rate);
            invcXY=abs(ifft(psdata));
            op.lags=lags./sample_rate;
            op.cc=invcXY;

            nFFT = 2^Power;
            n1=length(xc);
            nbins = floor(n1/nFFT);
            dFFT = n1 - nFFT*nbins;
            dbins = floor((n1-nFFT)/dFFT);
            f = zeros(nFFT/2,1);
            pd = zeros(nFFT/2,1);

                for i=1:dbins
                    data = xc((i-1)*dFFT + 1 : (i-1)*dFFT + nFFT);
                    [f_buff,pd_buff] = DoPSD(data,sample_rate);
                    f = f + (1/dbins).*f_buff;
                    pd = pd + (1/dbins).*pd_buff;
                end
            op.frequency = f;
            op.binPSD = pd;
        end

 %% Function Potential       
        function op =Potential(data,kappa,bins)
 % Description: This function computes Potential of a given unscalled data. 
 % It also calibrates the data (unscaled--> scalled). To do this, it fits 
 % the potential with a parabola and determines the stiffness. It compares 
 % the fitted stiffness with input stiffness (which can be a result from 
 % another fittype,eg, AC fitting, PSD fitting) to calculate the conversion 
 % factor g. Naturally, the output stiffness from this fitting is same as
 % kappa from input fittype.
 % @p ConvFactor(g): unscalled data/scalled data
 % @p bins: No. of bins to plot potential (similar to histogram bins)
            kB = 1.38066e-23;
            T = 300;
            d= data - mean(data);
            [counts, centers] = hist(d,bins);
            p1=-log(counts)*(kB*T);
            p1(isinf(p1)|isnan(p1)) = 0;
            IMfitP= polyfit(centers,p1,2);
            op.ConvFactor = sqrt(0.5*kappa/IMfitP(1));
            op.fitdata_xc = centers./op.ConvFactor;
            op.fitP = polyfit(op.fitdata_xc,p1,2);
            op.POT = p1;
            op.fitdata_POT = polyval(op.fitP,op.fitdata_xc);
        end
        
 %% Function Potential_n
        function op =Potential_n(data,bins)
 % Description: This function computes Potential of a given scalled data.
 % Additionally, it fits the potential with a parabola to determine the 
 % stiffness of the potential. The output stiffness from this fitting is 
 % independent of other fittypes.
 % See also Potential
 % @p bins: No. of bins to plot potential (similar to histogram bins)
            kB = 1.38066e-23;
            T = 300;
            d= data - mean(data);
            [counts, centers] = hist(d,bins);
            p1=-log(counts)*(kB*T); 
            p1(isinf(p1)|isnan(p1)) = 0;
            op.POT = p1;
            op.xc = centers;
            op.fitP = polyfit(centers,p1,2);
            op.fitdata_POT = polyval(op.fitP,centers);
            op.kappa = 2*op.fitP(1);
        end
        
%% Function lor     
        function y =lor(x,k)
% Description: This function calculates lorengian of the input data
% @p x= frequency array, 
% @p k(1)= Amplitude of Lor, 
% @p k(2)= corner frequency^2
  
            y = k(1)./(x .^2 + k(2));
        end
        
  %% Function DoPSD
        function op=DoPSD(d,SamplingRate)
% Description: This function computes Power spectral density of a given
% data. The outputs are frequency (@p freq) and PSD (@p psdata) arrays. 
% The length of output arrays are n/2, where n is the length of input data.  
% @p d: input data
% @p SamplingRate: sampling rate of the data
            n=length(d);
            wd = hann(n).*d*2;                %multiply with window function
            tf = (n/SamplingRate);          %tf: duration of data record
            freq1 = (0:n)./tf;              %frequencies
            freq1 = transpose(freq1);
            L = 1:floor(n/2);               %only take first half of data
            op. freq = freq1(L);
            xft = fft(wd)/n;
            psd = xft.*conj(xft);
            op.psdata = psd(L);
            % a_x=abs(xft/n);                 % normalised amp of fourier coefficients of data
            % ax=a_x(1:n/2 +1);                % Selection of positive half of the amp of data
            % ax(2:end-1)=2*ax(2:end-1);   
            % psdata= ax.^2;
        end
 
   %% Function binningPSD
        function op = binningPSD(d,Power,SamplingRate)
% Description: This function computes Power spectral density of a given
% time series with binning. It distributes the input data into n bins, where 
% length of each bin is 2^power. For optimum result, power should be chosen 
% in such a way that the total length of data > 2n or more. It computes the
% PSD of each bin separatesly using DoPSD function and take and average afterwrds. 
% The outputs are binned frequency and binned PSD arrays.The length of output 
% arrays are 2^(power-1).
% @p d: input data
% @p SamplingRate: sampling rate of the data
% See also DoPSD.
            nFFT = 2^Power;
            n1=length(d);
            nbins = floor(n1/nFFT);
            dFFT = n1 - nFFT*nbins;
            dbins = floor((n1-nFFT)/dFFT);
            f = zeros(nFFT/2,1);
            pd = zeros(nFFT/2,1);

                for i=1:dbins
                    data = d((i-1)*dFFT + 1 : (i-1)*dFFT + nFFT);
                    [f_buff,pd_buff] = DoPSD(data,SamplingRate);
                    f = f + (1/dbins).*f_buff;
                    pd = pd + (1/dbins).*pd_buff;
                end
            op.frequency = f;
            op.binPSD = pd;
        end
 
   %% Function psdfit
        function op =psdfit(freq,psdata,fmin,fmax,ip,visc_drag)
% Description:  This function fits a given data with Lorentzian
% functionusing the function lor (See lor).
% @p freq: input frequency array (eg, output of binningPSD)
% @p psdata: input PSD array (eg, output of binningPSD). length should be
% same as freq array.
% @p fmin: Lower cutoff frequency for Lorentzian fitting
% @p fmax: Upper cutoff frequency for Lorentzian fitting
% @p ip: Innitial guess array for Lorentzian fitting, 
%    ip: [Lor amplitude, Lor corner freq]
% @p visc_drag: 6*pi*viscous coff*particle radious
           
            f_d = freq(freq>fmin & freq<fmax);
            pd_d = psdata(freq>fmin & freq<fmax);
            options = optimset('MaxIter',5000,'MaxFunEvals',5000);
            op.fitparam = fminsearch(@(k) norm(pd_d - lor(f_d,k)),ip,options);
            op.CornerFreq = sqrt(op.fitparam(2));
            op.k_psd = 2*pi*visc_drag*sqrt(op.fitparam(2));
        end


  function op = lorFit(F, P, w, ip,f)
%CREATEFIT(F,P)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input: F
%      Y Output: P
%      Weights: F
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 04-Nov-2022 15:53:03


%% Fit: 'untitled fit 1'.
%remove nans
Nans=isnan(P);
P=P(~Nans);
F=F(~Nans);
w=w(~Nans);
[xData, yData, weights] = prepareCurveData( F, P, w );

% Set up fittype and options.
ft = fittype( '2*a./((2*pi*x).^2+b)+c', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.StartPoint = ip;
opts.Display = 'Off';
opts.Weights = weights;
opts.Lower = [0 0 0];
opts.MaxIter = 1000;
excludedPoints = (xData < f(1)) | (xData > f(2));
opts.Exclude = excludedPoints;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
op.fitparam = coeffvalues(fitresult);
op.CornerFreq = sqrt(op.fitparam(2));
op.D=op.fitparam(1);     
op.gof= gof;

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, xData, yData, excludedPoints );
set(gca,'xscale','log');set(gca,'yscale','log');
legend( h, 'P vs. F with w', 'Excluded P vs. F with w', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'F', 'Interpreter', 'none' );
ylabel( 'P', 'Interpreter', 'none' );
grid on

end

        %% Function logbin

        function op =logbin(x,ydata,n)
              binbounds=logspace(log10(x(2)),log10(x(end)),n);
            op.PSD=1:n-1;
            op.FREQ=1:n-1;
            for i=1:n-1
                xbin=(x<binbounds(i+1) & x>=binbounds(i));
                op.PSD(i)=mean(ydata(xbin));
                op.FREQ(i)=mean(x(xbin));
                op.nbin(i)=sum(xbin);
            end
           
    

           
        end   

         %% Function plot timeseries

        function plot_timeseries(tview,tdata,dataview,AveData,lim1,lim2,i)
            color=['r','g','b'];
            name=['X','Y','Z'];
            

            plot(tview,dataview(:,i),color(i));hold on;
            plot(tdata,AveData(:,i),'k-'); hold off;
            ylim([lim1 lim2]);
            yt = get(gca, 'YTick');set(gca, 'YTick',yt, 'YTickLabel',yt);
            xlabel('Time (sec)');ylabel(strcat(name(i),' (??m)'));
            grid on;legend(strcat(name(i),'-signal'));title(strcat(name(i),'-signal'));
           
    

           
        end   
        %% Plot 2Dhistogram
        function plot_2Dhist(newD,lim1,lim2,i)

            
            
            name=['X','Y','Z'];
            j=mod(i,3)+1;
            k=mod(i+1,3)+1;
            if i==2
                k=j;
                j=1;
            end
            histogram2(newD(:,j),newD(:,k),[80 80],'facecol','flat');grid on;axis equal;view(2);colormap jet;
            xlim([lim1 lim2]);ylim([lim1 lim2]);
            xt = get(gca, 'XTick');set(gca, 'XTick',xt, 'XTickLabel',xt);
            yt = get(gca, 'YTick');set(gca, 'YTick',yt, 'YTickLabel',yt);
            xlabel(strcat(name(j),'(??m)'));ylabel(strcat(name(k),'(??m)'));
            title(strcat(name(j),'-',name(k),' trajectory'));
          
        end   

        %% Plot PSD
        function plot_PSD(PSDmean,PSDFIT,kBT,i)
    
            
            hold on;
            color=['r','g','b'];
            name=['X','Y','Z'];
            

                     
            title(strcat(name(i),' PSD'));
            plot(PSDmean(i).FREQ,PSDmean(i).PSD,strcat(color(i),'-'),'LineWidth',0.5);
            set(gca,'xscale','log');set(gca,'yscale','log');
            loglog(PSDmean(i).FREQ,lor(PSDmean(i).FREQ*2*pi,PSDFIT(i).fitparam),'black','LineWidth',0.5); grid on;
            
            
  
            xlabel('Frequency (Hz)');ylabel('PSD (??m^2/Hz)');grid on;

            str1=strcat('k_',name(i));
            
            str2= sprintf(': %2.2f pN/??m ', PSDFIT(i).CornerFreq*kBT/PSDFIT(i).D*1e6 );
            
            str_k=strcat(str1,str2);


            str1=strcat('w_c_',name(i));
            str2= sprintf(': %2.2f Hz ', PSDFIT(i).CornerFreq );

            str_w=strcat(str1,str2);

            str_x = [str_k newline str_w];
            

            h= text(PSDmean(i).FREQ(1)*10,PSDmean(i).PSD(1)/1000,str_x); 
            h.FontWeight = 'bold';
        end   


        %% Plot AC
        function plot_AC(ACFIT,mask,i)
    
            sp=subplot(3,2,1+2*(i-1));

            hold on;
            color=['r','g','b'];
            name=['X','Y','Z'];
            

                     
            title(strcat('Autocorrelation-',name(i)));
       
            
            
            plot(ACFIT(i).lags(mask(:,i)),ACFIT(i).ACdata(mask(:,i)),strcat(color(i),'-'),'linewidth',0.5);
            plot(ACFIT(i).lags(mask(:,i)), ACFIT(i).ExpFitData(mask(:,i)),'k--','linewidth',1.5);
            hold off;grid on;axis tight;set(gca,'yscale','log');
            legend('Auto-Y','Efit-Y');
            
            xt = get(gca, 'XTick');set(gca, 'XTick',xt, 'XTickLabel',xt*1e3);
            xlabel('Delay (ms)');ylabel('Autocorrelation (??m)');
            try 
                str_l= sprintf('??_{l}: %2.2f pN/??m, ??_{l}: %2.2f ms, w_{cl}: %2.1f Hz \n',...
                   ACFIT(i).kappa_linear*1e6, ACFIT(i).tau_linear*1e3,  ACFIT(i).wc_linear);
            catch
                str_l= " ";
            end

            try 
                str_e2= sprintf('??_{e2}: %2.2f pN/??m, ??_{e2}: %2.2f ms, w_{ce2}: %2.1f Hz \n',...
                    ACFIT(i).kappa_exp2*1e6,  ACFIT(i).tau_exp2*1e3,  ACFIT(i).wc_exp2);
            catch
                str_e2= " ";
            end

            str_e= sprintf('??_{e}: %2.2f pN/??m, ??_{e}: %2.2f ms, w_{ce}: %2.1f Hz',...
                    ACFIT(i).kappa_exp*1e6,  ACFIT(i).tau_exp*1e3,  ACFIT(i).wc_exp);
            str_x= str_l + str_e2 + str_e;
            dim = [.71 .07 .15 .15];
            h=annotation('textbox',dim,'String',str_x,'Position',sp.Position,'Hori','left','Vert','bottom','FitBoxToText','on');hold off;
            h.FontWeight = 'bold';
        end   

       function [fitresult, gof] = doubleExpAC(lags, AC, weight_time,i,kBT,gamma)
%CREATEFIT(lags,AC,WEIGHTS)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input: lags
%      Y Output: AC
%      Weights: weights
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 17-Oct-2022 11:23:23


%% Fit: 'untitled fit 1'.
weights=ones(length(lags),1);
[xData, yData, weights] = prepareCurveData( lags, AC,weights);

% Set up fittype and options.

heavy= 1:ceil(weight_time/lags(2)); %number of points that are getting a heigher weight

weights(heavy)=10;
excludedPoints = yData < yData(1)/6;
npoints=ceil(sum(~excludedPoints*1.5));

ft = fittype( 'kBT/k_1*exp(-k_1*x/gamma1)+kBT/k_2*exp(-k_2*x/gamma2)', 'independent', 'x', 'dependent', 'y' );

opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Algorithm = 'Levenberg-Marquardt';
opts.Display = 'Off';
opts.Lower = [0, 0, kBT,0,0];
opts.Upper = [inf, inf, kBT,inf,inf];
opts.MaxIter = 1000;
opts.Robust = 'LAR';
opts.StartPoint = [gamma, gamma, kBT, kBT/yData(1)/2, kBT/yData(1)*2];
opts.Weights = weights;
opts.Exclude = excludedPoints;


opts2 = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts2.Algorithm = 'Levenberg-Marquardt';
opts2.Display = 'Off';
opts2.Lower = [gamma, gamma, kBT,0,0];
opts2.Upper = [gamma, gamma, kBT,inf,inf];
opts2.MaxIter = 1000;
opts2.Robust = 'LAR';
opts2.StartPoint = [gamma, gamma, kBT, kBT/yData(1)/2, kBT/yData(1)*2];
opts2.Weights = weights;
opts2.Exclude = excludedPoints;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
[fitresult2, gof2] = fit( xData, yData, ft, opts2 );

% Plot fit with data.

sp=subplot(3,2,1+2*(i-1));

hold on;
plot(xData(heavy),yData(heavy),Marker="o",LineStyle="none");
plot( fitresult, xData, yData, excludedPoints );
plot(xData,kBT/fitresult2.k_1*exp(-fitresult2.k_1*xData/fitresult2.gamma1)+kBT/fitresult2.k_2*exp(-fitresult2.k_2*xData/fitresult2.gamma2))

xlim([0 lags(npoints)]);
legend('weights*10','AC', 'Excluded AC ', 'unknown gamma a*exp(-b*x)+c*exp(-d*x)','known gamma a*exp(-b*x)+c*exp(-d*x)', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xt = get(gca, 'XTick');set(gca, 'XTick',xt, 'XTickLabel',xt*1e3);
xlabel('Delay (ms)');ylabel('Autocorrelation (??m)');

coeffs = coeffnames(fitresult);
coeffvals= coeffvalues(fitresult);
ci = confint(fitresult,0.95);
str1 = sprintf('\n %s = %0.3f pN/??m','\kappa_1',fitresult.k_1*1e6);
str2 = sprintf('\n %s = %0.3f pN/??m','\kappa_2',fitresult.k_2*1e6);
str3 = sprintf('\n %s = %0.3f *10^-9','\gamma_1',fitresult.gamma1*1e9);
str4 = sprintf('\n %s = %0.3f *10^-9 ','\gamma_2',fitresult.gamma2*1e9);
str5 = sprintf('\n \n %s = %0.3f pN/??m ','\kappa_1',fitresult2.k_1*1e6);
str6 = sprintf('\n %s = %0.3f pN/??m ','\kappa_2',fitresult2.k_2*1e6);

annotation('textbox',[.6 .5 .3 .2],'String',['Coefficients: ', str1, str2,str3, str4,str5,str6],'Position',sp.Position,'Hori','left','Vert','bottom','FitBoxToText','on');

grid on
end

    %% plot_CC
    function plot_cc(CC,mask,i)
        hold on;
        color=['r','g','b'];
        name=['X','Y','Z'];

        j=mod(i,3)+1;
        k=mod(i+1,3)+1;
        plot(CC(i).lags(mask(:,j) | mask(:,k)),CC(i).cc(mask(:,j) | mask(:,k)),color(i),'linewidth',1);
        grid on;axis tight;set(gca,'yscale','log');
        xt = get(gca, 'XTick');set(gca, 'XTick',xt, 'XTickLabel',xt*1e3);
        xlabel('Delay (ms)');
        ylabel(strcat('Cross-corr ',name(j),'-',name(k)));
        title(strcat('Cross-corr ',name(j),'-',name(k)));
    end

%% lowfreqfilter
function filtdata  = lowfreqfilt(data,filtfac,range, n)
    %This FIlter was designed to filter out the noise generated by a slowly oscillation
    %disturbance.
    % It first calculates the PSD. The first datapoints are used
    %to calculate the mean of the plateau in the beginning of the data.
    %Then any values higher by the filtfactor than the plateau are discarded.
    %Instead the mean of the neighbourhoods data (range) is  used.
    %This procedure is repeated once again, as the first time, the
    %disturbances can compromise the mean value ofthe plateau.
    % to get back to a filtered position data, the backtransorm of the FFT is used. 
    % THe PSD is just used to locate the problematic frequency domains,
    % while the FFT is used to go back to the position data.
    %calc fourier
    
        a=data;

        

        FFT=fft(a);
        

        %PSD
        PSD_=abs(FFT(1:ceil(length(FFT)/2))).^2;

        %remove noise
        % filter in fourier domain
        plateau=1:length(PSD_) < ceil(length(PSD_/1e3));
        a_start=mean(PSD_(plateau));
        excludedPoints = PSD_ > a_start*filtfac;
        a_start=mean(PSD_(~excludedPoints & plateau'));
        excludedPoints = PSD_ > a_start*filtfac;
        a_start=mean(PSD_(~excludedPoints & plateau'));
        
        N_exclude=sum(excludedPoints);
        flipped=flip(excludedPoints);
        excludedPointsfromFFT=[excludedPoints;flipped(1:end-1)];
        %input the value of the neighboursdata where the errorness data is
        ind=find(excludedPointsfromFFT);
        for k=1:sum(excludedPointsfromFFT)
            
            minind=ind(k)-range;
            maxind=ind(k)+range;
            if maxind > n
                maxind=n;
            end
            if minind < 1
                minind=1;
            end

            FFT(ind(k))=nanmean(FFT(minind:maxind));
               
        end
        
        %backtransform

        filtdata=real(ifft(FFT))-mean(real(ifft(FFT)));
end

        function [fitresult, gof] = weight_lor_fit(x, y, weights,lim,wc_start)
%CREATEFIT(X,Y,WEIGHTS)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input: x
%      Y Output: y
%      Weights: weights
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 14-Sep-2022 08:15:55


%% Fit: 'Lorentian fit with weigths'.
[xData, yData, weights_1] = prepareCurveData( x, y, weights );

% Set up fittype and options.
ft = fittype( 'a./(x.^2+b.^2)+c', 'independent', 'x', 'dependent', 'y' );
excludedPoints = (xData < lim(1)) | (xData > lim(2));
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
%starting point
D2=y(~excludedPoints);
D2=nanmean(D2(1:7));
opts.StartPoint = [D2*wc_start.^2 wc_start y(end)*0.5];
opts.MaxIter = 4000;
opts.Weights = weights_1;
opts.Exclude = excludedPoints;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, xData, yData, excludedPoints );set(gca,'xscale','log');set(gca,'yscale','log'); 
legend( h, 'y vs. x with weights', 'Excluded y vs. x with weights', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'x', 'Interpreter', 'none' );
ylabel( 'y', 'Interpreter', 'none' );
grid on
        

        end
    end
end

