%% 
clear; close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% RC signal response lab example %%%
%%%%%%%%%% Ryan Haseman %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% Import lvm file from myRIO data capture %%%%%%%
data_struct1 = lvm_import('datafile030.lvm'); %import 1st set of data

% parse capacitor voltage data
vc = data_struct1.Segment1.data; 
% selecting data at t=0 (when circut was closed) to t = 5RC
vc = vc(125:590);

%%%%%%%% Setup time array %%%%%%%%%

dt = 10e-3;         % data was sampled every 1ms 
n = length(vc);     % number of samples
time = 0:dt:((n-1)*dt); % time array
scale = 1e3;

%% Setup theoretical curve
Vs = 9.15;              % input voltage from battery (volts)
R = 9.86e3;             %  measured resistance from resitor (ohms)
C = 100e-6;             % measured capacitance from capacitor (F)
tau = R*C;              % time constant (seconds)

vc_func = @(t) ...
  Vs * ( 1 - exp(-t./tau) );  % theoretical curve

%% error analysis
%%%%%%%%% taking the difference between the observed data and the %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% theoretical curve %%%%%%%%%%%%%%%%%%%%%%%%%%%%

ER =NaN*ones(1,n);        % allocate mem
  for j = 1:n
      ER(j) = vc(j)-vc_func(time(j));    
  end
  
  myMean = mean(ER);     % some stats ...
  myVar = var(ER);
  mySigma = std(ER);

 K = ceil(1+3.322* log10(n)); % sturgis method for # of bins rounded up 
 
%% Uncertainty
%%%%%%% using measurement uncertianty equations to find the uncertianty of
%%%%%%% Vc due to our R, C and Vs measurements with the multimeter

% component uncertainties
u_c = C*.01 + .5e-6;   % uncertianty of capacitance as per Fluke manual
u_r = R*.0005 + .002;  % uncertianty of resistance as per Fluke manual
u_vs = Vs*.0003 + .003;  % uncertianty of DC Vs as per Fluke manual

% Sensitivity coefficients
theta_vs = 1-exp(-1/(R*C));                % Vs coeff
theta_R = -Vs*exp(-1/(R*C))*(1/(C*R^2));   % R coeff
theta_C = -Vs*exp(-1/(R*C))*(1/(R*C^2));   % C coeff

% uncertainty of Vc 
u_vc = sqrt((theta_vs*u_vs)^2 + (theta_R*u_r)^2 + (theta_C*u_c)^2);

%% Plotting
%%%%%%%  Plot of observed Vc vs theoretical curve %%%%%%%%%
figure;
    plot(time*scale, vc,'.b');
    hold on;                  
    plot(time*1e3, vc_func(time),'-r');
    grid on;                  
    xlabel('time (ms)');       
    ylabel('V_c (V)');         
    title('Capacitor Voltage vs time');
    legend('observed V_c','theoretical V_c','Location','southeast');
    save2pdf('capCharge',gcf,300);

%%%%%%%%%% Distribution of Errors  %%%%%%%%%%%%%%
figure; 
    [a,b] = hist(ER,K);
    F = bar(b,a/n);
    hold on;
    title('Distribution of Errors');
    xlabel(['x-values','(bin width, K = ',num2str(K),')']);
    ylabel('Frequency');
    save2pdf('freqDist',gcf,300);

%%%%%%%%%%% Linearizing data %%%%%%%%%%%%%%%

lnData = log(Vs-vc);   % linearizing as deltaV
% linear regression function  (order 1)
p = polyfit(time,transpose(lnData),1);
pcurve = p(1).*time + p(2); % regression line

%%%%%%%%%%%%%%%% linear plot of data with regression line %%%%%%%%%%%%%%
figure; 
    plot(time,-lnData, '.b')
    hold on;
    grid on;
    plot(time,-pcurve,'r');
    xlabel('time (s)')
    ylabel('ln(\DeltaV)')
    legend('observed data','regression line','Location','southeast')
    title('linearized \DeltaV for charging capacitor')
    text(2.5,0,['\tau = ',num2str(-p(1))]);
    save2pdf('linear',gcf,300);

%%%%%%%%%%%%%% plot showing quantization error %%%%%%%%%%%%%%%%%%
figure; 
    subplot(2,1,1)
    plot(time,-lnData, '.b')
    hold on;
    grid on;
    plot(time,-pcurve,'r');
    xlim([4,4.75]);
    xlabel('time (s)')
    ylabel('ln(\DeltaV)')
    legend('observed data','regression line','Location','southeast')
    title('linearized \DeltaV for charging capacitor (close up)')
   
    subplot(2,1,2)
    plot(time*scale, vc,'.b');
    hold on;                 
    plot(time*1e3, vc_func(time),'r');
    grid on;   
    xlim([4000,4750]);
    xlabel('time (ms)');       
    ylabel('v_c (V)');         
    title('Capacitor voltage vs time (close up)');
    legend('recorded voltages','theoretical values',...
        'Location','southeast');
    save2pdf('quant',gcf,300);

%% Error amplitude vs sample 
% stemplot(ER,time)
% ylabel('Error')

%%%%%%%% Original plot with uncertainty bounds %%%%%%%%%%%%%
figure;
    plot(time*scale, vc,'.b');
    hold on;                 
    plot(time*1e3, vc_func(time),'r');
    plot(time*1e3, vc_func(time) + u_vc,'--k'); % upper bound
    plot(time*1e3, vc_func(time) - u_vc,'--k'); % lower bound
    grid on;                  
    xlabel('time (ms)');       
    ylabel('v_c (V)');         
    title('Capacitor voltage vs time');
    legend('recorded voltages','theoretical values',...
        'uncertainty bounds','Location','southeast');
    save2pdf('capBounds',gcf,300);

%%%%%%%%%%% close up section of original plot to show bounds %%%%%%%%%%%%
figure;
    plot(time*scale, vc,'.b');
    hold on;                 
    plot(time*1e3, vc_func(time),'r');
    plot(time*1e3, vc_func(time) + u_vc,'--k');
    plot(time*1e3, vc_func(time) - u_vc,'--k');
    grid on;   
    xlim([3000,4750]);
    xlabel('time (ms)');       
    ylabel('v_c (V)');         
    title('Capacitor Voltage vs time (close up)');
    legend('recorded voltages','theoretical values',...
        'uncertainty bounds','Location','southeast');
    save2pdf('closeUp',gcf,300);
