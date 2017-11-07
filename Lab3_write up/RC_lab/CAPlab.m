clear; close all
%% declare vars
data_struct1 = lvm_import('datafile030.lvm'); %import 1st set of data
data_struct2 = lvm_import('datafile031.lvm'); %import 2nd set of data
data_struct3 = lvm_import('datafile032.lvm'); %import 3rd set of data
vc1 = data_struct1.Segment1.data; 

vc1 = vc1(125:590);
vc2 = data_struct2.Segment1.data;
vc3 = data_struct3.Segment1.data;
dt = 10e-3; %sample period 
n = length(vc1); %number of samples
time = 0:dt:((n-1)*dt); %time array
scale = 1e3;
K = ceil(1+3.322* log10(n)); % sturgis method for number of bins rounded

%% setup analytical curve
v0 = 9.15;           % V ... your measured input voltage
R = 9.86e3;         %  9.86K... your measured resistance
C = 100e-6;       % 10 micoF ... your measured capacitance
tau = R*C;        % sec ... time constant
time_offset = 0;%1.25;  % sec ... time offset to match data
vc_func = @(t, tau) ...
  v0 * ( 1 - exp(-(t-time_offset)./tau) ) .* (t>time_offset);

%% error analysis
  ER =[1:n];
  for j = 1:n
      ER(j) = vc1(j)-vc_func(time(j),tau);
  end
  myMean = mean(ER);
  myVar = var(ER);
  mySigma = std(ER);
  mySkew = skewness(ER);
  mykurt = kurtosis(ER);
  
%% Uncertainty
u_c = C*.01 + .5e-6;
u_r = R*.0005 + .002;
u_vs = v0*.0003 + .003;
theta_vs = 1-exp(-1/(R*C));
theta_R = -v0*exp(-1/(R*C))*(1/(C*R^2));
theta_C = -v0*exp(-1/(R*C))*(1/(R*C^2));
u_vc = sqrt((theta_vs*u_vs)^2 + (theta_R*u_r)^2 + (theta_C*u_c)^2);
%sigmaTau = sqrt((u_c*R)^2 + (u_r*C)^2);

%% plot data with analytical curve

figure;
plot(time*scale, vc1,'.b');
hold on;                  % allows us to keep plotting on the same axes
plot(time*1e3, vc_func(time,tau),'-r', 'LineWidth', 1); % plot the analytic result
grid on;                  % turn the grid on
xlabel('time (ms)');       % label the x-axis
ylabel('v_c (V)');         % label the y-axis
title('voltage vs time');
legend('recorded voltages','theoretical values','Location','southeast');
save2pdf('capCharge',gcf,300);

figure;
[a,b] = hist(ER,K);
F = bar(b,a/n);
hold on;
title('Frequency Dist. of Errors');
xlabel(['x-values','(bin width, K = ',num2str(K),')']);
ylabel('Frequency');
text(-.3,.5,['mean error = ',num2str(round(myMean,3)),'v']);
text(-.3,.45,['variance = ',num2str(round(myVar,3)),'v^2']);
text(-.3,.4,['standard deviation  = ',num2str(round(mySigma,3)),'v']);
save2pdf('freqDist',gcf,300);


lnData = log(v0-vc1);

theoData = -(1/(R*C)).*time + log(v0);
p = polyfit(time,transpose(lnData),1);
pcurve = p(1).*time + p(2);

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
title('linearized \DeltaV for charging capacitor')
text(2.5,0,['\tau = ',num2str(-p(1))]);

subplot(2,1,2)
plot(time*scale, vc1,'.b');
hold on;                 
plot(time*1e3, vc_func(time,tau),'r');
plot(time*1e3, vc_func(time,tau) + u_vc,'--k');
plot(time*1e3, vc_func(time,tau) - u_vc,'--k');
grid on;   
xlim([4000,4750]);
xlabel('time (ms)');       
ylabel('v_c (V)');         
title('voltage vs time');
legend('recorded voltages','theoretical values','uncertainty bounds',...
    'Location','southeast');
save2pdf('quant',gcf,300);

% 
stemplot(ER,time)
ylabel('Error')

figure;
plot(time*scale, vc1,'.b');
hold on;                 
plot(time*1e3, vc_func(time,tau),'r');
plot(time*1e3, vc_func(time,tau) + u_vc,'--k');
plot(time*1e3, vc_func(time,tau) - u_vc,'--k');
grid on;                  
xlabel('time (ms)');       
ylabel('v_c (V)');         
title('voltage vs time');
legend('recorded voltages','theoretical values','uncertainty bounds',...
    'Location','southeast');
save2pdf('capBounds',gcf,300);

figure;
plot(time*scale, vc1,'.b');
hold on;                 
plot(time*1e3, vc_func(time,tau),'r');
plot(time*1e3, vc_func(time,tau) + u_vc,'--k');
plot(time*1e3, vc_func(time,tau) - u_vc,'--k');
grid on;   
xlim([3000,4750]);
xlabel('time (ms)');       
ylabel('v_c (V)');         
title('voltage vs time');
legend('recorded voltages','theoretical values','uncertainty bounds',...
    'Location','southeast');
save2pdf('closeUp',gcf,300);
