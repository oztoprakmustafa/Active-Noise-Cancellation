
clear all;

load('TF.mat');
load('SEC18R.mat');
load('SEC13R.mat');
T=5000; 
L=50;

% white noise signal,
x_noise=randn(1,T); 
y_d=filter(S_z,S_p, x_noise);
%  Off-line training
Shx=zeros(1,L);     % the state of Sh(z)
Shz=zeros(1,L);     % the weight of Sh(z)
e_z=zeros(1,T);   % data buffer for the identification error
%Applying least mean square algorithm
 
mu=0.001;                      % learning rate                       
for k=1:T                     % discrete time k
    Shx=[x_noise(k) Shx(1:L-1)];  % update the state
    Shy=Shx*Shz';	        % calculate output of Sh(z)
    e_z(k)=y_d(k)-Shy;    % calculate error         
    Shz=Shz+mu*Shx*e_z(k);   % adjust the weight
end

% Lets check the result

X=audioread('WashingMachine-16-8-mono-1000secs.wav');

% and measure the arriving noise at the sensor position,
d=X;
% Initiate the system,
Xh=zeros(size(X));
Cyx=zeros(1,L);
Z = 5000;
xn = zeros(1,L);
y=zeros(1,L);
W = zeros(1,L);
e = zeros(1,Z);
 X = filter(Shz,1,X) ;  %xhat
mu1 =0.00000006;         %step size
%Real ANC part
for k=2:Z                      
    % update the state of xhat
   Xh(k)=e(k-1)+sum(Cyx.*Shz);
   xn = [Xh(k), xn(1:L-1)];
   y=xn*W'; % calculate output of yhat
e(k) = d(k)-y;% calculate error
W =W +mu1*xn*e(k); % adjust the weight

end

figure
plot(e)
title('e')
 figure
 freqz(d(1:1000),1)
 hold on
 freqz(e(1:1000),1)
 lines =findall(gcf,'type','line');
set(lines(1),'color','red');
  set(lines(2),'color','blue');
  legend('e','d')
 figure
plot(d(1:T),'y')
hold on
plot(e, 'r')
%%% Zoom in
title('ANC')
