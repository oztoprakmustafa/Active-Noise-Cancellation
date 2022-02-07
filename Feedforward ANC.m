%%FxLMS Algorithm
%Mustafa Oztoprak
clear all;

load('TF.mat');
load('SEC18R.mat');
load('SEC13R.mat');
T=1000; 
R=500;
num=P_z;
den=P_p;
h=tf(num',den');
num1=S_z;
den1=S_p;
v=tf(num1',den1');
[x ,y] = tfdata(h, 'v');
Pz=[x ,y];
[c,d]=tfdata(v, 'v');
 
Sz=[c,d];
L=100;

mu=0.005;  



% white noise signal,
x_noise=randn(1,T); 
y_d=filter(Sz, 1, x_noise);

%  the identification process
Shx=zeros(1,L);     % the state of Sh(z)
Shz=zeros(1,L);     % the weight of Sh(z)
e_z=zeros(1,R);   % data buffer for the identification error

% and apply least mean square algorithm
                       % learning rate
for k=1:R                     % discrete time k
    Shx=[x_noise(k) Shx(1:L-1)];  % update the state
    Shy=Shx*Shz';	        % calculate output of Sh(z)
    e_z(k)=y_d(k)-Shy;    % calculate error         
    Shz=Shz+mu*e_z(k)*Shx;   % adjust the weight
end

% Lets check the result
figure
  plot(e_z)
  title('e')
 figure
 freqz(Sz,1)
 hold on
 freqz(Shz,1)
 lines =findall(gcf,'type','line');
set(lines(1),'color','red');
  set(lines(2),'color','blue');
  legend('Sz','Shz')
 




X=SEC18R;
% and measure the arriving noise at the sensor position,
d=filter(Pz, 1, X);
% Initiate the system,
xn = zeros(1,L);
y=zeros(1,L);
W = zeros(1,L);
e = zeros(1,R);

 Sx=zeros(1,L);     
                          
for k=1:R                       
     xn = [X(k), xn(1:L-1)];
    y=xn*W';
    Sx=[y Sx(1:L-1)]; 
    yn_hat = Sx.*Shz;
    e(k) = d(k)-yn_hat(k);

     
    xn_hat=filter(Shz,1,xn);
    %Wz =  Wz + mu/norm(xn_hat,2)*err(k)*xn_hat;
    W =  W + mu*e(k)*xn_hat;

end
figure
  plot(e)
  title('e')
 figure
 freqz(d,1)
 hold on
 freqz(e,1)
 lines =findall(gcf,'type','line');
set(lines(1),'color','red');
  set(lines(2),'color','blue');
  legend('d','e')

