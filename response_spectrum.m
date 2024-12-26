clear all;
clc;

% Matlab Code to create the response spectrum for the given earthquake
% response

% Read Earthquake Data

[filename,pathname] = uigetfile({''},'Pick an excel file of ground acceleration');
fullname=fullfile(pathname,filename);
eq_data=xlsread(fullname);

% Mass of the system

m=input("Enter the mass of the system (in kg):");

% Damping of the system

zai = input("Enter damping of the system (in percentage):");
zeta=0.01*zai;

% eq_data = xlsread("earthquake.xlsx");
t=eq_data(:,1);
ugg=eq_data(:,2);

% Scale of Acceleration in terms of g

scale = input('Enter scale of the acceleration data in terms of <g>: ');

ug=ugg./scale;

% Plot between SA vs Frequency or SA vs Time period

choice = input('Plot needed: [1: SA vs Frequency, 2: SA vs Time Period ] Choose (1 or 2): ');


a=zeros(size(ug));
v=zeros(size(ug));
u=zeros(size(ug));

u(1)=input("Enter initial displacement: ");
v(1)=input("Enter initial velocity: ");


Tn=(0.01:0.01:6)';
w=(2*pi)./Tn;
f=1./Tn;

a_max=zeros(size(w));

P=zeros(length(ug));
P(1)=-m*ug(1);
for j=1:1:length(w)
    wn=w(j);
    
    % Stiffness of the system
    k=m*(wn)^2;
    dt=t(2)-t(1);
    wd=wn*(sqrt(1-zeta^2));

    a(1)=-(2*m*zeta*wn*v(1)+k*u(1))/m;

    A = exp(-zeta*wn*dt)*((zeta/sqrt(1-zeta^2))*sin(wd*dt)+cos(wd*dt));
    B = exp(-zeta*wn*dt)*((1/wd)*sin(wd*dt));
    C = (1/k)*(((2*zeta)/(wn*dt))+exp(-zeta*wn*dt)*((((1-2*zeta^2)/(wd*dt))-(zeta/sqrt(1-zeta^2)))*sin(wd*dt)-((1+((2*zeta)/(wn*dt))))*cos(wd*dt)));
    D = (1/k)*(1-((2*zeta)/(wn*dt))+(exp(-zeta*wn*dt))*(((2*zeta^2-1)/(wd*dt))*sin(wd*dt)+((2*zeta)/(wn*dt))*cos(wd*dt)));
    
    A_prime = -exp(-zeta*wn*dt) * (wn/sqrt(1-zeta^2)*sin(wd*dt));
    B_prime = exp(-zeta*wn*dt) * (cos(wd*dt) - zeta/sqrt(1-zeta^2)*sin(wd*dt));
    C_prime = (1/k) * (-1/dt + exp(-zeta*wn*dt) * ((wn/sqrt(1-zeta^2) + zeta/(dt*sqrt(1-zeta^2)))*sin(wd*dt) + (1/dt)*cos(wd*dt)));
    D_prime = (1/(k*dt)) * (1 - exp(-zeta*wn*dt) * (zeta/sqrt(1-zeta^2)*sin(wd*dt) + cos(wd*dt)));
    
    %P=zeros(size(ug));

    for i=2:length(ug)
        
        P(i)=-m*(ug(i));
        u(i)=A*u(i-1)+B*v(i-1)+C*P(i-1)+D*P(i);
        v(i)=A_prime*u(i-1)+B_prime*v(i-1)+C_prime*P(i-1)+D_prime*P(i);

        a(i)=-(2*m*zeta*wn*v(i)+k*u(i))/m;
        
    end
    a_max(j)=abs(max(a));
 
end

if choice == 1
    
    figure;
    semilogx(f,a_max);
    title('Response Spectrum (SA vs Frequncy');
    ylabel('SA');
    xlabel('log(Frequency in Hz)');

elseif choice == 2
    
    
    figure;
    semilogx(Tn,a_max);
    title('Response Spectrum (SA vs Time Period');
    ylabel('SA');
    xlabel('log(Time Period');
end