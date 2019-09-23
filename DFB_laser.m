%%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%%%%!!!!!!!!!!!!!!!!!!j!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Code develloped  by   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  BOUCHENE Mohammed Mehdi   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Laboratoire des Tï¿½lï¿½communications (LT), Universitï¿½ 8 mai 1945 Guelma   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Please cite as in your publications if it helps your research.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   This code  allow to simulate      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     Light-Current of  different structures of        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    bulk or, or QW , or MQW   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    DFB laser......      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 


%%%%%:::::::: any want who wants to investigate the code must be familiar
%%%%%with :  TDTW approach,, how to solve time-dependant coupled wave equations using Split-step method
%%%%%::::::::  and to undrestand some physics behind DFB laser 

%%%%%:::::: The code takes  about one hour using INTEL® CORE™ i5-2 2.6 GHz
%%%%% PROCESSOR and 4 GB RAM under Windows 10

%%                                 intialize all parameters to zero on Matlab Buffer

clear all
% clc
% close all
 format long e
%% Define laser Material and structure parameters

c=3e10;                                                                    % velocity of light in free space
q=1.602e-19;                                                               % electron charge
h0=6.626e-34;                                                              % planck 
wl=1.55e-4;                                                                % reference lasing wavelength
f0=c/wl;
gratingperiod=242.1875e-7;                                                 % Grating periodï¿½ï¿½cmï¿½ï¿½ 2*3.2*242.1875e-3um=1.55um
L=600e-4;                                                                  %Laser lengthï¿½ï¿½cmï¿½ï¿½
d=0.2e-4;                                                                  %Thickness (cm)
w=1.5e-4;                                                                  % Width (cm)
V=L*d*w;                                                                   % active region  Volume (cm^3)
gama=0.3;                                                                  % Confinement factor
neff=3.2;                                                                  %Effective index without injection
ng=3.6;                                                                    % Group index
vg=c/ng;                                                                   % group velocity
a=10;                                                                      % Internal loss(cm^-1)
gN=2.5e-16;                                                                % Differential gain(cm^2)
NO=1e18;                                                                   % Transparent carrier density(cm^-3)
e=5e-17;                                                                   % Non-linear gain saturation coefficient(cm^3)
am=4;                                                                      % Linewidth enhancement factor
A=0.1e9;                                                                   % Lineear recombination cofficient(s^-1)
B=1e-10;                                                                   % Bimolecular radiation coefficient(cm^3*s^-1)
C=7.5e-29;                                                                 % Auger coefficient(cm^-6*s^-1)
b=5e-5;                                                                    % Spontaneous emission coefficient
K=1;                                                                       % Transverse Peterman factor
L1=10;L2=10;L3=10;L4=10;
%% Left and right facets parameters

% Left facet parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%left facet power reflectivity
left = 0.97;                                                               % Amplitude of reflectivity
%left = 0.32;
%left = 0.1;
                                                              
Rleft =sqrt(left);
phil=0; 

rl = Rleft*exp(i*phil);                                                     %field reflectivity at left facet with phase;

% Right facet parameters  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%AR facet
right=0.1;    
%right=0.32;

                                                            
Rright=sqrt(right);                                                        % Amplitude of reflectivity
phir=0; 
rr = Rright*exp(i*phir);
%% Computing parameters 

M=60;                                                                      %Number of sections
% M=L;

T=10e-9;                                                                   %Total simulation time
deltz=L/M;                                                                 % space step size
deltt=deltz/vg;
tnum=round(T/deltt);                                                       % total time step 

sinto=ones(1,M);
costo=ones(1,M); 


 
 
 kapa=30*ones(1,M);                                                        % Coupling coefficient                                  
 %kapa=60*ones(1,M);                                                       % Coupling coefficient                                  
 %kapa=80*ones(1,M);                                                       % Coupling coefficient                                  


%% set initial field and power and carrier density

F=zeros(1,M+1); Ft=zeros(1,tnum);
 R=zeros(1,M+1); Rt=zeros(1,tnum); 
S=zeros(1,M);  
S1=zeros(3,M);
N=ones(1,M)*NO;  
PO=zeros(1,tnum);POO=zeros(1,tnum);
lamda_t=zeros(1,tnum);
sn1=1;



%% % MAIN program
D=0:1:150;                                                               % Drive current        
qq=length(D);
I=zeros(qq,1);
power=zeros(qq,1);
for hh=1:qq
for tn=(1:tnum)                                                            %building up fields over tnum time steps;	
 S1(1,:)=S((1:M));
 G=0.5*(gama*gN*(N-NO)./(1+e*S)-a);
 Detu=2*pi*(neff-am*gama*gN*(N-NO)*wl./(4*pi*(1+e*S)))/wl-pi/gratingperiod;

 Spn=sqrt(abs(b*gama*K*B*N.^2/(deltz*vg)));
 a1=randn(1,M);  b1=2*pi*rand(1,M);
 Spn1=sn1*Spn.*a1.*exp(1i*b1);                                             %set up for random forward and backward spontaneous input;
 
 
 % MAIN GUTS FOR UPDATING FIELDS and adding in spontaneous emission;

 % % %First step% % %
 temp_f=exp((G(1:M)-1i*Detu(1:M))*deltz).*F(1:M)+deltz* Spn1(1:M);
 temp_r=exp((G(1:M)-1i*Detu(1:M))*deltz).*R(2:M+1)+deltz*Spn1(1:M);
  F(1)=rl*R(1);    R(M+1)=rr*F(M+1);
  
 % % %Second step% % %

F(2:M+1)=sech(kapa*deltz).*temp_f.*costo+1i*tanh(kapa*deltz).*temp_r.*costo;

 R(1:M)=1i*tanh(kapa*deltz).*temp_f.*sinto+sech(kapa*deltz).*temp_r.*sinto;

 % % %photon density% % %
 S(1:M)=0.5*((abs(F(1:M))).^2+(abs(R(1:M))).^2)+0.5*((abs(F(2:M+1))).^2+(abs(R(2:M+1))).^2);

 S1(2,:)=S(1:M);
 S1(3,:)=0.5*(S1(1,:)+S1(2,:));

 % % %carrier distribution% % %
 I(hh)=D(hh)*0.001;
 N=N+(I(hh)/(q*V)-A*N-B*N.^2-C*N.^3-vg*gN*(N-NO).*S./(1+e*S))*deltt;
 
 % OUTPUT RELATIONSHIPS
 PO(tn)=(abs(F(M+1)))^2*V*vg*h0*c*(1-Rright^2)/(wl*L*gama);                % right hand side's output power
 Ft(tn)=F(M+1);
 
 
 POO(tn)=(abs(R(M+1)))^2*V*vg*h0*c*(1-Rleft^2)/(wl*L*gama);                % Left hand side's output power
 Rt(tn)=R(M+1);
 
 
 lamda_t(tn)=2*gratingperiod*(-wl*am*gama*gN*((N(end)-NO))/4/pi);
end

for xx=(1:tnum-1)
PO(xx)=0.5*(PO(xx)+PO(xx+1));
end

ORight=PO(1:tnum-1)*1000;

power(hh)=sum(PO(5:tnum-1)*1000)/length(PO(5:tnum-1)*1000);

end



xlswrite('power.xlsx',power)
xlswrite('drive.xlsx',D)

 plot(D,power);                                
xlabel('Current(mA)');ylabel('Output power (mW)')

