%**************************************************************************%
       %  Studying two terminal Graphene nanoribbon based FET with defect %

      % Matlab code for Bandgap and Transmission calculation for increasing %
                      % number of vacancies along length %
                      
                     

%**************************************************************************%

clear all; clc;

% Constants
q = 1.6e-19; % electron charge C
hbar = 1.06e-34; % h/2pi
h = 2*pi*hbar;% planck's constant
m = 9.1e-31;% mass of electron
zplus=i*1e-3;
qh=q/hbar;BB=0;a=1e-9;
t = -3;% t is tight binding parameter
f1 =1; f2 = 0;% fermilevel
fprintf(1,'Calculating suband and Transmission of defected AGNR:\n');
fprintf(1,'Dimension of the GNR:\n');


% Dimension of GNR
NW =50; % number of Unit Cells in the width
fprintf('Width in unit cells is %d\n', NW);
NL = 5; % number of Unit Cells in the length 
fprintf('Width in unit cells is %d\n', NL);
NU = 4; % number of atoms in an unit cell

NT = NL*NW*NU; % size of hamiltonian(H)

% Convergence criteria
emax = 1e-2;

%% NEGF calculation 

% Hamiltonians with defect

% alpha, beta matrices

alphau = t*diag(ones(1,NU-1),-1)+t*diag(ones(1,NU-1),1);
betaw = zeros(NU);
betal = zeros(NU);
% Condition for Armchair Structure for GNR
betaw(2,1) = t; betaw(3,4) = t;
betal(4,1) = t;
alpha = kron(diag(ones(1,NW)),alphau)+kron(diag(ones(1,NW-1),1),betaw)+kron(diag(ones(1,NW-1),-1),betaw');
%beta = kron(diag(ones(1,NW)),betal);
beta = kron(diag(exp(i*qh*BB*a*a*[1:1:NW])),betal);
H = kron(diag(ones(1,NL)),alpha)+kron(diag(ones(1,NL-1),1),beta)+kron(diag(ones(1,NL-1),-1),beta')

% Incorporating defects
%{
NV=input("Enter number of Vacancy:\n");
fprintf('Total number of vacancies are %d\n',NV);
M=NW*NU;
J=zeros(M,M);
B=zeros(M,M);


NoL=input("Enter number of Vacancy:\n");
m=1;
for m=1:NoL
    db1=input('db1');
    p1=input('p1');
    p2=input('p2');
    p3=input('p3');
    p4=input('p4');
for s= 1:NV
    r=rem(s,2);
   
A1 = zeros(NL);
DU1= input('Enter position of first Defected Unit Cell in the system: ');
A1(DU1,DU1)=1;
B1=zeros(NW*NU);
DB1=input('Enter position of first Defected basis in the system: ');

% In alpha making all the neighbouring elements defected site to zero
   r=rem(DB1,2);
   if r==1
   B1(DB1,DB1-1)=3;
   B1(DB1-1,DB1)=3;    
   B1(DB1,DB1+1)=3;
   B1(DB1+1,DB1)=3;
   
   elseif r==0
   B1(DB1,DB1-1)=3;
   B1(DB1-1,DB1)=3;
   
   end
   % Defective alpha
   Halpha1=kron(A1,B1);
   
 % In beta making all the neighbouring elements defected site to zero

 ni1 = input('Enter position of neighbouring of first Defected basis in the system: ');
 Dbetaw1=zeros(NU);
 
 Dbetaw1(DB1-db1,ni1)=3;
 C1=zeros(NW);
   if r==0
    C1(p1,p2)=1
    
   elseif r==1
    C1(p3,p4)=1 
       
   end
 J1=kron(C1,Dbetaw1);
 R1=zeros(NL);
 R1(DU1,DU1)=1;
 
 % Defective beta
 Hbeta1=kron(R1,J1); 
 
% Combining all defective positions 
 J=J+J1;
 
 B=B+B1;

% Defective hamiltoninan

H = H+Halpha1+Hbeta1+Hbeta1';
   
end
end
%}
%H = kron(diag(ones(1,NL)),alpha)+kron(diag(ones(1,NL-1),1),beta)+kron(diag(ones(1,NL-1),-1),beta');%+Halpha1+Hbeta1+Hbeta1';
%}
% Energy grid
E = linspace(-10,10,501);

%Define the matrices Tr(transmission),gam1,gam2
Tr = zeros(1,length(E));
gam1 = inv(E(1)*eye(NW*NU)-alpha);
gam2 = inv(E(1)*eye(NW*NU)-alpha);
Es = zeros(NT);
Esin = zeros(NT);

%% NEGF Calculations

% Energy loop
for z = 1:length(E)
    
% Surface Green Function calculation 
 e = 100;
 
 while(e > emax)
 gam1new = inv((E(z)+zplus)*eye(NW*NU)-alpha-beta'*gam1*beta);
 e = (sum(sum(abs(gam1new-gam1))))/(sum(sum(abs(gam1new+gam1))));
 gam1 = (gam1+gam1new)/2;
 end
 
 sig1 = beta'*gam1*beta;
 
 e=100;
 while(e > emax  )
 gam2new = inv((E(z)+zplus)*eye(NW*NU)-alpha-beta*gam2*beta');
 e = (sum(sum(abs(gam2new-gam2))))/(sum(sum(abs(gam2new+gam2))));
 gam2 = (gam2+gam2new)/2;
 end
 
 sig2 = beta*gam2*beta';
 
% Self energy matrices
 SelfE1 = kron(diag([1 zeros(1,NL-1)]),sig1);
 SelfE2 = kron(diag([zeros(1,NL-1) 1]),sig2);
 
% Broadening, Gr1 and Gr2
 Gr1 = (1i)*(SelfE1-SelfE1'); Gr2 = (1i)*(SelfE2-SelfE2');

% including the phase breaking processes, D
 d = 1e-12; % phase breaking potential
 e = 100;
 while(e > emax)
 G = inv((E(z)+zplus)*eye(NT)-H-SelfE1-SelfE2-Es);
 Esnew = d*G;
 e = sum(sum(abs(Esnew-Es)))/sum(sum(abs(Esnew+Es)));
 Es = (Es+Esnew)/2;
 end
 
e = 100;
while(e > emax)
    Gn = f1*(G*Gr1*G')+f2*(G*Gr2*G')+(G*Esin*G');
    Esinnew = d*Gn;
    e = sum(sum(abs(Esinnew-Esin)))/sum(sum(abs(Esinnew+Esin)));
    Esin = Esinnew;
end

 Tr(z) = real(trace(Gr1*G*Gr2*G'));

end
% Position of maximum transmission
[max,idx]= max(Tr)

ka = linspace(-pi,pi,501);
eigE = zeros(length(ka),NW*NU);
%Subband calculation
for z = 1:length(ka)
    
 % when having defect
 %[V,D] = eig(alpha+B + (beta+J)*exp(i*ka(z))+ (beta'+J')*exp(-i*ka(z)));
 
  %When having no defect
  [V,D] = eig(alpha + (beta)*exp(i*ka(z))+ (beta')*exp(-i*ka(z)));
 
 eigE(z,:) = diag(D);
 
end
fprintf('Maximum transmission %d', max);
%% Plot Figures

%Figure 1: Transmission plot

hold on
figure(1);
subplot(1,2,1);
h=plot(Tr,E,'r');
xlabel('Transmission(E)');
ylabel('Energy(eV)');
y = [E;Tr];
fileID = fopen('3n+3_pristine.txt','w');
fprintf(fileID,'%f %f\n',y);
fclose(fileID);
title('Energy V Transmission');
%xlim([0 max(Tr)+0.5]);
%ylim([min(E) max(E)]);


% Figure 2: Subband plot


subplot(1,2,2);
plot(ka/pi,eigE,'b');
xlim([-2 1])
ylim([-2 2])
xlabel('ka/\pi');
ylabel('Energy[eV]');
title('Sub-bands');
%ylim([min(E) max(E)]);
