clear   
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%		The business cycle, J.-O. Hairault, F. Portier 	%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%  Choix du Pays %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cntrymap = containers.Map([1,2],{'France','USA'}); % Map the selected number with country name (To display country name later)
prompt= 'Type the number of the country you want to simulate and press Enter: \n 1:France[default] \n 2:USA \n'; % Ask the user to choose a country

% Stock the selected country and check that the user has entered a legit value :
cntryid = input(prompt);
if isempty(cntryid) % test if empty
    cntryid = 1;
elseif (cntryid~=1) || (cntryid~=2) % test if != 1 or 2
    cntryid = 1;
end
country = cntrymap(cntryid); % Stock the country name
fprintf('The model will be simulated for %s\n ',country);


%**********************
% Calibration and steady state
%**********************

C_us = [ 0.95 ; 0.95 ; 0.009  ; 0.008  ;  1.014 ; 0.377 ; 0.009   ; 0.2  ; 0.3   ;  0.025 ; - 1/3 ; 1/9 ; 0.197 ; 0.1  ; 1 ; 0.01  ; 0.988]
C_fr = [ 0.98 ; 0.99 ;  0.009 ; 0.0058 ; 1.028  ;  0.50 ; 0.0098  ; 0.2  ; 0.335 ; 0.0125 ; - 1/3 ; 1/9 ; 0.225 ; 0.25 ; 1 ; 0.007 ; 0.988]


if cntryid == 1 
	C = [C_fr];
else
	 C = [C_us]; 
end	

rho_an = C(1)
rho_ap = C(2)
epsilon_an = C(3)
epsilon_ap = C(4)
g_barre = C(5)
rho_g = C(6)
epsilon_g = C(7)
H = C(8)
alpha = C(9)
delta = C(10)	
a = C(11)
sigma = C(12)
gamma = C(13)
MPY = C(14) 
phi_us = C(15) 
MC = C(16)
beta = C(17)

A_barre = 0.2;
P = 0.3; 

z = (1/beta) - 1 + delta; 

ksy = (alpha / 1 + gamma ) * (1/z) ; 
csy = 1 - delta * ksy  ; 
ysh = (A_barre)^(1/1-sigma) * ksy^(alpha/1-alpha) * H ;

y = ysh * H;
c = y - delta * ksy * y ;
k = ksy * y ; 

w = (1-alpha/ 1+gamma)* ysh ;
gamma_2 = - (1- g_barre * (1+z - delta)) * (MPY * (1/csy))^(1/sigma) ;


vsp = gamma / (1 + gamma) ; 
v = vsp * P ;

MP = MPY * y ;

disp(' ');
param1=[alpha beta delta H rho_ap rho_g gamma sigma];
param2=[y     c     k     w   v   z   MP ];
disp('    alpha     beta  	delta 	   H   	   rho_a    rhog     gamma      sigma'); 
disp(' ');
disp(param1);
disp(' ');
disp('    y         c   	k  	  w         v  	      z	       MP ');
disp(param2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%		   	Linearized Model in Matrix Form 				  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


M1=zeros(7,7); 
M2=zeros(7,5);
M3=zeros(7,2);


M4 = zeros(5,5); 
M5 = zeros(5,7); 
M6 = zeros(2,1); 


