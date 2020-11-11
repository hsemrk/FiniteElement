clear
clc
function E = Ee(x1,x2,n) %retourner n valeur de E sur entre x1 et x2
  h=x2-x1;
  x=x1:h/n:x2;
  for i=1:n+1
    E(i)=200e9; %expression de E
  end
endfunction

function Igz = Igze(x1,x2,n) %retourner n valeur de Igz sur entre x1 et x2
  h=x2 - x1;
  x=x1:h/n:x2;
  for i=1:n+1
    Igz(i)=29e-6; %expression de Igz
  end
endfunction

function q = qe(x1,x2,n) %retourner n valeur de q sur entre x1 et x2
  h=x2-x1;
  x=x1:h/n:x2;
  for i=1:n+1
    q(i)=2400*x(i); %expression de q
  end
endfunction

function Ke = Kele(x1,x2)
  he= x2 - x1; %longeur de l'element 
  n=1000; %nombre de valeurs des fonctions (precision) 
  x=0:he/n:he; % x bar
  %calculer les dérivées secondes des fonctions de forme d'Hermite P3
  for i=1:n+1
    psi_sec(1,i)=-6/(he*he)+12*x(i)/(he*he*he);
    psi_sec(2,i)=-4/he+6*x(i)/(he*he);
    psi_sec(3,i)=6/(he*he)-12*x(i)/(he*he*he);
    psi_sec(4,i)=-2/he+6*x(i)/(he*he);
  end
  E=Ee(x1,x2,n); %obtenir les valeurs de E
  Igz=Igze(x1,x2,n); %obtenir les valeurs de Igz
  %calculer les élements de la matrice Ke
  for i=1:4
    for j=i:4
      Ke(i,j)=trapz(x, E.*Igz.*psi_sec(i,:).*psi_sec(j,:));
      %utiliser le faite que Ke est symetrique pour reduire le nbr d'itérations 
      Ke(j,i)=Ke(i,j); 
    end
  end
endfunction

function Fe = Fele(x1,x2)
  he= x2 - x1; %longeur de l'element 
  n=1000; %nombre de valeurs des fonctions (precision) 
  x=0:he/n:he; % x bar
  %calculer les fonctions de forme d'Hermite P3
  for i=1:n+1
    psi(1,i)=1-3*x(i).*x(i)./(he*he)+2*x(i).*x(i).*x(i)./(he*he*he);
    psi(2,i)=x(i)-2*x(i).*x(i)./he+x(i).*x(i).*x(i)./(he*he);
    psi(3,i)=3*x(i).*x(i)./(he*he)-2*x(i).*x(i).*x(i)./(he*he*he);
    psi(4,i)=-x(i).*x(i)./he+x(i).*x(i).*x(i)./(he*he);
  end
  q=qe(x1,x2,n); %obtenir les valeurs de q
  %calculer les élements du vecteur Fe
  for i=1:4
    Fe(i,1)=trapz(x,q.*psi(i,:));
  end
endfunction

function K = assK(maill) %Assemblage de la matrice de rigidité
  nh=length(maill)-1; %nbt d'elements
  K=zeros(2*nh+2,2*nh+2);
  for i=1:nh 
    K(2*i-1:2*i+2,2*i-1:2*i+2)=K(2*i-1:2*i+2,2*i-1:2*i+2)+Kele(maill(i),maill(i+1));
  end
endfunction

function F = assF(maill) %Assemblage du vecteur force
  nh=length(maill)-1; %nbt d'elements
  F=zeros(2*nh+2,1);
  for i=1:nh
    F(2*i-1:2*i+2,1)=F(2*i-1:2*i+2,1)+Fele(maill(i),maill(i+1));
  end
endfunction

%---Résolution----%
n=10; %nbr d'elements finis
L=1; %longeur de la poutre
f0= 60000; %force appliquée à l'extremité de la poutre
%preparation des vecteur U et Q et application des condition au limites
U=zeros(2*n+2,1);
Q=zeros(2*n+2,1);
Q(2*n+1,1)=f0;

mail=2; %selectionner le maillage
if mail==1 
  m=0:1/n:L; % maillage 1 
else
  j=1:n+1;
  m= L/2 * (1 - cos(pi*(j-1)/n)); % maillage 2
endif
%préparation de la matrice de rigidité et du vecteur force
K=assK(m);
F=assF(m);
%resolution du systeme [K]*U=F+Q
U(3:2*n+2,1)=K(3:2*n+2,3:2*n+2)\(F(3:2*n+2,1)+Q(3:2*n+2,1));
Q(1:2)=K(1:2,1:4)*U(1:4)-F(1:2);
%Regroupement de la déformée
for i=1:n+1
  Vef(i)=U(2*i-1);
end
for i=1:n+1
  Van(i)=f0*m(i)^2 * (3*L-m(i)) /(6*200e9*29e-6) + 2400* L^4 * (20*m(i)^2/L^2-10*m(i)^3/L^3+m(i)^5/L^5)/(120*200e9*29e-6);
end
%traçage
figure (1);
plot(m,Vef,'x',m,Van,'r')
title ("Déformé de la poutre");
xlabel("x(m)");ylabel("v(x)[m]");
legend ("élements finis", "Solution analytique","location","west");
%construction du vecteur Qe des variables secondaires
Qe=zeros(4*n,1);
Qe(1:2)=Q(1:2);
for i=2:n
  Qe(4*i-3:4*i) = Kele(m(i),m(i+1))*U(2*i-1:2*i+2) - Fele(m(i),m(i+1));
end
Qe(3)=-Qe(5);Qe(4)=-Qe(6);
%Regroupement des efforts
Tef(1)=-Qe(1);
for i=2:n+1
  Tef(i)=Qe(4*i-5);
end
%solution analytique
for i=1:n+1
  Tan(i)=f0+2400*(1-m(i)^2/L^2)/2;
end
%traçage
figure (2);
plot(m,Tef,'x',m,Tan,'r')
title ("Effort tranchant");
xlabel("x(m)");ylabel("Ty(x)[N]");
legend ("élements finis", "Solution analytique","location","west");
%Regroupement des moments
Mef(1)=-Qe(2);
for i=2:n+1
  Mef(i)=Qe(4*i-4);
end
%solution analytique
for i=1:n+1
  Man(i)=f0*(L-m(i))+2400*L^2*(2-3*m(i)/L+m(i)^3/L^3)/6;
end
%traçage
figure (3);
plot(m,Mef,'x',m,Man,'r')
title ("Moment fléchissant");
xlabel("x(m)");ylabel("Mfz(x)[N.m]");
legend ("élements finis", "Solution analytique","location","west");
