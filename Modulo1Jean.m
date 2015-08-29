%% Este arquivo contém a inicialização de todas as variáveis
%% utilizadas nas simulações do sistema de quatro tanques.


% área da base dos tanques

A1 = 28;
A2 = 32;
A3 = 28;
A4 = 32;

% área dos orifícios de saída dos tanques

a1 = 0.071;
a2 = 0.057;
a3 = 0.071;
a4 = 0.057;

% constante gravitacional  

g = 981;

% constante do medidor

kc = 0.5;

% matriz C dos estados para saída

C=[kc 0 0 0;0 kc 0 0];


% constantes das bombas

%k1 = 3.14;
%k2 = 3.29;
k1 = 3.33;
k2 = 3.35;

% valores de divisão dos fluxos das bombas

gama1 = 0.7;
%gama2 = 0.5;
gama2 = 0.6;

% condições iniciais do sistema 

h1t0 = 0;
h2t0 = 0;
h3t0 = 0;
h4t0 = 0;



%2. Como todos os termos são positivos definidos, T1,T2,T3 e T4 são
%positivos definidos, logo isso implica em valores positivos.


%5. a abertura do orificio do tanque 3 não influencia, pois este ja se
%encontra em equilibrio, ou seja, sua vazão é igual a vazão de entrada
%desse mesmo tanque, que assim define indiretamente o equilibrio do tanque
%1.

%--------------------------------1-------------------------------------
h1eq = 10;
h2eq = 10;

A = [gama1*k1/sqrt(a1^2*2*g) (1-gama2)*k2/sqrt(a1^2*2*g);(1-gama1)*k1/sqrt(a2^2*2*g) gama2*k2/sqrt(a2^2*2*g)];

b = [sqrt(h1eq);sqrt(h2eq)];

v = inv(A)*b

h3eq = (1-gama2)^2*k2^2*v(2)^2/(a3^2*2*g)

h4eq = (1-gama1)^2*k1^2*v(1)^2/(a4^2*2*g)


%--------------------------------4---------------------------------------

%vazão de saida do tanque 3:

vouth3 = a3*sqrt(2*g*h3eq)

 %--------------------------------------3--------------------------------
%Como todos os termos são positivos definidos, não há como h3 ser maior que
%h1 e h4 ser maior que h2. pelas proprias equações de equilibrio.
 
 %---------------------------------------6--------------------------------
 T1 = A1/a1*sqrt(2*h1eq/g);
 T2 = A2/a2*sqrt(2*h2eq/g);
 T3 = A1/a3*sqrt(2*h3eq/g);
 T4 = A1/a4*sqrt(2*h4eq/g);
 
 A = [-inv(T1) 0 A3/(A1*T3) 0;0 -inv(T2) 0 A4/(A2*T4);0 0 -inv(T3) 0; 0 0 0 -inv(T4)]
 B = [gama1*k1/A1 0;0 gama2*k2/A2;(1-gama1)*k1/A3 0;0 (1-gama2)*k2/A4]
 C = [kc 0 0 0;0 kc 0 0]
 D = zeros(2,2);
 I = eye(4);

 
% s=sym('s');
% R = [(s*I-A) -B;C D];
% zerotrans = solve(det(C*inv(s*I - A)*B +D),s)
% zeroinv - solve(det(R),s)
 
 

 %---------------------------------------7--------------------------------
 sys = ss(A,B,C,D)
 G = tf(sys)
poles = pole(G)
 zeroinv =zero(sys)
 zerotrans =zero(G)
 
 %--------------------------------------8-------------------------------

C2 = [kc 0 0 0;0 0 0 kc];

sys2 = ss(A,B,C2,D)

zeroinv2 = zero(sys)

G2 = tf(sys2)

zerotrans2 = zero(G2)

poles2 = pole(G2)

autovalores2 = eig(A)

pole(sys)
