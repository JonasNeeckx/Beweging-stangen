%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Kinematica en werkuigendynamica.
%
% Voorbeeldanalyse van een vierstangenmechanisme.
%
% Bram Demeulenaere <bram.demeulenaere@mech.kuleuven.be>
% Maarten De Munck <maarten.demunck@mech.kuleuven.be>
% Johan Rutgeerts <johan.rutgeerts@mech.kuleuven.be>
% Wim Meeussen <wim.meeussen@mech.kuleuven.be>
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function F=loop_closure_eqs(phi_init,phi1,r2,r3,r5,r6,r7,r9,r10,r11,gamma,r14x,r14y,r47y,r18x,r18y,r811y)

% first argument: the initial values of the unknown angles phi3 and phi4
% argument phi2: input angle phi2 for which we want to calculate the unknown angles phi3 and phi4
% arguments a1 ... phi1: constants

% copy initial values of unknown angles phi3 and phi4
phi2=phi_init(1);
phi3=phi_init(2);
phi5=phi_init(3);
phi6=phi_init(4);
phi8=phi_init(5);
phi10=phi_init(6);
phi11=phi_init(7);
r13=phi_init(8);
r4=phi_init(9);
r8 = phi_init(10);
phi4 = phi3 + gamma;
phi9 = phi8 + 2*pi - gamma;


% loop closure equations:
F(1)=r2*cos(phi1)+r3*cos(phi2);
F(2)=r2*sin(phi1)+r3*sin(phi2)+r13;

F(3)=-r4*cos(phi3)+r14x;
F(4)=-r14y+r13-r4*sin(phi3);

F(5)=-r8*cos(phi8)-r18x;
F(6)=-r8*sin(phi8)+r13-r18y;

F(7)=r5*cos(phi3 + gamma)+r6*cos(phi5)+r7*cos(phi6);
F(8)=r5*sin(phi3 + gamma)+r6*sin(phi5)+r7*sin(phi6)+r47y;

F(9)=r9*cos(phi8 + 2*pi -gamma)+r10*cos(phi10)+r11*cos(phi11);
F(10)=r9*sin(phi8 + 2*pi - gamma)+r10*sin(phi10)+r11*cos(phi11)+r811y;


