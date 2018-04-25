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


function F=loop_closure_eqs(phi_init,phi2,r2,r3,r5,r6,r7,r9,r10,r11,gamma,r14x,r14y,r47y,r18x,r18y,r811y)

% first argument: the initial values of the unknown angles phi3 and phi4
% argument phi2: input angle phi2 for which we want to calculate the unknown angles phi3 and phi4
% arguments a1 ... phi1: constants

% copy initial values of unknown angles phi3 and phi4
phi3=phi_init(1);
phi4=phi_init(2);
phi6=phi_init(3);
phi7=phi_init(4);
phi8=phi_init(5);
phi10=phi_init(6);
phi11=phi_init(7);
r13=phi_init(8);
r4=phi_init(9);
r8 = phi_init(10);
phi5 = phi4 + gamma;
phi9 = phi8 + 2*pi - gamma;


% loop closure equations:
F(1)=r2*cos(phi2)-r3*cos(phi3);
F(2)=r2*sin(phi2)-r3*sin(phi3)+r13;

F(3)=-r4*cos(phi4)+r14x;
F(4)=-r14y+r13-r4*sin(phi4);

F(5)=-r8*cos(phi8)-r18x;
F(6)=-r8*sin(phi8)+r13-r18y;

F(7)=r5*cos(phi4 + gamma)+r6*cos(phi6)+r7*cos(phi7);
F(8)=r5*sin(phi4 + gamma)+r6*sin(phi6)+r7*sin(phi7)+r47y;

F(9)=+r9*cos(phi8 + 2*pi -gamma)+r10*cos(phi10)+r11*cos(phi11);
F(10)=+r9*sin(phi8 + 2*pi - gamma)+r10*sin(phi10)+r11*sin(phi11)+r811y;


