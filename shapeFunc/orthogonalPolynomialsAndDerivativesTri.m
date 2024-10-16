function [P,dPdxi,dPdeta]=orthogonalPolynomialsAndDerivativesTri(degree,Xi)

u = Xi(:,1); v=Xi(:,2);

U = ones(size(u)); O = zeros(size(u));

switch degree
    case 1
        P=[U, 1/2+(3/2)*v, 1/2+u+(1/2)*v];
        dPdxi=[O, O, U];
        dPdeta=[O, U*3/2, U/2];
    case 2
        P=[ U, 1/2+(3/2)*v, -1/2+v+(5/2)*v.^2, 1/2+u+(1/2)*v, (3/2+(5/2)*v).*u+3/4+(5/4)*v.^2+2*v, (3/2).*u.^2+((3/2)*v+3/2).*u+(1/4)*v.^2+v+1/4];
        dPdxi=[ O, O, O, U, 3/2+(5/2)*v, 3*u+(3/2)*v+3/2];
        dPdeta=[ O, U*3/2, 1+5*v, U/2, (5/2).*u+(5/2)*v+2, (3/2).*u+(1/2)*v+1];
    case 3
        P=[ U, 1/2+(3/2)*v, -1/2+v+(5/2)*v.^2, -3/8-(15/8)*v+(15/8)*v.^2+(35/8)*v.^3, 1/2+u+(1/2)*v, (3/2+(5/2)*v).*u+3/4+(5/4)*v.^2+2*v, (1/4+(9/2)*v+(21/4)*v.^2).*u+1/8+(21/8)*v.^3+(39/8)*v.^2+(19/8)*v, (3/2).*u.^2+((3/2)*v+3/2).*u+(1/4)*v.^2+v+1/4, (15/4+(21/4)*v).*u.^2+(15/4+(21/4)*v.^2+9*v).*u+5/8+(7/8)*v.^3+(33/8)*v.^2+(27/8)*v, (5/2).*u.^3+((15/4)*v+15/4).*u.^2+((3/2)*v.^2+(9/2)*v+3/2).*u+(1/8)*v.^3+(9/8)*v.^2+(9/8)*v+1/8];
        dPdxi=[O, O, O, O, U, 3/2+(5/2)*v, 1/4+(9/2)*v+(21/4)*v.^2, 3*u+(3/2)*v+3/2, (15/2+(21/2)*v).*u+15/4+(21/4)*v.^2+9*v, (15/2).*u.^2+((15/2)*v+15/2).*u+(3/2)*v.^2+(9/2)*v+3/2 ];
        dPdeta=[ O, U*3/2, 1+5*v, -15/8+(15/4)*v+(105/8)*v.^2, U/2, (5/2)*u+(5/2)*v+2, (9/2+(21/2)*v).*u+(63/8)*v.^2+(39/4)*v+19/8, (3/2)*u+(1/2)*v+1, (21/4)*u.^2+((21/2)*v+9).*u+(21/8)*v.^2+(33/4)*v+27/8, (15/4)*u.^2+(3*v+9/2).*u+(3/8)*v.^2+(9/4)*v+9/8];
    case 4
        P=[ U, 1/2+(3/2)*v, -1/2+v+(5/2)*v.^2, -3/8-(15/8)*v+(15/8)*v.^2+(35/8)*v.^3, 3/8-(3/2)*v-(21/4)*v.^2+(7/2)*v.^3+(63/8)*v.^4, 1/2+u+(1/2)*v, (3/2+(5/2)*v).*u+3/4+(5/4)*v.^2+2*v, (1/4+(9/2)*v+(21/4)*v.^2).*u+1/8+(21/8)*v.^3+(39/8)*v.^2+(19/8)*v, (-1+(21/2)*v.^2+(21/2)*v.^3).*u-1/2+(21/4)*v.^4+(21/2)*v.^3+(21/4)*v.^2-(1/2)*v, (3/2).*u.^2+((3/2)*v+3/2).*u+(1/4)*v.^2+v+1/4, (15/4+(21/4)*v).*u.^2+(15/4+(21/4)*v.^2+9*v).*u+5/8+(7/8)*v.^3+(33/8)*v.^2+(27/8)*v, ((27/2)*v.^2+15*v+3).*u.^2+(3+(27/2)*v.^3+(57/2)*v.^2+18*v).*u+1/2+(9/4)*v.^4+(23/2)*v.^3+(51/4)*v.^2+(9/2)*v, (5/2).*u.^3+((15/4)*v+15/4).*u.^2+((3/2)*v.^2+(9/2)*v+3/2).*u+(1/8)*v.^3+(9/8)*v.^2+(9/8)*v+1/8, (35/4+(45/4)*v).*u.^3+(105/8+(135/8)*v.^2+30*v).*u.^2+(21/4+(27/4)*v.^3+(51/2)*v.^2+(45/2)*v).*u+7/16+(9/16)*v.^4+(11/2)*v.^3+9*v.^2+(9/2)*v, (35/8).*u.^4+((35/4)*v+35/4).*u.^3+((45/8)*v.^2+15*v+45/8).*u.^2+((5/4)*v.^3+(15/2)*v.^2+(15/2)*v+5/4).*u+(1/16)*v.^4+v.^3+(9/4)*v.^2+v+1/16];
        dPdxi=[ O, O, O, O, O, U, 3/2+(5/2)*v, 1/4+(9/2)*v+(21/4)*v.^2, -1+(21/2)*v.^2+(21/2)*v.^3, 3*u+(3/2)*v+3/2, (15/2+(21/2)*v).*u+15/4+(21/4)*v.^2+9*v, (27*v.^2+30*v+6).*u+3+(27/2)*v.^3+(57/2)*v.^2+18*v, (15/2).*u.^2+((15/2)*v+15/2).*u+(3/2)*v.^2+(9/2)*v+3/2, (105/4+(135/4)*v).*u.^2+(105/4+(135/4)*v.^2+60*v).*u+21/4+(27/4)*v.^3+(51/2)*v.^2+(45/2)*v, (35/2).*u.^3+((105/4)*v+105/4).*u.^2+((45/4)*v.^2+30*v+45/4).*u+(5/4)*v.^3+(15/2)*v.^2+(15/2)*v+5/4];
        dPdeta=[ O, U*3/2, 1+5*v, -15/8+(15/4)*v+(105/8)*v.^2, -3/2-(21/2)*v+(21/2)*v.^2+(63/2)*v.^3, U/2, (5/2).*u+(5/2)*v+2, ...
            (9/2+(21/2)*v).*u+(63/8)*v.^2+(39/4)*v+19/8, (21*v+(63/2)*v.^2).*u+21*v.^3+(63/2)*v.^2+(21/2)*v-1/2, (3/2).*u+(1/2)*v+1, ...
            (21/4).*u.^2+((21/2)*v+9).*u+(21/8)*v.^2+(33/4)*v+27/8, (27*v+15).*u.^2+((81/2)*v.^2+57*v+18).*u+9*v.^3+(69/2)*v.^2+(51/2)*v+9/2, ...
            (15/4).*u.^2+(3*v+9/2).*u+(3/8)*v.^2+(9/4)*v+9/8, ...
            (45/4).*u.^3+((135/4)*v+30).*u.^2+((81/4)*v.^2+51*v+45/2).*u+(9/4)*v.^3+(33/2)*v.^2+18*v+9/2, ...
            (35/4).*u.^3+((45/4)*v+15).*u.^2+((15/4)*v.^2+15*v+15/2).*u+(1/4)*v.^3+3*v.^2+(9/2)*v+1];
    case 5
        P=[ U, 1/2+(3/2)*v, -1/2+v+(5/2)*v.^2, -3/8-(15/8)*v+(15/8)*v.^2+(35/8)*v.^3, 3/8-(3/2)*v-(21/4)*v.^2+(7/2)*v.^3+(63/8)*v.^4, 5/16+(35/16)*v-(35/8)*v.^2-(105/8)*v.^3+(105/16)*v.^4+(231/16)*v.^5, 1/2+u+(1/2)*v, (3/2+(5/2)*v).*u+3/4+(5/4)*v.^2+2*v, (1/4+(9/2)*v+(21/4)*v.^2).*u+1/8+(21/8)*v.^3+(39/8)*v.^2+(19/8)*v, (-1+(21/2)*v.^2+(21/2)*v.^3).*u-1/2+(21/4)*v.^4+(21/2)*v.^3+(21/4)*v.^2-(1/2)*v, (-3/8-(11/2)*v-(9/4)*v.^2+(45/2)*v.^3+(165/8)*v.^4).*u-3/16+(165/16)*v.^5+(345/16)*v.^4+(81/8)*v.^3-(31/8)*v.^2-(47/16)*v, (3/2).*u.^2+((3/2)*v+3/2).*u+(1/4)*v.^2+v+1/4, (15/4+(21/4)*v).*u.^2+(15/4+(21/4)*v.^2+9*v).*u+5/8+(7/8)*v.^3+(33/8)*v.^2+(27/8)*v, ((27/2)*v.^2+15*v+3).*u.^2+(3+(27/2)*v.^3+(57/2)*v.^2+18*v).*u+1/2+(9/4)*v.^4+(23/2)*v.^3+(51/4)*v.^2+(9/2)*v, ((495/16)*v.^3+(675/16)*v.^2+(189/16)*v-15/16).*u.^2+(-15/16+(495/16)*v.^4+(585/8)*v.^3+54*v.^2+(87/8)*v).*u-5/32+(165/32)*v.^5+(885/32)*v.^4+(141/4)*v.^3+(59/4)*v.^2+(43/32)*v, (5/2).*u.^3+((15/4)*v+15/4).*u.^2+((3/2)*v.^2+(9/2)*v+3/2).*u+(1/8)*v.^3+(9/8)*v.^2+(9/8)*v+1/8, (35/4+(45/4)*v).*u.^3+(105/8+(135/8)*v.^2+30*v).*u.^2+(21/4+(27/4)*v.^3+(51/2)*v.^2+(45/2)*v).*u+7/16+(9/16)*v.^4+(11/2)*v.^3+9*v.^2+(9/2)*v, ((275/8)*v.^2+(175/4)*v+95/8).*u.^3+(285/16+(825/16)*v.^3+(1875/16)*v.^2+(1335/16)*v).*u.^2+(57/8+(165/8)*v.^4+(705/8)*v.^3+(213/2)*v.^2+(381/8)*v).*u+19/32+(55/32)*v.^5+(565/32)*v.^4+(143/4)*v.^3+(107/4)*v.^2+(241/32)*v, (35/8).*u.^4+((35/4)*v+35/4).*u.^3+((45/8)*v.^2+15*v+45/8).*u.^2+((5/4)*v.^3+(15/2)*v.^2+(15/2)*v+5/4).*u+(1/16)*v.^4+v.^3+(9/4)*v.^2+v+1/16, (315/16+(385/16)*v).*u.^4+(315/8+(385/8)*v.^2+(175/2)*v).*u.^3+(405/16+(495/16)*v.^3+(1725/16)*v.^2+(1575/16)*v).*u.^2+(45/8+(55/8)*v.^4+(375/8)*v.^3+75*v.^2+(325/8)*v).*u+9/32+(11/32)*v.^5+(185/32)*v.^4+(135/8)*v.^3+(125/8)*v.^2+(155/32)*v, (63/8).*u.^5+((315/16)*v+315/16).*u.^4+((35/2)*v.^2+(175/4)*v+35/2).*u.^3+((105/16)*v.^3+(525/16)*v.^2+(525/16)*v+105/16).*u.^2+((15/16)*v.^4+(75/8)*v.^3+(75/4)*v.^2+(75/8)*v+15/16).*u+(1/32)*v.^5+(25/32)*v.^4+(25/8)*v.^3+(25/8)*v.^2+(25/32)*v+1/32];
        dPdxi=[O, O, O, O, O, O, U, 3/2+(5/2)*v, 1/4+(9/2)*v+(21/4)*v.^2, -1+(21/2)*v.^2+(21/2)*v.^3, -3/8-(11/2)*v-(9/4)*v.^2+(45/2)*v.^3+(165/8)*v.^4, 3*u+(3/2)*v+3/2, (15/2+(21/2)*v).*u+15/4+(21/4)*v.^2+9*v, (27*v.^2+30*v+6).*u+3+(27/2)*v.^3+(57/2)*v.^2+18*v, ((495/8)*v.^3+(675/8)*v.^2+(189/8)*v-15/8).*u-15/16+(495/16)*v.^4+(585/8)*v.^3+54*v.^2+(87/8)*v, (15/2).*u.^2+((15/2)*v+15/2).*u+(3/2)*v.^2+(9/2)*v+3/2, (105/4+(135/4)*v).*u.^2+(105/4+(135/4)*v.^2+60*v).*u+21/4+(27/4)*v.^3+(51/2)*v.^2+(45/2)*v, ((825/8)*v.^2+(525/4)*v+285/8).*u.^2+(285/8+(825/8)*v.^3+(1875/8)*v.^2+(1335/8)*v).*u+57/8+(165/8)*v.^4+(705/8)*v.^3+(213/2)*v.^2+(381/8)*v, (35/2).*u.^3+((105/4)*v+105/4).*u.^2+((45/4)*v.^2+30*v+45/4).*u+(5/4)*v.^3+(15/2)*v.^2+(15/2)*v+5/4, (315/4+(385/4)*v).*u.^3+(945/8+(1155/8)*v.^2+(525/2)*v).*u.^2+(405/8+(495/8)*v.^3+(1725/8)*v.^2+(1575/8)*v).*u+45/8+(55/8)*v.^4+(375/8)*v.^3+75*v.^2+(325/8)*v, (315/8).*u.^4+((315/4)*v+315/4).*u.^3+((105/2)*v.^2+(525/4)*v+105/2).*u.^2+((105/8)*v.^3+(525/8)*v.^2+(525/8)*v+105/8).*u+(15/16)*v.^4+(75/8)*v.^3+(75/4)*v.^2+(75/8)*v+15/16 ];
        dPdeta=[O, U*3/2, 1+5*v, -15/8+(15/4)*v+(105/8)*v.^2, -3/2-(21/2)*v+(21/2)*v.^2+(63/2)*v.^3, 35/16-(35/4)*v-(315/8)*v.^2+(105/4)*v.^3+(1155/16)*v.^4, U/2, (5/2).*u+(5/2)*v+2, (9/2+(21/2)*v).*u+(63/8)*v.^2+(39/4)*v+19/8, (21*v+(63/2)*v.^2).*u+21*v.^3+(63/2)*v.^2+(21/2)*v-1/2, (-11/2-(9/2)*v+(135/2)*v.^2+(165/2)*v.^3).*u+(825/16)*v.^4+(345/4)*v.^3+(243/8)*v.^2-(31/4)*v-47/16, (3/2).*u+(1/2)*v+1, (21/4).*u.^2+((21/2)*v+9).*u+(21/8)*v.^2+(33/4)*v+27/8, (27*v+15).*u.^2+((81/2)*v.^2+57*v+18).*u+9*v.^3+(69/2)*v.^2+(51/2)*v+9/2, ((1485/16)*v.^2+(675/8)*v+189/16).*u.^2+((495/4)*v.^3+(1755/8)*v.^2+108*v+87/8).*u+(825/32)*v.^4+(885/8)*v.^3+(423/4)*v.^2+(59/2)*v+43/32, (15/4).*u.^2+(3*v+9/2).*u+(3/8)*v.^2+(9/4)*v+9/8, (45/4).*u.^3+((135/4)*v+30).*u.^2+((81/4)*v.^2+51*v+45/2).*u+(9/4)*v.^3+(33/2)*v.^2+18*v+9/2, ((275/4)*v+175/4).*u.^3+((2475/16)*v.^2+(1875/8)*v+1335/16).*u.^2+((165/2)*v.^3+(2115/8)*v.^2+213*v+381/8).*u+(275/32)*v.^4+(565/8)*v.^3+(429/4)*v.^2+(107/2)*v+241/32, (35/4).*u.^3+((45/4)*v+15).*u.^2+((15/4)*v.^2+15*v+15/2).*u+(1/4)*v.^3+3*v.^2+(9/2)*v+1, (385/16).*u.^4+((385/4)*v+175/2).*u.^3+((1485/16)*v.^2+(1725/8)*v+1575/16).*u.^2+((55/2)*v.^3+(1125/8)*v.^2+150*v+325/8).*u+(55/32)*v.^4+(185/8)*v.^3+(405/8)*v.^2+(125/4)*v+155/32, (315/16).*u.^4+(35*v+175/4).*u.^3+((315/16)*v.^2+(525/8)*v+525/16).*u.^2+((15/4)*v.^3+(225/8)*v.^2+(75/2)*v+75/8).*u+(5/32)*v.^4+(25/8)*v.^3+(75/8)*v.^2+(25/4)*v+25/32 ];
    otherwise
        error('Degree not implemented in orthogonalPolynomialsAndDerivativesTri.m')
end
