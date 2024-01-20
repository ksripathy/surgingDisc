

function [force]=F_cal_CMI(coordinates,velocity,pressure,normaldirection)

% coordinates=[xline;yline;zline];
% velocity=UC_VEC;
% pressure=PC;
% normaldirection=NORVEC;


N=3; % dimensions of space

% normal direction
n1 = normaldirection(1,:);
n2 = normaldirection(2,:);
n3 = normaldirection(3,:);

% coordinates
x = coordinates(1,:);
y = coordinates(2,:);
z = coordinates(3,:);

% velocity terms
u = velocity(1,:);
v = velocity(2,:);
w = velocity(3,:);

% scalar pressure 
p=pressure; 


f1 = [ -n1.*(u.^2)-n2.*u.*v-n3.*u.*w  ;...
       -n1.*u.*v-n2.*(v.^2)-n3.*v.*w ;...
         -n1.*u.*w-n2.*v.*w-n3.*(w.^2)]; 
     
f2 = [ -n1.*p  ;...
       -n2.*p ;...
       -n3.*p];   
   
   % % % total force 
   force=f1+f2; 
     