

function [force]=F_cal_CMI_unsteady(coordinates,velocity,pressure,normaldirection,velocity_tp1,velocity_tm1,DT)

% coordinates=[xline;yline;zline];
% velocity=UC_VEC;
% pressure=PC;
% normaldirection=NORVEC;
% velocity_tp1=UC_VEC_tp1;
% velocity_tm1=UC_VEC_tm1;
% DT=DT;


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


% % % velocity at time t+1
u_tp1 = velocity_tp1(1,:);
v_tp1 = velocity_tp1(2,:);
w_tp1 = velocity_tp1(3,:);

% % % velocity at time t+1
u_tm1 = velocity_tm1(1,:);
v_tm1 = velocity_tm1(2,:);
w_tm1 = velocity_tm1(3,:);

% scalar pressure 
p=pressure; 


f1 = [ -n1.*(u.^2)-n2.*u.*v-n3.*u.*w  ;...
       -n1.*u.*v-n2.*(v.^2)-n3.*v.*w ;...
         -n1.*u.*w-n2.*v.*w-n3.*(w.^2)]; 
     
f2 = [ -n1.*p  ;...
       -n2.*p ;...
       -n3.*p];   
   
   
% % % % force contribution from the unsteady term 
f3_tp1= [ -n1.*(x.*u_tp1)-n2.*x.*v_tp1-n3.*x.*w_tp1  ;...
       -n1.*u_tp1.*y-n2.*(v_tp1.*y)-n3.*y.*w_tp1 ;...
         -n1.*u_tp1.*z-n2.*v_tp1.*z-n3.*(w_tp1.*z)]; 
 
f3_tm1= [ -n1.*(x.*u_tm1)-n2.*x.*v_tm1-n3.*x.*w_tm1  ;...
       -n1.*u_tm1.*y-n2.*(v_tm1.*y)-n3.*y.*w_tm1 ;...
         -n1.*u_tm1.*z-n2.*v_tm1.*z-n3.*(w_tm1.*z)];      
     
f3= (f3_tp1-f3_tm1)/DT;    
   
   
   
   % % % total force 
   force=f1+f2+f3; 
     