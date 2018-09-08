function [J]= pid_optim(x)

s=tf('s');
plant=1/(s^3+6*(s^2)+5*s);

Kp=x(1)
Kd=x(2)
Ki=x(3)

cont = Kp+Ki/s+Kd*s;
%step(feedback(plant*cont,1))
dt=0.01;
t=0:dt:1;
e=1-step(feedback(plant*cont,1),t);
%ITAE
J=sum(t'.*abs(e)*dt);
