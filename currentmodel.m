%toy model to simulate and demonstrate electric current in a metal
%assume 2D to keep it simple to visualise
%notation may be a bit harder to read as all particles are kept in a long
%vector (eg. x is the vector of x coordinates for all charge carriers in the simulation)
%this allows an important speed gain in the simulation (more optimisation
%is definitely possible)

%exercises, play with number of particles and study current behavior
%increase field Ex and see what happens with current
%does ohms law hold? (also for small nr of particles?, also for very high
%field strengths?)

%work out mobility and compare to prediction in course notes
%and compare to the statistical result for J

%what does it say about measuring really small currents in small devices?
%or really high field strengths?

close all
clear

%create a metal rectangle
w=5e-6; %width in [m]
h=w;% height in [m]

%create a set of charge carriers
q=-1.6e-19; %charge of electrons [C]
m=9.1e-31; %electron rest mass [kg]
np=50; %number of particles in the box
vth=1e5; %thermal speed of electrons (assume constant here, in reality it is a stochastic variable [m/s], relate to temp
l=100e-9; %free path length [m]
tau=l/vth; %free time [s]

%define electric field over the box in x direction
Ex=1000000; %electric field in x direction [V/m]
Ey=0; %electric field in y direction [V/m]

%create random initial positions of particles
xinit=rand([np,1])*w;
yinit=rand([np,1])*h;
thetainit=rand([np,1])*2*pi; %random angle in xy plane

%plot the particles in the box
plot(xinit,yinit,'r.')
axis image

%run over time (take enough time points to keep accuracy, do convergence test, this depends also on E field)
tpoints=100;
tmax=10*tau;
t=ones([np,1])*linspace(0,tmax,tpoints); %each particle has its own timescale that gets reset when a collision occurs
x=xinit;
y=yinit;
theta=thetainit;
vx=vth*cos(theta);
vy=vth*sin(theta);
dt=tmax/tpoints;

for tid=1:tpoints,
    %calculate new point some time later   
    x=x+vx*dt;
    y=y+vy*dt;
    %and acceleration term due to E field
    vx=vx+dt*q*Ex/(m);
    vy=vy+dt*q*Ey/(m);
    
    %check if they hit the wall and keep them inside box by modulo
    x=mod(x,w); %particle that dissapear through top wall will reappear bottom to keep nr of particles constant
    y=mod(y,h); %particle that dissapear through right wall will reappear left to keep nr of particles constant
     
    %reset speed and direction in case a collision with impurity occurs
    collid=find(rand([np,1])>(exp(-t(:,tid)/tau))); %a vector with indices for those particles that have a collision now (statistical model for exponential chance of colliding)
    %collid=find(t(:,tid)>tau*ones([np,1])); %same but classical collision all at tau
     
    t(collid,:)=t(collid,:)-t(collid,tid)*ones([1,tpoints]); %reset time to zero for those particles that had a collision
    theta(collid)=rand(size(theta(collid)))*2*pi; %new angle for those particles that had a collision
    vx(collid)=vth.*cos(theta(collid)); %in case of collision take new speed, otherwise keep old speed
    vy(collid)=vth.*sin(theta(collid));
    
    %calculate current
    Jx(tid)=np*q*mean(vx); %instantenous current (note the random fluctuation which should get less for more particles)
    Jy(tid)=np*q*mean(vy); %instantenous current (note the random fluctuation which should get less for more particles)
    
    %show
    plot(x,y,'r.',x(collid),y(collid),'b.','MarkerSize',4); %plot particle positions, plot the ones that collided as blue cross (or otherwise in corner)
    hold on %trick to keep previous plot and plot over it
    axis image
    pause(0.01); %this forces a refresh of the plot during the loop
end


J=mean(Jx)
sigma=mean(Jx)/Ex %conductivity estimate (with stochastic noise)
mu=sigma/(np*q) %mobility estimate

mus=q*tau/m %result for statistical model in course notes
sigmas=np*q^2*tau/m %result for statistical model in course notes

figure
plot(t(1,:),Jx,'b',t(1,:),sigmas*Ex*ones(size(Jx)),'r')
title('instantenous current')
xlabel('t');
ylabel('J')
