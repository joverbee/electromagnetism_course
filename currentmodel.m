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
%switch models between 'classic' and 'statistical' to see differences,
%judge which one would be closer to reality especially when studying the
%instanteneous current

%what does it say about measuring really small currents in small devices (low np)?
%or really high field strengths?

%Jo Verbeeck, EMAT, University of Antwerp, Nov 2017
close all
clear

rng('shuffle'); %make sure random gen does not repeat itself

%model='classic';
model='statistical';

%create a square 2D conductor
w=5e-6; %width in [m]
h=w;% height in [m]

%create a set of charge carriers
q=-1.6e-19; %charge of electrons [C]
m=9.1e-31; %electron rest mass [kg]
np=50; %number of particles in the box
kB=1.38e-23;%Boltzman constant [m^2kgs^-2K-1]
T=300; %absolute temperature [K]
%vth=1e5; %thermal speed of electrons (assume constant here, in reality it is a stochastic variable [m/s], todo: relate to temp 
vth=sqrt(3*kB*T/m); %3D electron gas assumed, to be strict there are only 2 degrees of freedom which would turn the '3' into a '2'
l=100e-9; %free path length [m] (depends on purity of materials and density of defects, remember: no material is ever really perfect)
tau=l/vth; %free time [s]

%define electric field over the box in x direction
Ex=100000; %electric field in x direction [V/m] (note above 1e5 very few materials would survive, but parabolic trajectory becomes hard to notice)
Ey=0; %electric field in y direction [V/m]

%create random initial positions of particles
xinit=rand([np,1])*w;
yinit=rand([np,1])*h;
thetainit=rand([np,1])*2*pi; %random angle in xy plane

%plot the particles in the box
plot(xinit,yinit,'r.')
axis image

%run over time (take enough time points to keep accuracy, do convergence test, this depends also on E field)
tpoints=1000;
tmax=10*tau;
t=ones([np,1])*linspace(0,tmax,tpoints); %each particle has its own timescale that gets reset when a collision occurs
x=xinit;
y=yinit;
theta=thetainit;

switch (model)
    case 'statistical',
        vx=vth*cos(theta)+tau*q*Ex/m; %assume the particle has a past, it already gained an average acceleration over a time tau (approximate but better than assuming 0)
        vy=vth*sin(theta)+tau*q*Ey/m;
    case 'classic',
        vx=vth*cos(theta)+tau*q*Ex/(2*m); %assume the particle has a past, it already gained an average acceleration over a time tau (approximate but better than assuming 0)
        vy=vth*sin(theta)+tau*q*Ey/(2*m);
        
end

dt=tmax/tpoints;

for tid=1:tpoints,
    %calculate new point some time later   
    x=x+vx*dt;
    y=y+vy*dt;
    %and acceleration term due to E field
    vx=vx+ones(size(vx))*dt*q*Ex/m;
    vy=vy+ones(size(vy))*dt*q*Ey/m;
    
    %calculate current
    Jx(tid)=np*q*mean(vx); %instantenous current (note the random fluctuation which should get less for more particles)
    Jy(tid)=np*q*mean(vy); %instantenous current (note the random fluctuation which should get less for more particles)
    
    %check if they hit the wall and keep them inside box by modulo
    x=mod(x,w); %particle that dissapear through top wall will reappear bottom to keep nr of particles constant
    y=mod(y,h); %particle that dissapear through right wall will reappear left to keep nr of particles constant
     
    %reset speed and direction in case a collision with impurity occurs
    switch (model)
        case 'statistical',
            %careful with this formula, contains statistical difficulties,
            %as we do a statistical draw in each time point which needs to
            %be compensated (otherwise more time points would lead to
            %faster collision)        
            collid=find(rand([np,1])>(exp(-t(:,tid)./(tau*(t(:,tid)/dt+ones([np,1])))))); %a vector with indices for those particles that have a collision now (statistical model for exponential chance of colliding)
        case 'classic',
            collid=find(t(:,tid)>tau*ones([np,1])); %classical collision all when t>tau
        otherwise,   
            %unknown model
            break
    end
    mean(t(collid,tid)) %should be close to tau if statistics is correct
    t(collid,:)=t(collid,:)-t(collid,tid)*ones([1,tpoints]); %reset time to zero for those particles that had a collision
    theta(collid)=rand(size(theta(collid)))*2*pi; %new angle for those particles that had a collision
    vx(collid)=vth.*cos(theta(collid)); %in case of collision take new speed, otherwise keep old speed
    vy(collid)=vth.*sin(theta(collid));
    
    
    
    %show
    plot(x,y,'r.',x(collid),y(collid),'b.','MarkerSize',4); %plot particle positions, plot the ones that collided as blue cross (or otherwise in corner)
    hold on %trick to keep previous plot and plot over it
    axis image
    pause(0.01); %this forces a refresh of the plot during the loop
end


J=mean(Jx)
sigma=mean(Jx)/Ex %conductivity estimate (with stochastic noise)
mu=sigma/(np*q) %mobility estimate

%analytical results from course notes
mus=q*tau/m 
sigmas=np*q^2*tau/m 
muc=q*tau/(2*m) 
sigmac=np*q^2*tau/(2*m) 

figure
plot(t(1,:),Jx,'b',t(1,:),sigmas*Ex*ones(size(Jx)),'r',t(1,:),sigmac*Ex*ones(size(Jx)),'g')
title('instantenous current')
xlabel('t');
ylabel('J');
legend('Jx','statistical prediction','classical prediction')
