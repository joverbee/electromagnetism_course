function wos()
%walk on spheres test program
%aim: understand the usefulness and performance and details of walk on
%spheres for solving Laplace equation with boundary conditions

%we walk over a grid of points in which we want to know V, but note that
%this grid is only for display reasons and the method works just as well
%for a single point (as opposed to other numerical solutions to Laplace)

%essence is based on the mean value theorem stating that the potential in
%the center of sphere is equal to the average potential over the sphere
%wall as long as there are no charges inside the sphere.
%we choose a point in space
%then we choose the largest sphere from there, just touching a conductor
%we take a arbitrary step onto the sphere surface
%if this step brings us close to boundary than we take the potential of
%that boundary as an approximation to V in the center of the sphere (one
%vote for e.g. V1)
%if the step takes us to a point that is not on the boundary, we create a
%new sphere around that new point and take a next step from there until
%after a few more spheres we hit a boundary and this gives a vote for the
%potential of that boundary and the votex automatically takes into account
%how likely it was that we did this sequence of steps
%the more walks we make, the more accurate we sample the statistics of the
%problem and the closer our result will be to the 'real' result.
%Jo Verbeeck, EMAT, University of Antwerp 2017.

%TO DO: can we expand the method to give directly the field E instead of
%the potential? For trajectory calculation we need the force in a specific
%point and doing a numerical gradient seems to contradict the philosophy of
%this method.

close all
clear

%for demo case: only 2D


%define a circular disc electrode
V1=10; %potential
r1=1; %radius of disc
p1=[1,1]; %position of center

%define another disc electrode
V2=2;
r2=0.5;
p2=[-2,-1];


%show electrodes to give impression of situation (note the grid is ONLY for
%visualisation
xmax=3;
xpoints=50;
x=linspace(-xmax,xmax,xpoints);
y=x;
[X,Y]=meshgrid(x,y);
elec1=V1*indisc(X,Y,p1,r1);
elec2=V2*indisc(X,Y,p2,r2);

figure
imagesc(x,y,elec1+elec2)
title('position of electrodes')
hold on

kpoints=100; %number of walks per point (higher is better and more important than increasing the recursion level)
eps=xmax/xpoints;

Vgrid=elec1+elec2; %fill in potential that is already known
vac=(~indisc(X,Y,p1,r1)).*(~indisc(X,Y,p2,r2));
idvac=find(vac); 
gridpoints=size(idvac)-1;

for l=1:gridpoints,
        %choose a point in space
        pn=[X(idvac(l)),Y(idvac(l))];       
        Vn=0; %original estimate of V
        %find the largest sphere just touching an electrode from here
        rn=min(dist(pn,p1,r1),dist(pn,p2,r2));
        %and draw this circle
        %viscircles(pn,rn);
        
        %make some steps from this point and collect potential on surface of this
        %circle        
        for k=1:kpoints
            i=0;
            Vn=makestep(Vn,pn,rn,p1,r1,V1,p2,r2,V2,i,eps);
        end
        Vgrid(idvac(l))=Vn/kpoints; %store potential
end

figure
imagesc(x,y,Vgrid);

[Ex,Ey]=gradient(-Vgrid,x,y);
laplace=divergence(X,Y,Ex,Ey);
figure
imagesc(x,y,laplace);
title('laplacian(V) should be 0')
colorbar


%and find its potential with walk on spheres random walks
end
%-------------------------
function hit = indisc(x,y,p,r)
%returns true if coordinate x,y is inside the circular electrode position
%at xpos,ypos and with radius r1
hit=sqrt( (p(:,1)-x).^2+ (p(:,2)-y).^2 )<r;
end

function r= dist(rn,p,r)
%returns distance to electrode from x,y to circular electrode position
%at xpos,ypos and with radius r1
r=sqrt( (rn(:,1)-p(:,1)).^2+  (rn(:,2)-p(:,2)).^2 )-r;
end

function Vn=makestep(Vn,pn,rn,p1,r1,V1,p2,r2,V2,i,eps)
maxrecurse=7;
theta=rand(1)*2*pi; %arbitrary direction of walk
    pn2=pn+rn.*[cos(theta),sin(theta)]; %new position somewhere on circle
    if (dist(pn2,p1,r1)<eps)
        %we reached a point very close to boundary 1, assign potential of
        %that boundary
        Vn=Vn+V1;
    elseif (dist(pn2,p2,r2)<eps)
        %we reached a point very close to boundary 2, assign potential of
        %that boundary
        Vn=Vn+V2;
    else
        rn2=min(dist(pn2,p1,r1),dist(pn2,p2,r2));
        %viscircles(pn2,rn2); %show the new circle 
        %if no boundary found, take new circle from here (recursion)
        i=i+1; %keep track of recursion depth
        if (i<maxrecurse) 
            Vn=makestep(Vn,pn2,rn2,p1,r1,V1,p2,r2,V2,i,eps);
        end
    end
end



