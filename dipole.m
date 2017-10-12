%visualisation of dipole field made with 2 charged conducting spheres
%show E field as vectors, as E field strenght |E|, as streamlines
%show potential and equipotential surface
%Jo Verbeeck, EMAT, University of Antwerp,  oct 2017
clear all
close all

%create grid to plot (choose 2D plane for visualisation cutting through charge centers , but calculation is correct for 3D)
xpoints=256; %nr of grid points in 1 direction
xmax=1; %extension of grid [m]
pref=9e9; %1/(4pi eps0)
x=linspace(-xmax,xmax,xpoints);
y=x;
[X,Y]=meshgrid(x,y); %2D matrices holding X or Y coordinate for each point on the grid
R=sqrt(X.^2+Y.^2);

%define 2 conducting spheres with opposite charge 
q=1e-10; %charge on spheres [C]
radius=0.2*xmax; %radius of spheres [m]

q1=q; %charge on sphere 1 [C]
x1=-0.5*xmax; %x position [m]
y1=0; %y position [m]

R1=sqrt((X-x1).^2+(Y-y1).^2); %distance to charge 1
in1=R1<radius; %logical function 1 if inside sphere, 0 if outside
out1=R1>=radius; %logical function 0 if inside sphere, 1 if outside

E1x=pref*q1*((X-x1)./(R1).^3); %x component of E field valid outside sphere
E1y=pref*q1*((Y-y1)./(R1).^3); %y component of E field valid outside sphere
V1=pref*q1*(ones(size(R1))./R1); % potential 1/r valid outside sphere


q2=-q; %charge on sphere 2 [C]
x2=+0.5*xmax; %x position [m]
y2=0; %y position [m]
R2=sqrt((X-x2).^2+(Y-y2).^2); %distance to charge
out2=R2>=radius; %logical function 1 if inside sphere, 0 if outside
in2=R2<radius; %logical function 0 if inside sphere, 1 if outside
E2x=pref*q2*(X-x2)./(R2).^3; %x component of E field valid outside sphere
E2y=pref*q2*(Y-y2)./(R2).^3; %x component of E field valid outside sphere
V2=pref*q2*ones(size(R2))./R2; % potential 1/r valid outside sphere

%and apply superposition of both fields (but take care of special cases
%inside spheres)
Ex=(E1x+E2x).*(out1.*out2); %E field only if outside both spheres (logical and), inside sphere E=0 because conductor [V/m]
Ey=(E1y+E2y).*(out1.*out2); 
E=sqrt(Ex.^2+Ey.^2); %and size of E
V=(V1+V2).*(out1.*out2)+in1*pref*q1/radius+in2*pref*q2/radius; %potential, inside spheres is constant (conductors)

%and compare to approximate formula of dipole field (only valid far away
%from charges), note how it really only starts to compare to real dipole
%when far away from the dipole
px=q*(x1-x2); %x component of dipole moment [Cm]
py=q*(y1-y2); %x component of dipole moment [Cm]
Vdip=pref*((px*X+py*Y)./R.^3);
Exdip=-pref*(px./R.^3-3*(px*X+py*Y).*X./R.^5);
Eydip=-pref*(py./R.^3-3*(px*X+py*Y).*Y./R.^5);
Edip=sqrt(Exdip.^2+Eydip.^2);

%and showtime!
figure
%show vector plot, but limit number of points to keep the number of vector
%reasonable
scale=1;
nsteps=25; %total number of vector points in each direction
step=round(xpoints/nsteps);
quiver(X(1:step:xpoints,1:step:xpoints),Y(1:step:xpoints,1:step:xpoints),Ex(1:step:xpoints,1:step:xpoints),Ey(1:step:xpoints,1:step:xpoints),scale);
title('electric field')
xlabel('x');
ylabel('y');
axis image

figure
imagesc(x,y,E)
title('electric field and fieldlines')
xlabel('x');
ylabel('y');
axis image
hold on
hLines = streamslice(X,Y,Ex,Ey);
set(hLines,'Color','r');
colorbar

figure
imagesc(x,y,Edip,[min(E(:)), max(E(:))])
title('electric field and fieldlines approximate dipole')
xlabel('x');
ylabel('y');
axis image
hold on
hLines = streamslice(X,Y,Exdip,Eydip);
set(hLines,'Color','r');
colorbar


figure
imagesc(x,y,V)
title('V')
xlabel('x');
ylabel('y');
axis image
colorbar


figure
nlines=50;
contour(x,y,V,nlines)
title('equipotential surfaces')
xlabel('x');
ylabel('y');
axis image
colorbar