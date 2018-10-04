%demonstration of electrostatic shielding in a multiconductor case
%make use of numerical solution of Laplace equation with boundary
%conditions, to show redistribution of charges for a collection of
%conductors and how inside a shielding box (electrode 2) the fields and charges inside
%are unaffected by the outside charges. Note that the surface charges on
%the inside of the shielding box DO change in order to keep the internal
%field and the charges on electrode 1 from changing when external changes (V3 and V4) are applied
%Jo Verbeeck, EMAT, University of Antwerp, oct 2017.
clear all
close all

%define shape of 4 electrodes, work in 2D to keep it simple (3D takes much
%longer to calculate)

%define grid
xmax=1; %extension of grid [m]
xpoints=256; %number of grid points to calculate
x=linspace(-xmax,xmax,xpoints);
y=x;
[X,Y]=meshgrid(x,y);
R=sqrt(X.^2+Y.^2);


%define 4 electrodes with funny shapes
gap=xmax/8; %defines gap
w=xmax/10;
%logic functions being 1 where a certain electrode is defined and zero
%elsewhere
electrode1=(R<gap/2).*(X<0); %half disc inside shielding ring
electrode2=(R>gap).*(R<(gap+w)); %shielding ring
electrode3=sqrt((X-4*gap).^2+Y.^2)<gap; %disc outside shielding
electrode4=sqrt((Y-4*gap).^2+X.^2)<gap; %another disc outside shielding
%and the potential on each electrode
V1=10;
V2=0; %assune shielding at gnd potential (arbitrary, but this is convention)
V3=1;
V4=10;


vac=(~electrode1).*(~electrode2).*(~electrode3).*(~electrode4); %logic function defining the vacuum region where no electrodes are present
V=V1*electrode1+V2*electrode2+V3*electrode3+V4*electrode4; %potential boundary conditions=known constant potential on the electrodes
idelectrode=find(~vac); %index for all points that contain electrodes
idvac=find(vac); %index for all points in vacuum

figure
imagesc(x,y,~vac)
axis image
title('position of electrodes')

figure
imagesc(x,y,V)
title('potential on electrodes')
colorbar
axis image

%solve the Laplace equation for vacuum with boundary conditions
%first rough attempt: all vacuum points to average potential on electrodes
%(this is clearly too crude but we need to start somewhere)
V(idvac)=mean(V(idelectrode)); %replave vacuum potential by the mean potential on the electrodes

%second step: relaxation method
Vn=V; %starting condition is first rough estimate
kmax=100000; %number of itteration steps

K=[1/4,0,1/4;0,0,0;1/4,0,1/4]; %Laplace kernel replacing the central point by the average of 4 neighbouring points. 
% This is a crude numerical approximation to the mean value theorem stating
% that the potential inside a sphere is equal to the average potential on
% the surface of that sphere if no charges are present inside the sphere

%faster algorithms exist that start with a rough grid and refine the grid
%in the process, also the first guess for potential could be a lot smarter
%when e.g. taking a weighted average depending on distance

for k=1:kmax,
    Vold=Vn;
    Vn=conv2(Vold,K,'same'); %apply kernel to old estimate
    Vn(idelectrode)=V(idelectrode);%and reinforce the boundary conditions
    dif=Vold-Vn;%difference between old and new estimate
    error(k)=sum(dif(:).^2); %square difference as indicator of convergence
   
    %imagesc(Vn); %see how the potential evolves while calculating (but
    %takes a lot of time that should really go to the calculation and not
    %to plotting
    %pause(0.1); %to force refreshing everytime we come here, otherwise the
    %image is not updated
end

figure
semilogy(error)
title('mean squared correction')

figure
imagesc(x,y,Vn)
colorbar
title('estimated V')
axis image

%get E field from this
[Ex,Ey]=gradient(-Vn,x,y);
E=sqrt(Ex.^2+Ey.^2);

%figure
%quiver(X,Y,Ex,Ey);

figure
imagesc(x,y,E)
colorbar();
title('|E|')
axis image

figure
eps0=8.854e-12;
rho=eps0*divergence(X,Y,Ex,Ey);
imagesc(x,y,rho)
colorbar();
title('charge density \rho')
axis image


