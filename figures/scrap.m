clear
L = 200E-9; %x direction
W = 100E-9; %y direction
Vo = 0.1;

nx = 20;
ny = 20;

V = zeros(nx,ny);
G = sparse(nx*ny,nx*ny);
F = zeros(nx*ny,1);

case1 = 1; %question 1a) boundary conditions 
%case1 = 0; %question 1b) boundary conditions

m0 = 9.1093837015E-31;
mn = 0.26*m0;
T = 300;
kB = 1.38064852E-23;
tmn = 0.2E-12;

vth = ((2*kB*T)/mn)^0.5; %thermal velocity

% mfp = vth*tmn; %mean free path
% 
% V = vth;

xmax = L; %max positions
ymax = W;


Np = 1000; % # particles, want 1000-10000
Nplot = 25;

% %injection = 1; %injection is ON
% injection = 0; %injection is OFF
% 
% if injection == 1
%     Px = zeros(Np,1); %initial positions shifted so they come from the left side of the screen
%     Vx = V*(abs(randn(Np,1))); %took away the -0.5 and changed randn back to rand so that all the inital x velocities are positive
% else
%     Px = xmax*rand(Np,1);
%     Vx = V*(randn(Np,1)-0.5);
% end

Py = ymax*rand(Np,1);
Vy = V.*(randn(Np,1)-0.5); %initial velocities

dt = 0.01*(ymax/V); %time step
a = 1; %no acceleration

tstop = 500; %simulation time

Ppx = Px; %previous postions
Ppy = Py;

c = hsv(Nplot);

box1 = 0; %turning on box 1
box2 = 0; %turning on box 2

if box1 == 1
b1 = rectangle('Position',[0.8E-7,0,0.4E-7,0.4E-7]); %box1 position (bottom)
end

if box2 == 1
b2 = rectangle('Position',[0.8E-7,0.6E-7,0.4E-7,0.4E-7]); %box2 position (top)
end

inbox = Px > 0.8E-7 & Px < 1.2E-7 & (Py < 0.4E-7 | Py > 0.6E-7);
while sum(inbox) > 0
    Px(inbox) = rand(sum(inbox),1)*xmax;
    Py(inbox) = rand(sum(inbox),1)*ymax;
    inbox = Px > 0.8E-7 & Px < 1.2E-7 & (Py < 0.4E-7 | Py > 0.6E-7);
end


for i = 1:nx
    for j = 1:ny
        n = j + (i-1)*ny; %mapping equation
        
        if i == 1
                F(n) = Vo;
                G(n,n) = 1;
        elseif i == nx
            if case1 == 1
                F(n) = 0;
                G(n,n) = 1;
            else
                F(n) = Vo;
                G(n,n) = 1;
            end
        elseif j == 1
            if case1 == 1
                F(n) = 0;
                nxm = j + ((i-1)-1)*ny; %(i-1,j)
                nxp = j + ((i+1)-1)*ny; %(i+1,j)
                nyp = (j+1) + (i-1)*ny; %(i,j+1)
                G(n,n) = -3;
                G(n,nxm) = 1;
                G(n,nxp) = 1;
                G(n,nyp) = 1;
            else
                F(n) = 0;
                G(n,n) = 1;
            end
        elseif j == ny
            if case1 == 1
                F(n) = 0;
                nxm = j + ((i-1)-1)*ny; %(i-1,j)
                nxp = j + ((i+1)-1)*ny; %(i+1,j)
                nym = (j-1) + (i-1)*ny; %(i,j-1)
                G(n,n) = -3;
                G(n,nxm) = 1;
                G(n,nxp) = 1;
                G(n,nym) = 1;
            else
                F(n) = 0;
                G(n,n) = 1;
            end
            
        else
            nxm = j + ((i-1)-1)*ny; %(i-1,j)
            nxp = j + ((i+1)-1)*ny; %(i+1,j)
            nym = (j-1) + (i-1)*ny; %(i,j-1)
            nyp = (j+1) + (i-1)*ny; %(i,j+1)
            
            G(n,n) = -4;
            G(n,nxm) = 1;
            G(n,nxp) = 1;
            G(n,nym) = 1;
            G(n,nyp) = 1;
            
        end
    end
end

axis([0 L 0 W]);

M = G\F;
for i = 1:nx
    for j = 1:ny
        n = j + (i-1)*ny;
        V(i,j) = M(n);
    end
end

figure(1)
surf(V)
title('Surface Potential Plot (n = 20)')
xlabel('Y')
ylabel('X')
zlabel('V(x,y)')


[Ey, Ex] = gradient(V); %plotting gradient
    
figure(2)
quiver(-Ex',-Ey',10) %vector arrows
title('Electric Field Gradient (n = 20)')
xlabel('Y')
ylabel('X')




for i = 1:tstop
    Ppx = Px;
    Ppy = Py;
    
    Px = Px + Vx*dt;
    Py = Py + Vy*dt;
   
    
    %xBC = 1; %x boundary conditions ON
    xBC = 0; %x boundary conditions OFF
   
    if xBC == 1
        ix1 = Px < 0;
        Px(ix1) = Px(ix1) + xmax;
        Ppx(ix1) = Ppx(ix1) + xmax;
       
        ix2 = Px > xmax;
        Px(ix2) = Px(ix2) - xmax;
        Ppx(ix2) = Ppx(ix2) - xmax;
    elseif xBC == 0
        ix1 = Px < 0;
        Px(ix1) = Px(ix1);
        Ppx(ix1) = Ppx(ix1);
        
        ix2 = Px > xmax;
        Px(ix2) = Px(ix2);
        Ppx(ix2) = Ppx(ix2);
    end
    
    iy1 = Py < 0 | Py > ymax;
    Vy(iy1) = -Vy(iy1);

    
    scatter = 1; %scattering is ON
    %scatter = 0; %scattering is OFF
    
    if scatter == 1
        Psc = 1 - exp(-(dt/tmn));
    elseif scatter == 0
        Psc = 0;
    end
    
    std = sqrt((kB*T)/mn);
    
    isc = Psc > rand(Np,1);
    Vx = Vx + a*dt;
    Vy = Vy + a*dt;
    Px = Px + Vx*dt + 0.5*a*(dt)^2;
    Py = Py + Vy*dt + 0.5*a*(dt)^2;
    Vx(isc) = randn(sum(isc),1)*std;
    Vy(isc) = randn(sum(isc),1)*std;
    
    
    if box1 == 1
        inbox1 = Px > 0.8E-7 & Px < 1.2E-7 & Py < 0.4E-7; %in box1
        LorR = Ppx < 0.8E-7 | Ppx > 1.2E-7; %coming from L or R
        Px(inbox1 & LorR) = Ppx(inbox1 & LorR);
        Vx(inbox1 & LorR) = -Vx(inbox1 & LorR);
        Py(inbox1 & ~LorR) = Ppy(inbox1 & ~LorR);
        Vy(inbox1 & ~LorR) = -Vy(inbox1 & ~LorR);
    end
    
    if box2 == 1
        inbox2 = Px > 0.8E-7 & Px < 1.2E-7 & Py > 0.6E-7; %in box2
        LorR = Ppx < 0.8E-7 | Ppx > 1.2E-7; %coming from L or R
        Px(inbox2 & LorR) = Ppx(inbox2 & LorR);
        Vx(inbox2 & LorR) = -Vx(inbox2 & LorR);
        Py(inbox2 & ~LorR) = Ppy(inbox2 & ~LorR);
        Vy(inbox2 & ~LorR) = -Vy(inbox2 & ~LorR);
    end
    
    
    figure(1)
    xlabel('X (m)')
    ylabel('Y (m)')
    hold on
    axis([0 xmax 0 ymax]);
    
    for j = 1:Nplot
    plot([Ppx(j),Px(j)]',[Ppy(j),Py(j)]','color',c(j,:));
    end
    
    
    vavg = mean(sqrt(Vx.^2 + Vy.^2)); %average velocity
    
    TSi = ((vavg.^2)*mn)/(2*kB); %temperature of the Si semiconductor
    
    pause(0.001)
end

hold off

figPlot = figure (2);
PxPy = [Px,Py];
hist3(PxPy,'CDataMode','auto','FaceColor','interp'); %electron density histogram, colour coated to show density in 2D as well
title('Electron Density Map')
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Number of Electrons')


Xbins = discretize(Px, 20); %bins for the temp map
Ybins = discretize(Py, 10); 

Tmap = zeros(20,10); %for loops checking the temp of the particles in each of the x & y bins
for x = 1:20
    for y= 1:10
    Vmean = mean(sqrt(Vx(Xbins == x & Ybins == y).^2 + Vy(Xbins == x & Ybins == y).^2));
    Tmap(x,y) = ((Vmean).^2*mn)/(2*kB);
    end 
end

figure(3)
surf(Tmap)
title('Temperature Map')
xlabel('Y (m)')
ylabel('X (m)')
zlabel('Temperature (K)')
            
            