%% FDTD file edited and modified by Raul Barroso (Simon Bolivar University 2025) from the one created by Nathan Zechar Wright of State University 2021, 
%called "Intro to FDTD", which is published in the MATLAB webpage

%% In this file, the electric field inside the water volume to be heated is calculated as a function of time in a coordinate for both cases, 
%% with and without the designed metamaterial matching system. Once this algorithm is finished, it is necessary to clear the graphics that appear there and run the file "Efficiency1216.m" and from a previous estimation of the gain between the two mentioned systems, the efficiency of the matching system is finally estimated, assuming that the power efficiency of the system without metamaterial slab is 36%. The final obtained graphic is the power efficiency as a function of frequency in hertz.

clear all; close all;
%% Define Simulation Based off Source and Wavelength
f0 = 2.45e9; % Frequency of Source [Hertz]
Lf = round(800); % Divisions per Wavelength in vaccum at the Frequency of Source [unitless]
[Lx,Ly,Lz] = deal(0.028*19/20,5,0.028*19/20); % Wavelengths x,y,z [unitless] %0.084
nt = round(15000); % Number of time steps [unitless]
%nt =5;
arista = round(12);% This is the value of "b" given in FDTD cells number. For the paper "a" was fixed in 16 FDTD cells. As "b=12", then b/a = 0.75.
GapY=0;
t = 0;
% nt =1800 %
%% Simulation Performance
single = 0; % Use Single Precision | 1 = yes , 0 = no |
usegpu = 0; % Use gpuArray | 1 = yes , 0 = no |
%% Spatial and Temporal System
e0 = 8.854*10^-12; % Permittivity of vacuum [farad/meter]
u0 = 4*pi*10^-7; % Permeability of vacuum [henry/meter]
c0 = 1/(e0*u0)^.5; % Speed of light [meter/second]
sgCu= 5.8*10^7; % Copper Conductivity;
sgH2O = 0.7;% Water Conductivity;
eCu = e0; %permittivity of copper
eH2O = 80*e0; %permittivity of water
eta0= (u0/e0)^.5; 
etaH2O= real((u0/(eH2O*(1-1i*(sgH2O/(eH2O*2*3.1415923*f0)))))^.5)+1i*imag((u0/(eH2O*(1-1i*(sgH2O/(eH2O*2*3.1415923*f0)))))^.5); %Intrinsic impedance of water
etaH2OReal = real((u0/(eH2O*(1-1i*(sgH2O/(eH2O*2*3.1415923*f0)))))^.5); 
etaH2OImag = imag((u0/(eH2O*(1-1i*(sgH2O/(eH2O*2*3.1415923*f0)))))^.5);

etaCu= (u0/(eCu*(1-1i*(sgCu/(eCu*2*3.1415923*f0)))))^.5;
etametateorico = (etaH2O*eta0)^.5;
L0 = c0/f0; % Freespace Wavelength [meter]
t0 = 1/f0; % Source Period [second]
[Nx,Ny,Nz] = deal(round(Lx*Lf),round(Ly*Lf),round(Lz*Lf)); % Points in x,y,z [unitless]
x = linspace(0,Lx,Nx+1)*L0; % x vector [meter]
y = linspace(0,Ly,Ny+1)*L0; % y vector [meter]
z = linspace(0,Lz,Nz+1)*L0; % z vector [meter]
[dx,dy,dz] = deal(x(2),y(2),z(2)); % x,y,z increment [meter]
dt = (dx^-2+dy^-2+dz^-2)^-.5/c0*.99; % Time step CFL condition [second]

%% Initialize Incident Magnetic and Electric Field Vectors and Materials

[Exi] = deal(zeros(Nx,Ny+1,Nz+1)); % Ex cells
[Eyi] = deal(zeros(Nx+1,Ny,Nz+1)); % Ey cells
[Ezi] = deal(zeros(Nx+1,Ny+1,Nz)); % Ez cells
[Hxi] = deal(zeros(Nx+1,Ny,Nz)); % Hx cells
[Hyi] = deal(zeros(Nx,Ny+1,Nz)); % Hy cells
[Hzi] = deal(zeros(Nx,Ny,Nz+1)); % Hz cells
[udx,udy,udz] = deal(dt/(u0*dx),dt/(u0*dy),dt/(u0*dz)); % Hi vaccum coeffcients
[edx,edy,edz] = deal(dt/(e0*dx),dt/(e0*dy),dt/(e0*dz)); % Ei vaccum coeffcients

[udxs,udys,udzs] = deal(dt/(u0*dx),dt/(u0*dy),dt/(u0*dz)); % Hs vaccum coeffcients
[edxs,edys,edzs] = deal(dt/(e0*dx),dt/(e0*dy),dt/(e0*dz)); % Es vaccum coeffcients
[edxCu,edyCu,edzCu] = deal(2*dt/((2*e0+sgCu*dt)*dx),2*dt/((2*e0+sgCu*dt)*dy),2*dt/((2*e0+sgCu*dt)*dz)); % Es curl material Coeffcients
[edxH2O,edyH2O,edzH2O] = deal(dt/((eH2O+sgH2O*dt)*dx), dt/((eH2O+sgH2O*dt)*dy), dt/((eH2O+sgH2O*dt)*dz)); % Es curl material Coeffcients


eta0= (u0/e0)^.5;
etaH2O= real(u0/(eH2O*(1-1i*(sgH2O/(eH2O*2*3.1415923*f0)))))^.5;
etaCu= real(u0/(e0*(1-1i*(sgCu/(e0*2*3.1415923*f0)))))^.5;

CexeCu = eCu/(eCu+sgCu*dt);
CeicCu = (e0-eCu-sgCu*dt)/(eCu+sgCu*dt);
CeicCu2 = (etaCu-eta0)/(etaCu+eta0);
CeipCu = -(e0-eCu)/(eCu+sgCu*dt);

CexeH2O = eH2O/(eH2O+sgH2O*dt);
CeicH2O = (e0-eH2O-sgH2O*dt)/(eH2O+sgH2O*dt);
CeipH2O = -(e0-eH2O)/(eH2O+sgH2O*dt);


[Material] = deal(zeros(Nx+1,Ny+1,Nz+1)); % Materials, conductors or vaccum

Material(:,Ny+1,:) = 2; % Este es el final del sistema en Y, con PEC
%Material(:,1,:) = 3; % Este es el inicio del sistema en Y, con PMC
%Material(:,Ny,:) = 3; % Este es el final del sistema en Y, con PMC
%Material(:,Ny-1,:) = 3; % Este es el final del sistema en Y, con PMC
Material(1,:,:) = 3; %Pared izquierda
Material(2,:,:) = 3; %Pared izquierda
Material(Nx,:,:) = 3; %Pared derecha
Material(Nx+1,:,:) = 3; %Pared derecha
Material(:,:,1) = 2; %Pared inferior
Material(:,:,2) = 2; %Pared inferior
Material(:,:,Nz-1) = 2; %Pared superior
Material(:,:,Nz) = 2; %Pared superior
Material(:,:,Nz+1) = 2; %Pared superior

Nobj =10;
Dsepx = arista;
Dsepz = arista;
Dsepy = arista;
amenosb=(Nz-Dsepz-4.5);
a = Nx-5;
b = Dsepx;

etametadiseno = 377*((a-b)*(a+b)*(a-b)/(a*a*a))^.5;
impevista = etametadiseno*etametadiseno/etaH2O;
relab=b/a;
vdisenoNorm = (((a-b)*(a+b)/(a*a))*(a/(a-b)))^(-0.5);
Ndivmetadiseno = Lf*(1/4)*vdisenoNorm;
modcoefrefldiseno = abs((impevista-etametadiseno)/(impevista+etametadiseno));
roediseno = (1+modcoefrefldiseno)/(1-modcoefrefldiseno);
CeldaY=round(Dsepy+amenosb-1);
Ezi0(nt)=0;
Ezt(nt)=0;
MaxEz(Ny)=0;
MinEz(Ny)=0;

%% Routine to build the joint of conductor cubes. Each conductor cube is performed with a joint of FDTD cells.

for M=1:Nobj;
for I=1:Nx+1;% Objeto conductor rectangular
    for J=1:Ny+1;
        for K=1:Nz+1;
                %( K > (Nz/2-(Dsepz/2))) && ( K < (Nz/2+(Dsepz/2)) )
            if ( ( I > (1+Nx/2-(Dsepx/2))) && ( I < (2+Nx/2+(Dsepx/2))) && ( J > ((2*Ny/3-Dsepy-GapY)-(M-1)*CeldaY)) && ( J < ((2*Ny/3)-GapY-(M-1)*CeldaY)) && ( K > (Nz/2-(Dsepz/2))) && ( K < (Nz/2+(Dsepz/2)) ));
            Material(I,J,K)=2; 
            end
           
        end
    end
end
end


for I=1:Nx+1;% Objeto conductor rectangular
    for J=1:Ny+1;
        for K=1:Nz+1;
           
            if ( ( I > (2+Nx/2-Dsepx/2)) && ( I < (1+Nx/2+Dsepx/2)) && ( J > (Ny/2)) && ( J < (Ny-3)) && ( K > (Nz/2-Dsepz/2)) && ( K < (Nz/2+Dsepz/2)) );
            %Material(I,J,K)=4;
            end
           
        end
    end
end

for I=1:Nx+1;% Medio volumen lleno de material con pérdidas
    for J=1:Ny+1;
        for K=1:Nz+1;
           
            if ( ( J > (Ny/2)) && (Material(I,J,K)==0) );
%            Material(I,J,K)=4;
            end
           
        end
    end
end


%% Initialize Scattered Magnetic and Electric Field Vectors

[Exs] = deal(zeros(Nx,Ny+1,Nz+1)); % Ex cells
[Eys] = deal(zeros(Nx+1,Ny,Nz+1)); % Ey cells
[Ezs] = deal(zeros(Nx+1,Ny+1,Nz)); % Ez cells
[Hxs] = deal(zeros(Nx+1,Ny,Nz)); % Hx cells
[Hys] = deal(zeros(Nx,Ny+1,Nz)); % Hy cells
[Hzs] = deal(zeros(Nx,Ny,Nz+1)); % Hz cells

%% Performance Enchancement
%if single == 1; vars2Single(); end
%if usegpu == 1; gpu2Single(); end
%% Start loop
while  t <= nt % Starting a loop of time
    tic
t = t + 1;

%% Incident Magnetic Field Update
Hxi = Hxi +0*udz.*diff(Eyi,1,3)- udy.*diff(Ezi,1,2);
%Hyi = Hyi+0*udz.*diff(Ezi,1,1)-0*udx.*diff(Exi,1,3);
Hzi = Hzi+0*udy.*diff(Exi,1,2)-0*udx.*diff(Eyi,1,1);

%for I=1:Nx+1;%PMC para Hx en y=0 y en y=Ly
%for J=1:Ny-1;    
%for K=1:Nz;  
       
%         if (Material(I,J,K)==3 && Material(I,J+1,K)==0) %  Reversa                                          
%         Hxi(I,J+1,K) = -Hxi(I,J,K);
%         else
%         end
         
%         if (Material(I,J,K)==0 && Material(I,J+1,K)==3) %  Directa                                          
%         Hxi(I,J,K) = -Hxi(I,J+1,K);
%         else
%         end

%end
%end
%end



%for I=1:Nx;%PMC para Hz en y=0 y en y=Ly
%for J=1:Ny-1;    
%for K=1:Nz+1;  
       
%         if (Material(I,J,K)==3 && Material(I,J+1,K)==0) %  Reversa                                          
%         Hzi(I,J+1,K) = -Hzi(I,J,K);
%         else
%         end
         
%         if (Material(I,J,K)==0 && Material(I,J+1,K)==3) %  Directa                                          
%         Hzi(I,J,K) = -Hzi(I,J+1,K);
%         else
%         end
%end
%end
%end


%% Incident Electric Field Update


Exianterior = Exi;
Eyianterior = Eyi;
Ezianterior = Ezi;



    


for I=1:Nx;
for J=2:Ny;    
for K=2:Nz;
   if ( J < (2*Ny/3)) 
    Exi(I,J,K) = Exi(I,J,K)+ edy*(Hzi(I,J,K)-Hzi(I,J-1,K))- edz*(Hyi(I,J,K)-Hyi(I,J,K-1));
   else
    Exi(I,J,K) = CexeH2O*Exi(I,J,K) + edyH2O*(Hzi(I,J,K)-Hzi(I,J-1,K))- edzH2O*(Hyi(I,J,K)-Hyi(I,J,K-1));       
   end
end
end
end

for I=2:Nx;
for J=1:Ny;    
for K=2:Nz;  
   if ( J < (2*Ny/3))
   Eyi(I,J,K) = Eyi(I,J,K) + edz*(Hxi(I,J,K)-Hxi(I,J,K-1))- edx*(Hzi(I,J,K)-Hzi(I-1,J,K));
   else
   Eyi(I,J,K) = CexeH2O*Eyi(I,J,K) + edzH2O*(Hxi(I,J,K)-Hxi(I,J,K-1))- edxH2O*(Hzi(I,J,K)-Hzi(I-1,J,K));   
   end
end
end
end    

for I=2:Nx;
for J=2:Ny;    
for K=1:Nz;
   if ( J < (2*Ny/3))    
   Ezi(I,J,K) = Ezi(I,J,K) + edx*(Hyi(I,J,K)-Hyi(I-1,J,K))- edy*(Hxi(I,J,K)-Hxi(I,J-1,K));
   else
   Ezi(I,J,K) = CexeH2O*Ezi(I,J,K) + edxH2O*(Hyi(I,J,K)-Hyi(I-1,J,K))- edyH2O*(Hxi(I,J,K)-Hxi(I,J-1,K));       
   end
end
end
end





%[Exi] = deal(zeros(Nx,Ny+1,Nz+1)); % Ex cells
%[Eyi] = deal(zeros(Nx+1,Ny,Nz+1)); % Ey cells
%[Ezi] = deal(zeros(Nx+1,Ny+1,Nz)); % Ez cells

%for I=1:Nx;%Perfect conductor at the end
%for J=1:Ny+1;    
%for K=1:Nz+1;  
       
%         if (Material(I,J,K)==2) %                                            
%         Exi(I,J,K) = -Exi(I,J,K);
%         else
%         end

%end
%end
%end

for I=1:Nx+1;%Perfect conductor at the end
for J=1:Ny+1;    
for K=1:Nz;  
       
         if (J==(Ny+1)); %                                            
         Ezi(I,J,K) = -Ezi(I,J,K);
         else
         end

         if (J==(Ny)); %                                            
%         Ezi(I,J,K) = -Ezi(I,J,K);
         else
         end        
end
end
end

%[Exi] = deal(zeros(Nx,Ny+1,Nz+1)); % Ex cells
%[Eyi] = deal(zeros(Nx+1,Ny,Nz+1)); % Ey cells
%[Ezi] = deal(zeros(Nx+1,Ny+1,Nz)); % Ez cells

%for I=1:Nx+1;%Perfect conductor at other border
%for J=1:Ny;    
%for K=1:Nz+1;  
       
%         if (Material(I,J,K)==2) %                                            
%         Eyi(I,J,K) = -Eyi(I,J,K);
%         else
%         end

%end
%end
%end

%% Scattered Magnetic Field Update
Hxs = Hxs+udz.*diff(Eys,1,3)-udy.*diff(Ezs,1,2);
Hys = Hys+udz.*diff(Ezs,1,1)-udx.*diff(Exs,1,3);
Hzs = Hzs+udy.*diff(Exs,1,2)-udx.*diff(Eys,1,1);

%[Hxs] = deal(zeros(Nx+1,Ny,Nz)); % Hx cells
%[Hys] = deal(zeros(Nx,Ny+1,Nz)); % Hy cells
%[Hzs] = deal(zeros(Nx,Ny,Nz+1)); % Hz cells

for I=1:Nx+1;%PMC para Hx en y=0 y en y=Ly
for J=1:Ny-1;    
for K=1:Nz;  
       
         if (Material(I,J,K)==3 && Material(I,J+1,K)==0) %  Reversa                                          
         Hxs(I,J,K) = -Hxi(I,J,K);
         Hxs(I,J+1,K) = -Hxi(I,J+1,K);        
         else
         end


         if (Material(I,J,K)==3 && Material(I,J+1,K)==3) %  Directa                                          
         Hxs(I,J+1,K) = -Hxi(I,J+1,K);
         Hxs(I,J,K) = -Hxi(I,J,K);
         else
         end        
         
         if (Material(I,J,K)==0 && Material(I,J+1,K)==3) %  Directa                                          
         Hxs(I,J+1,K) = -Hxi(I,J+1,K);
         Hxs(I,J,K) = -Hxi(I,J,K);
         else
         end

end
end
end

%[Hxs] = deal(zeros(Nx+1,Ny,Nz)); % Hx cells
%[Hys] = deal(zeros(Nx,Ny+1,Nz)); % Hy cells
%[Hzs] = deal(zeros(Nx,Ny,Nz+1)); % Hz cells

for I=1:Nx+1;%PMC para Hx en z=0 y en z=Lz
for J=1:Ny;    
for K=1:Nz-1;  
       
         if (Material(I,J,K)==3 && Material(I,J,K+1)==0) %  Reversa                                          
         Hxs(I,J,K) = -Hxi(I,J,K);
         Hxs(I,J,K+1) = -Hxi(I,J,K+1);        
         else
         end


         if (Material(I,J,K)==3 && Material(I,J,K+1)==3) %  Directa                                          
         Hxs(I,J,K+1) = -Hxi(I,J,K+1);
         Hxs(I,J,K) = -Hxi(I,J,K);
         else
         end        
         
         if (Material(I,J,K)==0 && Material(I,J,K+1)==3) %  Directa                                          
         Hxs(I,J,K+1) = -Hxi(I,J,K+1);
         Hxs(I,J,K) = -Hxi(I,J,K);
         else
         end

end
end
end




%[Hxs] = deal(zeros(Nx+1,Ny,Nz)); % Hx cells
%[Hys] = deal(zeros(Nx,Ny+1,Nz)); % Hy cells
%[Hzs] = deal(zeros(Nx,Ny,Nz+1)); % Hz cells

for I=1:Nx-1;%PMC para Hy en x=0 y en x=Ly
for J=1:Ny+1;    
for K=1:Nz;  
       
         if (Material(I,J,K)==3 && Material(I+1,J,K)==0) %  Reversa                                          
         Hys(I,J,K) = -Hyi(I,J,K);
         Hys(I+1,J,K) = -Hyi(I+1,J,K);
         else
         end

         if (Material(I,J,K)==3 && Material(I+1,J,K)==3) %  Directa                                          
         Hys(I,J,K) = -Hyi(I,J,K);        
         Hys(I+1,J,K) = -Hyi(I+1,J,K);
         else
         end        
         
         if (Material(I,J,K)==0 && Material(I+1,J,K)==3) %  Directa                                          
         Hys(I,J,K) = -Hyi(I,J,K);
         Hys(I+1,J,K) = -Hyi(I+1,J,K);
         else
         end

end
end
end

%[Hxs] = deal(zeros(Nx+1,Ny,Nz)); % Hx cells
%[Hys] = deal(zeros(Nx,Ny+1,Nz)); % Hy cells
%[Hzs] = deal(zeros(Nx,Ny,Nz+1)); % Hz cells

for I=1:Nx;%PMC para Hy en z=0 y en z=Lz
for J=1:Ny+1;    
for K=1:Nz-1;  
       
         if (Material(I,J,K)==3 && Material(I,J,K+1)==0) %  Reversa                                          
         Hys(I,J,K) = -Hyi(I,J,K);
         Hys(I,J,K+1) = -Hyi(I,J,K+1);
         else
         end

         if (Material(I,J,K)==3 && Material(I,J,K+1)==3) %  Directa                                          
         Hys(I,J,K) = -Hyi(I,J,K);
         Hys(I,J,K+1) = -Hyi(I,J,K+1);
         else
         end        
         
         if (Material(I,J,K)==0 && Material(I,J,K+1)==3) %  Directa                                          
         Hys(I,J,K) = -Hyi(I,J,K);
         Hys(I,J,K+1) = -Hyi(I,J,K+1);
         else
         end

end
end
end

%[Hxs] = deal(zeros(Nx+1,Ny,Nz)); % Hx cells
%[Hys] = deal(zeros(Nx,Ny+1,Nz)); % Hy cells
%[Hzs] = deal(zeros(Nx,Ny,Nz+1)); % Hz cells

for I=1:Nx;%PMC para Hz en y=0 y en y=Ly
for J=1:Ny-1;    
for K=1:Nz+1;  
       
         if (Material(I,J,K)==3 && Material(I,J+1,K)==0) %  Reversa                                          
         Hzs(I,J,K) = -Hzi(I,J,K);
         Hzs(I,J+1,K) = -Hzi(I,J+1,K);
         else
         end

         if (Material(I,J,K)==3 && Material(I,J+1,K)==3) %  Directa                                          
         Hzs(I,J,K) = -Hzi(I,J,K);
         Hzs(I,J+1,K) = -Hzi(I,J+1,K);
         else
         end        
         
         if (Material(I,J,K)==0 && Material(I,J+1,K)==3) %  Directa                                          
         Hzs(I,J,K) = -Hzi(I,J,K);
         Hzs(I,J+1,K) = -Hzi(I,J+1,K);
         else
         end
end
end
end


%[Hxs] = deal(zeros(Nx+1,Ny,Nz)); % Hx cells
%[Hys] = deal(zeros(Nx,Ny+1,Nz)); % Hy cells
%[Hzs] = deal(zeros(Nx,Ny,Nz+1)); % Hz cells


for I=1:Nx-1;%PMC para Hz en x=0 y en x=Lx
for J=1:Ny;    
for K=1:Nz+1;  
       
         if (Material(I,J,K)==3 && Material(I+1,J,K)==0) %  Reversa                                          
         Hzs(I,J,K) = -Hzi(I,J,K);
         Hzs(I+1,J,K) = -Hzi(I+1,J,K);        
         else
         end

         if (Material(I,J,K)==3 && Material(I+1,J,K)==3) %  Reversa                                          
         Hzs(I,J,K) = -Hzi(I,J,K);
         Hzs(I+1,J,K) = -Hzi(I+1,J,K);
         else
         end        
         
         if (Material(I,J,K)==0 && Material(I+1,J,K)==3) %  Directa                                          
         Hzs(I,J,K) = -Hzi(I,J,K);
         Hzs(I+1,J,K) = -Hzi(I+1,J,K);
         else
         end
end
end
end


%% Scattered Electric Field Update


 


for I=1:Nx;
for J=2:Ny;    
for K=2:Nz;
   if ( J < (2*Ny/3)) 
Exs(I,J,K) = Exs(I,J,K) + edys*(Hzs(I,J,K)-Hzs(I,J-1,K))- edzs*(Hys(I,J,K)-Hys(I,J,K-1));
   else
Exs(I,J,K) = CexeH2O*Exs(I,J,K) + edyH2O*(Hzs(I,J,K)-Hzs(I,J-1,K))- edzH2O*(Hys(I,J,K)-Hys(I,J,K-1));       
   end
end
end
end

for I=2:Nx;
for J=1:Ny;    
for K=2:Nz;  
   if ( J < (2*Ny/3))
Eys(I,J,K) = Eys(I,J,K) + edzs*(Hxs(I,J,K)-Hxs(I,J,K-1))- edxs*(Hzs(I,J,K)-Hzs(I-1,J,K));
   else
Eys(I,J,K) = CexeH2O*Eys(I,J,K) + edzH2O*(Hxs(I,J,K)-Hxs(I,J,K-1))- edxH2O*(Hzs(I,J,K)-Hzs(I-1,J,K));       
   end
end
end
end    

for I=2:Nx;
for J=2:Ny;    
for K=1:Nz;
   if ( J < (2*Ny/3))    
Ezs(I,J,K) = Ezs(I,J,K) + edxs*(Hys(I,J,K)-Hys(I-1,J,K))- edys*(Hxs(I,J,K)-Hxs(I,J-1,K));
   else
Ezs(I,J,K) = CexeH2O*Ezs(I,J,K) + edxH2O*(Hys(I,J,K)-Hys(I-1,J,K))- edyH2O*(Hxs(I,J,K)-Hxs(I,J-1,K));       
   end
end
end
end






%Exs(I,J,K) = Exs(I,J,K) + edys*(Hzs(I,J,K)-Hzs(I,J-1,K))- edzs*(Hys(I,J,K)-Hys(I,J,K-1));
%Exs(I,J,K) = CexeH2O*Exs(I,J,K)+ CeicH2O*Exi(I,J,K) + CeipH2O*Exianterior(I,J,K) + edyH2O*(Hzs(I,J,K)-Hzs(I,J-1,K))- edzH2O*(Hys(I,J,K)-Hys(I,J,K-1));

%Eys(I,J,K) = Eys(I,J,K) + edzs*(Hxs(I,J,K)-Hxs(I,J,K-1))- edxs*(Hzs(I,J,K)-Hzs(I-1,J,K));
%Eys(I,J,K) = CexeH2O*Eys(I,J,K) +CeicH2O*Eyi(I,J,K) + CeipH2O*Eyianterior(I,J,K) + edzH2O*(Hxs(I,J,K)-Hxs(I,J,K-1))- edxH2O*(Hzs(I,J,K)-Hzs(I-1,J,K));

%Ezs(I,J,K) = Ezs(I,J,K) + edxs*(Hys(I,J,K)-Hys(I-1,J,K))- edys*(Hxs(I,J,K)-Hxs(I,J-1,K));
%Ezs(I,J,K) = CexeH2O*Ezs(I,J,K) + CeicH2O*Ezi(I,J,K)+ CeipH2O*Ezianterior(I,J,K) + edxH2O*(Hys(I,J,K)-Hys(I-1,J,K))- edyH2O*(Hxs(I,J,K)-Hxs(I,J-1,K));


%% Scattered Electric Field Update in Water

for I=1:Nx;%Ex calculated in the volume and in the six faces
for J=1:Ny;    
for K=1:Nz+1;                

         if (Material(I,J,K)==0 && Material(I,J+1,K)==4) %  Directa                                          
         Exs(I,J,K) = CexeH2O*Exs(I,J,K) + CeicH2O*Exi(I,J,K) + CeipH2O*Exianterior(I,J,K) + edyH2O*(Hzs(I,J,K)-Hzs(I,J-1,K))- edzH2O*(Hys(I,J,K)-Hys(I,J,K-1));
%        Exs(I,J,K) = -Exi(I,J,K);
%        Exs(I,J+1,K) = CexeH2O*Exs(I,J+1,K) + CeicH2O*Exi(I,J+1,K) + CeipH2O*Exianterior(I,J,K) + edyH2O*(Hzs(I,J+1,K)-Hzs(I,J,K))- edzH2O*(Hys(I,J+1,K)-Hys(I,J+1,K-1));
%        Exs(I,J+1,K) = -Exi(I,J+1,K);
         else
         end
         if (Material(I,J,K)==4 && Material(I,J+1,K)==4) % vol                                          
 %       Exs(I,J,K) = -Exi(I,J,K);
         Exs(I,J,K) = CexeH2O*Exs(I,J,K) +CeicH2O*Exi(I,J,K) + CeipH2O*Exianterior(I,J,K) + edyH2O*(Hzs(I,J,K)-Hzs(I,J-1,K))- edzH2O*(Hys(I,J,K)-Hys(I,J,K-1));
         else
         end
         if ((Material(I,J,K)==4 && Material(I,J+1,K)==0)) %  Reversa                                          
         Exs(I,J,K) = CexeH2O*Exs(I,J,K)+ CeicH2O*Exi(I,J,K) + CeipH2O*Exianterior(I,J,K) + edyH2O*(Hzs(I,J,K)-Hzs(I,J-1,K))- edzH2O*(Hys(I,J,K)-Hys(I,J,K-1));
%        Exs(I,J,K) = -Exi(I,J,K);
         Exs(I,J+1,K) = CexeH2O*Exs(I,J+1,K) + CeicH2O*Exi(I,J+1,K) + CeipH2O*Exianterior(I,J,K) + edyH2O*(Hzs(I,J+1,K)-Hzs(I,J,K))- edzH2O*(Hys(I,J+1,K)-Hys(I,J+1,K-1));
%        Exs(I,J+1,K) = -Exi(I,J+1,K);
         else
         end
               
       
end
end
end
for I=1:Nx;
for J=1:Ny+1;    
for K=1:Nz;  
         if (Material(I,J,K)==4 && Material(I,J,K+1)==0) %  Reversa                                          
         Exs(I,J,K) = CexeH2O*Exs(I,J,K) + CeicH2O*Exi(I,J,K) + CeipH2O*Exianterior(I,J,K) + edyH2O*(Hzs(I,J,K)-Hzs(I,J-1,K))- edzH2O*(Hys(I,J,K)-Hys(I,J,K-1));
         Exs(I,J,K+1) = CexeH2O*Exs(I,J,K+1) + CeicH2O*Exi(I,J,K+1) + CeipH2O*Exianterior(I,J,K) + edyH2O*(Hzs(I,J,K+1)-Hzs(I,J-1,K+1))- edzH2O*(Hys(I,J,K+1)-Hys(I,J,K));
         else
         end
         if (Material(I,J,K)==4 && Material(I,J,K+1)==4) %  vol                                          
         Exs(I,J,K) = CexeH2O*Exs(I,J,K) +CeicH2O*Exi(I,J,K) + CeipH2O*Exianterior(I,J,K) + edyH2O*(Hzs(I,J,K)-Hzs(I,J-1,K))- edzH2O*(Hys(I,J,K)-Hys(I,J,K-1));
         else
         end        
         if (Material(I,J,K)==0 && Material(I,J,K+1)==4) %  Directa                                          
         Exs(I,J,K) = CexeH2O*Exs(I,J,K) + CeicH2O*Exi(I,J,K) + CeipH2O*Exianterior(I,J,K) + edyH2O*(Hzs(I,J,K)-Hzs(I,J-1,K))- edzH2O*(Hys(I,J,K)-Hys(I,J,K-1));
         %Exs(I,J,K+1) = CexeH2O*Exs(I,J,K+1) + CeicH2O*Exi(I,J,K+1) + CeipH2O*Exianterior(I,J,K) + edyH2O*(Hzs(I,J,K+1)-Hzs(I,J-1,K+1))- edzH2O*(Hys(I,J,K+1)-Hys(I,J,K));
         else
         end        
end
end
end




for I=1:Nx;%Ey calculated in the volume and in the six faces
for J=1:Ny;    
for K=1:Nz+1;                        
         if (Material(I,J,K)==4 && Material(I+1,J,K)==0) %  Reversa                                          
%         Eys(I,J,K) = -Eyi(I,J,K);
         Eys(I,J,K) = CexeH2O*Eys(I,J,K) +CeicH2O*Eyi(I,J,K) + CeipH2O*Eyianterior(I,J,K) + edzH2O*(Hxs(I,J,K)-Hxs(I,J,K-1))- edxH2O*(Hzs(I,J,K)-Hzs(I-1,J,K));
%         Eys(I+1,J,K) = -Eyi(I+1,J,K);
         Eys(I+1,J,K) = CexeH2O*Eys(I+1,J,K) +CeicH2O*Eyi(I+1,J,K) + CeipH2O*Eyianterior(I,J,K) + edzH2O*(Hxs(I+1,J,K)-Hxs(I+1,J,K-1))- edxH2O*(Hzs(I+1,J,K)-Hzs(I,J,K));
         else
         end
         if (Material(I,J,K)==4 && Material(I+1,J,K)==4) %  vol                                          
         Eys(I,J,K) = CexeH2O*Eys(I,J,K) +CeicH2O*Eyi(I,J,K) + CeipH2O*Eyianterior(I,J,K) + edzH2O*(Hxs(I,J,K)-Hxs(I,J,K-1))- edxH2O*(Hzs(I,J,K)-Hzs(I-1,J,K));
         else
         end                  
         if (Material(I,J,K)==0 && Material(I+1,J,K)==4) %  Directa                                          
%         Eys(I,J,K) = -Eyi(I,J,K);
         Eys(I,J,K) = CexeH2O*Eys(I,J,K) +CeicH2O*Eyi(I,J,K) + CeipH2O*Eyianterior(I,J,K) + edzH2O*(Hxs(I,J,K)-Hxs(I,J,K-1))- edxH2O*(Hzs(I,J,K)-Hzs(I-1,J,K));
%         Eys(I+1,J,K) = -Eyi(I+1,J,K);
         %Eys(I+1,J,K) = CexeH2O*Eys(I+1,J,K) +CeicH2O*Eyi(I+1,J,K) + CeipH2O*Eyianterior(I,J,K) + edzH2O*(Hxs(I+1,J,K)-Hxs(I+1,J,K-1))- edxH2O*(Hzs(I+1,J,K)-Hzs(I,J,K));
         else
         end          
end
end
end
for I=1:Nx+1;
for J=1:Ny;    
for K=1:Nz;          
         if (Material(I,J,K)==4 && Material(I,J,K+1)==0) %  Reversa                                          
         Eys(I,J,K) = CexeH2O*Eys(I,J,K) + CeicH2O*Eyi(I,J,K) + CeipH2O*Eyianterior(I,J,K) + edzH2O*(Hxs(I,J,K)-Hxs(I,J,K-1))- edxH2O*(Hzs(I,J,K)-Hzs(I-1,J,K));
         Eys(I,J,K+1) = CexeH2O*Eys(I,J,K+1) +CeicH2O*Eyi(I,J,K+1) + CeipH2O*Eyianterior(I,J,K) + edzH2O*(Hxs(I,J,K+1)-Hxs(I,J,K))- edxH2O*(Hzs(I,J,K+1)-Hzs(I-1,J,K+1));
         else
         end
         if (Material(I,J,K)==4 && Material(I,J,K+1)==4) %  vol                                          
         Eys(I,J,K) = CexeH2O*Eys(I,J,K) +CeicH2O*Eyi(I,J,K) + CeipH2O*Eyianterior(I,J,K) + edzH2O*(Hxs(I,J,K)-Hxs(I,J,K-1))- edxH2O*(Hzs(I,J,K)-Hzs(I-1,J,K));
         else
         end                  
         if (Material(I,J,K)==0 && Material(I,J,K+1)==4) %  Directa                                          
         Eys(I,J,K) = CexeH2O*Eys(I,J,K) +CeicH2O*Eyi(I,J,K) + CeipH2O*Eyianterior(I,J,K) + edzH2O*(Hxs(I,J,K)-Hxs(I,J,K-1))- edxH2O*(Hzs(I,J,K)-Hzs(I-1,J,K));
         %Eys(I,J,K+1) = CexeH2O*Eys(I,J,K+1) +CeicH2O*Eyi(I,J,K+1) + CeipH2O*Eyianterior(I,J,K) + edzH2O*(Hxs(I,J,K+1)-Hxs(I,J,K))- edxH2O*(Hzs(I,J,K+1)-Hzs(I-1,J,K+1));
         else
         end          
end
end
end





for I=1:Nx+1;%Ez calculated in the volume and in the six faces
for J=1:Ny;    
for K=1:Nz;  
         if (Material(I,J,K)==4 && Material(I,J+1,K)==0) %  Reversa                                          
         Ezs(I,J,K) = CexeH2O*Ezs(I,J,K) + CeicH2O*Ezi(I,J,K) + CeipH2O*Ezianterior(I,J,K) + edxH2O*(Hys(I,J,K)-Hys(I-1,J,K))- edyH2O*(Hxs(I,J,K)-Hxs(I,J-1,K));
         Ezs(I,J+1,K) = CexeH2O*Ezs(I,J+1,K) + CeicH2O*Ezi(I,J+1,K)+ CeipH2O*Ezianterior(I,J,K) + edxH2O*(Hys(I,J+1,K)-Hys(I-1,J+1,K))- edyH2O*(Hxs(I,J+1,K)-Hxs(I,J,K));

%         Ezs(I,J,K) = -Ezi(I,J,K);
%         Ezs(I,J+1,K) = -Ezi(I,J+1,K);
         else
         end
         if (Material(I,J,K)==4 && Material(I,J+1,K)==4) %  vol                                          
         Ezs(I,J,K) = CexeH2O*Ezs(I,J,K) + CeicH2O*Ezi(I,J,K)+ CeipH2O*Ezianterior(I,J,K) + edxH2O*(Hys(I,J,K)-Hys(I-1,J,K))- edyH2O*(Hxs(I,J,K)-Hxs(I,J-1,K));
         %Ezs(I,J+1,K) = CexeH2O*Ezs(I,J+1,K) + CeicH2O*Ezi(I,J+1,K)+ CeipH2O*Ezianterior(I,J,K) + edxH2O*(Hys(I,J+1,K)-Hys(I-1,J+1,K))- edyH2O*(Hxs(I,J+1,K)-Hxs(I,J,K));
         else
         end                  
         if (Material(I,J,K)==0 && Material(I,J+1,K)==4) %  Directa                                          
         Ezs(I,J,K) = CexeH2O*Ezs(I,J,K) +CeicH2O*Ezi(I,J,K)+ CeipH2O*Ezianterior(I,J,K) + edxH2O*(Hys(I,J,K)-Hys(I-1,J,K))- edyH2O*(Hxs(I,J,K)-Hxs(I,J-1,K));
         %Ezs(I,J+1,K) = CexeH2O*Ezs(I,J+1,K) +CeicH2O*Ezi(I,J+1,K)+ CeipH2O*Ezianterior(I,J,K) + edxH2O*(Hys(I,J+1,K)-Hys(I-1,J+1,K))- edyH2O*(Hxs(I,J+1,K)-Hxs(I,J,K));
         else
         end        
end
end
end
for I=1:Nx;
for J=1:Ny+1;    
for K=1:Nz;    
         if (Material(I,J,K)==4 && Material(I+1,J,K)==0) %  Reversa                                          
         Ezs(I,J,K) = CexeH2O*Ezs(I,J,K) + CeicH2O*Ezi(I,J,K)+ CeipH2O*Ezianterior(I,J,K) + edxH2O*(Hys(I,J,K)-Hys(I-1,J,K))- edyH2O*(Hxs(I,J,K)-Hxs(I,J-1,K));
         Ezs(I+1,J,K) = CexeH2O*Ezs(I+1,J,K) + CeicH2O*Ezi(I+1,J,K)+ CeipH2O*Ezianterior(I,J,K) + edxH2O*(Hys(I+1,J,K)-Hys(I,J,K))- edyH2O*(Hxs(I+1,J,K)-Hxs(I+1,J-1,K));

%         Ezs(I,J,K) = -Ezi(I,J,K);
%         Ezs(I+1,J,K) = -Ezi(I+1,J,K);
         else
         end
         if (Material(I,J,K)==4 && Material(I+1,J,K)==4) %  vol                                          
         Ezs(I,J,K) = CexeH2O*Ezs(I,J,K) + CeicH2O*Ezi(I,J,K)+ CeipH2O*Ezianterior(I,J,K) +edxH2O*(Hys(I,J,K)-Hys(I-1,J,K))- edyH2O*(Hxs(I,J,K)-Hxs(I,J-1,K));
         %Ezs(I+1,J,K) = CexeH2O*Ezs(I+1,J,K) + CeicH2O*Ezi(I+1,J,K)+ CeipH2O*Ezianterior(I,J,K) + edxH2O*(Hys(I+1,J,K)-Hys(I,J,K))- edyH2O*(Hxs(I+1,J,K)-Hxs(I+1,J-1,K));
         else
         end                  
         if (Material(I,J,K)==0 && Material(I+1,J,K)==4) %  Directa                                          
         Ezs(I,J,K) = CexeH2O*Ezs(I,J,K) + CeicH2O*Ezi(I,J,K) + CeipH2O*Ezianterior(I,J,K) + edxH2O*(Hys(I,J,K)-Hys(I-1,J,K))- edyH2O*(Hxs(I,J,K)-Hxs(I,J-1,K));
         %Ezs(I+1,J,K) = CexeH2O*Ezs(I+1,J,K) + CeicH2O*Ezi(I+1,J,K)+ CeipH2O*Ezianterior(I,J,K) + edxH2O*(Hys(I+1,J,K)-Hys(I,J,K))- edyH2O*(Hxs(I+1,J,K)-Hxs(I+1,J-1,K));
         else
         end        
end
end
end




%% Scattered Electric Field Update in PECs


for I=1:Nx;%Ex calculated in the volume and in the six faces
for J=1:Ny;    
for K=1:Nz+1;                
         if ((Material(I,J,K)==2 && Material(I,J+1,K)==0)) %  Reversa                                          
         Exs(I,J,K) = CeicCu*Exi(I,J,K);
         Exs(I,J+1,K) = CeicCu*Exi(I,J+1,K);
         else
         end
         if (Material(I,J,K)==2 && Material(I,J+1,K)==2) % vol                                          
         Exs(I,J,K) = CeicCu*Exi(I,J,K);
         else
         end                  
         if (Material(I,J,K)==0 && Material(I,J+1,K)==2) %  Directa                                          
         Exs(I,J,K) = CeicCu*Exi(I,J,K);
         Exs(I,J+1,K) = CeicCu*Exi(I,J+1,K);
         else
         end        
end
end
end
for I=1:Nx;
for J=1:Ny+1;    
for K=1:Nz;  
         if (Material(I,J,K)==2 && Material(I,J,K+1)==0) %  Reversa                                          
         Exs(I,J,K) = CeicCu*Exi(I,J,K);
         Exs(I,J,K+1) = CeicCu*Exi(I,J,K+1);
         else
         end
         if (Material(I,J,K)==2 && Material(I,J,K+1)==2) %  vol                                          
         Exs(I,J,K) = CeicCu*Exi(I,J,K);
         else
         end        
         if (Material(I,J,K)==0 && Material(I,J,K+1)==2) %  Directa                                          
         Exs(I,J,K) = CeicCu*Exi(I,J,K);
         Exs(I,J,K+1) = CeicCu*Exi(I,J,K+1);
         else
         end        
end
end
end




for I=1:Nx;%Ey calculated in the volume and in the six faces
for J=1:Ny;    
for K=1:Nz+1;                        
         if (Material(I,J,K)==2 && Material(I+1,J,K)==0) %  Reversa                                          
         Eys(I,J,K) = CeicCu*Eyi(I,J,K);
         Eys(I+1,J,K) = CeicCu*Eyi(I+1,J,K);
         else
         end
         if (Material(I,J,K)==2 && Material(I+1,J,K)==2) %  vol                                          
         Eys(I,J,K) = CeicCu*Eyi(I,J,K);
         else
         end                  
         if (Material(I,J,K)==0 && Material(I+1,J,K)==2) %  Directa                                          
         Eys(I,J,K) = CeicCu*Eyi(I,J,K);
         Eys(I+1,J,K) = CeicCu*Eyi(I+1,J,K);
         else
         end          
end
end
end
for I=1:Nx+1;
for J=1:Ny;    
for K=1:Nz;          
         if (Material(I,J,K)==2 && Material(I,J,K+1)==0) %  Reversa                                          
         Eys(I,J,K) = CeicCu*Eyi(I,J,K);
         Eys(I,J,K+1) = CeicCu*Eyi(I,J,K+1);
         else
         end
         if (Material(I,J,K)==2 && Material(I,J,K+1)==2) %  vol                                          
         Eys(I,J,K) = CeicCu*Eyi(I,J,K);
         else
         end                  
         if (Material(I,J,K)==0 && Material(I,J,K+1)==2) %  Directa                                          
         Eys(I,J,K) = CeicCu*Eyi(I,J,K);
         Eys(I,J,K+1) = CeicCu*Eyi(I,J,K+1);
         else
         end          
end
end
end


for I=1:Nx+1;%Ez calculated in the volume and in the six faces
for J=1:Ny;    
for K=1:Nz;  
         if (Material(I,J,K)==2 && Material(I,J+1,K)==0) %  Reversa                                          
         Ezs(I,J,K) = CeicCu*Ezi(I,J,K);
         Ezs(I,J+1,K) = CeicCu*Ezi(I,J+1,K);
         else
         end
         if (Material(I,J,K)==2 && Material(I,J+1,K)==2) %  vol                                          
         Ezs(I,J,K) = CeicCu*Ezi(I,J,K);
         else
         end                  
         if (Material(I,J,K)==0 && Material(I,J+1,K)==2) %  Directa                                          
         Ezs(I,J,K) = CeicCu*Ezi(I,J,K);
         Ezs(I,J+1,K) = CeicCu*Ezi(I,J+1,K);
         else
         end        
end
end
end
for I=1:Nx;
for J=1:Ny+1;    
for K=1:Nz;    
         if (Material(I,J,K)==2 && Material(I+1,J,K)==0) %  Reversa                                          
         Ezs(I,J,K) = CeicCu*Ezi(I,J,K);
         Ezs(I+1,J,K) = CeicCu*Ezi(I+1,J,K);
         else
         end
         if (Material(I,J,K)==2 && Material(I+1,J,K)==2) %  vol                                          
         Ezs(I,J,K) = CeicCu*Ezi(I,J,K);
         else
         end                  
         if (Material(I,J,K)==0 && Material(I+1,J,K)==2) %  Directa                                          
         Ezs(I,J,K) = CeicCu*Ezi(I,J,K);
         Ezs(I+1,J,K) = CeicCu*Ezi(I+1,J,K);
         else
         end        
end
end
end


%% Total Magnetic and Electric Field Update in the whole volume

Hxtotal(:,:,:)=Hxs(:,:,:)+Hxi(:,:,:);

Extotal(:,:,:)=Exs(:,:,:);
%Eytotal(:,:,:)=Eys(:,:,:);
Eztotal(:,:,:)=Ezs(:,:,:)+Ezi(:,:,:);
%Hxtotal(:,:,:)=Hxs(:,:,:)+Hxi(:,:,:);

Exgraficable(:,:)=Extotal(:,round(Ny/2),(1:end-1));
Ezgraficable(:,:)=Eztotal((1:end-1),round(Ny/2),:);

Ezt(t) = Eztotal(round(Nx/2),round(2*Ny/3)+round(Lf/4),round(Nz/2)); % This is the electric field in the given coordinates inside the water as a function of the time, with the matching metamaterial present in the interface air water.
%Ezt(t) = Eztotal(round(Nx/2),10,round(Nz/2));
%Ezgraficable(:,:)=Eztotal(round(Nx/2),:,:);

%% Point Source
%Ezs(round(Nx/2)+1,round(Ny/2)+1,round(Nz/2)) =...
%sin(2*pi*f0*dt*t)./(1+exp(-.2*(t-60)));
%% Plane Source

Ezi(:,1,:) =...
sin(2*pi*f0*dt*(t-600))./((1+exp(-.025*(t-300)))*(1+exp(.025*(t-900))));
Ezi0(t) = Ezi(round(Nx/2),round(2*Ny/3)+round(Lf/4),round(Nz/2)); % This is the electric field in the given coordinates inside the water as a function of the time, with the without the matching metamaterial present in the interface air water.
%Ezi0(t) = Ezi(round(Nx/2),10,round(Nz/2));
%% Plotting

%if t == 1
%[X,Y,Z] = meshgrid(y,x,z(1:end-1)); % for plotting Ez
%end
%slice(X,Y,Z,Ezs,y(end)/2,x(end)/2,z(end)/2);
%axis([0 y(end) 0 x(end) 0 z(end)]);
%view([0.3 1 0.3]); caxis([-.15 .15]); shading interp; drawnow;






if (t > (nt*6.8/8));

    for P=1:(Ny);
    if ( Eztotal(round(Nx/2),P,round(Nz/2)) >=  MaxEz(P) );  
    MaxEz(P) = Eztotal(round(Nx/2),P,round(Nz/2));
    end
    if ( MinEz(P) >= Eztotal(round(Nx/2),P,round(Nz/2)));  
    MinEz(P) = Eztotal(round(Nx/2),P,round(Nz/2));
    end    
    end
    
end 

zmax = 2;
ymax = 2;
xmax = 2;

xlim('manual');
xlim([-xmax xmax]);
ylim('manual');
ylim([-ymax ymax]);
zlim('manual');
zlim([-zmax zmax]);
grafico0 = Eztotal(round(Nx/2),:,round(Nz/2));
grafico1 = Ezs(round(Nx/2),:,round(Nz/2));
grafico2 = MaxEz;
grafico3 = Ezi0;
grafico4 = Ezt;
subplot(4,1,1); plot(grafico0);figure(gcf);
%title('Income')
subplot(4,1,2); plot(grafico3);figure(gcf);
%title('Income')
subplot(4,1,3); plot(grafico4);figure(gcf);
%title('Outgo')

%plot(Ezs(round(Nx/2),:,round(Nz/2))); figure(gcf);
%Evacionoacoplado(t)=(Ezi(1+round(Nx/2),round(275*1.25),round(Nz/2)));
%Evacioacoplado(t)=(Eztotal(1+round(Nx/2),round(202*1.25),round(Nz/2)));
%Edielacoplado(t)=(Eztotal(1+round(Nx/2),round(2*Ny/3-50),round(Nz/2)));
%Edielnoacoplado(t)=(Ezi(1+round(Nx/2),round(2*Ny/3-50),round(Nz/2)));

Evacionoacoplado(t)=(Ezi(1+round(Nx/2),round(269),round(Nz/2)));
Evacioacoplado(t)=(Eztotal(1+round(Nx/2),round(449),round(Nz/2)));
Edielacoplado(t)=(Eztotal(1+round(Nx/2),round(2*Ny/3+20),round(Nz/2)));
Edielnoacoplado(t)=(Ezi(1+round(Nx/2),round(2*Ny/3+20),round(Nz/2)));
iteracion = t

%Auxiliar matrices for optional graphics

if(rem(t,20)==0);
save('Cubo1216vpaper');
end
toc

end
for I=1:Nx;
    for K=1:Nz;
A(I,K) = Eztotal(I,round(2*Ny/3-3),K);
G(I,K) = Material(I,round(2*Ny/3-3),K);
    end
end

for J=1:Ny;
    for K=1:Nz;
B(J,K) = Eztotal(round(Nx/2),J,K);
D(J,K) = Material(round(Nx/2),J,K);
    end
end

for I=1:Nx;
    for J=1:Ny;
C(I,J) = Eztotal(I,J,round(Nz/2));
F(I,J) = Material(I,J,round(Nz/2));    
    end
end


clear Transmitida;
clear Incidente;
clear f
clear P1
clear P2
clear P12
clear P22
clear Y
clear Y2
clear Datacruda
clear Datacruda2

Inicio = 1;
reduccion = 1;
ReducSpam = 200;
Mmax = round((nt-Inicio)/reduccion);

Datacruda=Ezt(Inicio:end);
Transmitida((Mmax))=0;
for I=1:round(Mmax);
Transmitida(I)=Datacruda(I*reduccion);
end

T = reduccion*dt;
Fs = 1/T;


L = Mmax;
%t = (0:L-1)*T;
%S = (0.7)*sin(2*pi*50*t);
S = Transmitida;
Y = fft(S,round(nt*ReducSpam));
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1)=2*P1(2:end-1);
f= Fs*(0:(L/2))/(L*ReducSpam);
%plot(f,P1,'green');
subplot(4,1,4); plot(f,P1,'green');figure(gcf);
hold;




Datacruda2=Ezi0(Inicio:end);
Incidente(round(Mmax))=0;
for I=1:round(Mmax);
Incidente(I)=Datacruda2(I*reduccion);
end


%L = 1500;
%L = Mmax;
%t = (0:L-1)*T;
%S = (0.7)*sin(2*pi*50*t);
S2 = Incidente;
Y2 = fft(S2,nt*ReducSpam);
P22 = abs(Y2/L);
P12 = P22(1:L/2+1);
P12(2:end-1)=2*P12(2:end-1);
f= Fs*(0:(L/2))/(L*ReducSpam);
subplot(4,1,4); plot(f,P12,'red');figure(gcf);

% Once this code is run, please run the file Eficiency_vs_f_paper.m if you
% want to reproduce the paper`s Fig. 13.


