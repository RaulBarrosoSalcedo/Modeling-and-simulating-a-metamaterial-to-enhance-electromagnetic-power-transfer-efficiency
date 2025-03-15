%% Nathan Zechar Wright State University 2021 - Intro to FDTD
clear all; close all;
%% Define Simulation Based off Source and Wavelength
f0 = 2.45e9; % Frequency of Source [Hertz], do not change this value.
f0att = 2.45e9; % Frequency for attenuation estimation [Hertz], change if you want to calculate attenuation in a frequency that is different to f0=2.45GHz.
Lf = round(720); % Divisions per Wavelength [unitless]
[Lx,Ly,Lz] = deal(0.084/4,5,0.084/4); % Wavelengths x,y,z [unitless] %Son cinco longitudes de onda
nt = round(8000); % Number of time steps [unitless] 
aristareal = 4; %This is the cubic conductor edge lenght "b" in FDTD cells, being "a" equal to 10 FDTD cells. 
% b/a ratio is adjusted fixing "b" in this line. If aristareal = 4, then b/a = 0.4.
arista = round(aristareal+1);
GapY=0;
t = 0;
%% Simulation Performance
single = 0; % Use Single Precision | 1 = yes , 0 = no |
usegpu = 0; % Use gpuArray | 1 = yes , 0 = no |
%% Spatial and Temporal System
e0 = 8.854*10^-12; % Permittivity of vaMetalum [farad/meter]
u0 = 4*pi*10^-7; % Permeability of vaMetalum [henry/meter]
c0 = 1/(e0*u0)^.5; % Speed of light [meter/second]
sgMetal = 3.77e7; % Metal Conductivity;
sgH2O = 2*0;% Water Conductivity;
eMetal = e0;
eH2O = (81)*e0;
eta0= (u0/e0)^.5;
etaH2O= (u0/(eH2O*(1-1i*(sgH2O/(eH2O*2*3.1415923*f0)))))^.5;

etametateorico = (etaH2O*eta0)^.5;
L0 = c0/f0; % Freespace Wavelength [meter]
t0 = 1/f0; % Source Period [second]
[Nx,Ny,Nz] = deal(round(Lx*Lf),round(Ly*Lf),round(Lz*Lf)); % Points in x,y,z [unitless]
x = linspace(0,Lx,Nx+1)*L0; % x vector [meter]
y = linspace(0,Ly,Ny+1)*L0; % y vector [meter]
z = linspace(0,Lz,Nz+1)*L0; % z vector [meter]
[dx,dy,dz] = deal(x(2),y(2),z(2)); % x,y,z increment [meter]
dt = (dx^-2+dy^-2+dz^-2)^-.5/c0*.99; % Time step CFL condition [second]
sgFDTD= (1/dy)*((2*sgMetal/(2*3.1415927*f0att*u0))^.5); 
etaMetal= (u0/(eMetal*(1-1i*(sgFDTD/(eMetal*2*3.1415923*f0att)))))^.5;
%% Initialize Incident Magnetic and Electric Field Vectors and Materials

[Exi] = deal(zeros(Nx,Ny+1,Nz+1)); % Ex cells
[Eyi] = deal(zeros(Nx+1,Ny,Nz+1)); % Ey cells
[Ezi] = deal(zeros(Nx+1,Ny+1,Nz)); % Ez cells
[Hxi] = deal(zeros(Nx+1,Ny,Nz)); % Hx cells
[Hyi] = deal(zeros(Nx,Ny+1,Nz)); % Hy cells
[Hzi] = deal(zeros(Nx,Ny,Nz+1)); % Hz cells
[udx,udy,udz] = deal(dt/(u0*dx),dt/(u0*dy),dt/(u0*dz)); % Hi vacMetalm coeffcients
[edx,edy,edz] = deal(dt/(e0*dx),dt/(e0*dy),dt/(e0*dz)); % Ei vacMetalm coeffcients

[udxs,udys,udzs] = deal(dt/(u0*dx),dt/(u0*dy),dt/(u0*dz)); % Hs vacMetalm coeffcients
[edxs,edys,edzs] = deal(dt/(e0*dx),dt/(e0*dy),dt/(e0*dz)); % Es vacMetalm coeffcients
[edxMetal,edyMetal,edzMetal] = deal(2*dt/((2*eMetal+sgFDTD*dt)*dx),2*dt/((2*eMetal+sgFDTD*dt)*dy),2*dt/((2*eMetal+sgFDTD*dt)*dz)); % Es Metalrl material Coeffcients
[edxH2O,edyH2O,edzH2O] = deal(dt/((eH2O+sgH2O*dt)*dx), dt/((eH2O+sgH2O*dt)*dy), dt/((eH2O+sgH2O*dt)*dz)); % Es Metalrl material Coeffcients
%CexeMetal = (2*e0-sgFDTD*dt)/(2*e0+sgFDTD*dt);
%CexeH2O = (2*eH2O-sgH2O*dt)/(2*eH2O+sgH2O*dt);
%CeicMetal = (2*(e0-e0)-sgFDTD*dt)/(2*e0+sgFDTD*dt);
%CeicH2O = (2*(e0-eH2O)-sgH2O*dt)/(2*e0+sgH2O*dt);
%CeipMetal = -(2*(e0-e0)+sgFDTD*dt)/(2*e0+sgFDTD*dt);
%CeipH2O = -(2*(e0-eH2O)+sgH2O*dt)/(2*e0+sgH2O*dt);

eta0= (u0/e0)^.5;
etaH2O= real(u0/(eH2O*(1-1i*(sgH2O/(eH2O*2*3.1415923*f0)))))^.5;
%etaMetal= real(u0/(e0*(1-1i*(sgFDTD/(e0*2*3.1415923*f0)))))^.5;


CexeMetal = (2*eMetal-sgFDTD*dt)/(2*eMetal+sgFDTD*dt);
CeicMetal = (2*(e0-eMetal)-sgFDTD*dt)/(2*eMetal+sgFDTD*dt);
CeipMetal = -(2*(e0-eMetal)+sgFDTD*dt)/(2*eMetal+sgFDTD*dt);



[Material] = deal(zeros(Nx+1,Ny+1,Nz+1)); % Materials, conductors or vacMetalm

Material(:,Ny+1,:) = 2; % Este es el final del sistema en Y, con PEC
%Material(:,1,:) = 3; % Este es el inicio del sistema en Y, con PMC
%Material(:,Ny,:) = 3; % Este es el final del sistema en Y, con PMC
%Material(:,Ny-1,:) = 3; % Este es el final del sistema en Y, con PMC
Material(1,:,:) = 3; %Pared izquierda
Material(2,:,:) = 3; %Pared izquierda
Material(Nx,:,:) = 3; %Pared derecha
Material(Nx+1,:,:) = 3; %Pared derecha
Material(:,:,1) = 4; %Pared inferior
Material(:,:,2) = 4; %Pared inferior
Material(:,:,Nz-1) = 4; %Pared superior
Material(:,:,Nz) = 4; %Pared superior
Material(:,:,Nz+1) = 4; %Pared superior

Nobj =150; %*********************************************************** Ojo con esto, va sin /5 para una simulacion seria
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
Ezt1(nt)=0;
Ezt2(nt)=0;
MaxEz(Ny)=0;
MinEz(Ny)=0;

Eztlongitudinal(nt,((Ny/CeldaY)+4))=0;

for M=1:Nobj;
for I=1:Nx+1;% Objeto conductor rectangular
    for J=1:Ny+1;
        for K=1:Nz+1;
                %( K > (Nz/2-(Dsepz/2))) && ( K < (Nz/2+(Dsepz/2)) )
            if ( ( I > (1+Nx/2-(Dsepx/2))) && ( I < (2+Nx/2+(Dsepx/2))) && ( J > ((2*Ny/3-Dsepy-GapY)+(M-1)*CeldaY)) && ( J < ((2*Ny/3)-GapY+(M-1)*CeldaY)) && ( K > (Nz/2-(Dsepz/2))) && ( K < (Nz/2+(Dsepz/2)) ));
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
while  t <= nt
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



    

    %% Incident Electric Field Update in Vaccum

for I=1:Nx;
for J=2:Ny;    
for K=2:Nz;

    Exi(I,J,K) = Exi(I,J,K)+ edy*(Hzi(I,J,K)-Hzi(I,J-1,K))- edz*(Hyi(I,J,K)-Hyi(I,J,K-1));

end
end
end

for I=2:Nx;
for J=1:Ny;    
for K=2:Nz;  

   Eyi(I,J,K) = Eyi(I,J,K) + edz*(Hxi(I,J,K)-Hxi(I,J,K-1))- edx*(Hzi(I,J,K)-Hzi(I-1,J,K));

end
end
end    

for I=2:Nx;
for J=2:Ny;    
for K=1:Nz;
   
   Ezi(I,J,K) = Ezi(I,J,K) + edx*(Hyi(I,J,K)-Hyi(I-1,J,K))- edy*(Hxi(I,J,K)-Hxi(I,J-1,K));
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

%% Scattered Electric Field Update in Metamaterial and Vaccum
 


for I=1:Nx;
for J=2:Ny;    
for K=2:Nz;
    
if (Material(I,J,K)==2)
Exs(I,J,K) = CexeMetal*Exs(I,J,K) + edyMetal*(Hzs(I,J,K)-Hzs(I,J-1,K))- edzMetal*(Hys(I,J,K)-Hys(I,J,K-1)) + CeicMetal*Exi(I,J,K) + CeipMetal*Exianterior(I,J,K);
else
Exs(I,J,K) = Exs(I,J,K) + edys*(Hzs(I,J,K)-Hzs(I,J-1,K))- edzs*(Hys(I,J,K)-Hys(I,J,K-1));
end

end
end
end

for I=2:Nx;
for J=1:Ny;    
for K=2:Nz;  

if (Material(I,J,K)==2)
Eys(I,J,K) = CexeMetal*Eys(I,J,K) + edzMetal*(Hxs(I,J,K)-Hxs(I,J,K-1))- edxMetal*(Hzs(I,J,K)-Hzs(I-1,J,K)) + CeicMetal*Eyi(I,J,K) + CeipMetal*Eyianterior(I,J,K);
else
Eys(I,J,K) = Eys(I,J,K) + edzs*(Hxs(I,J,K)-Hxs(I,J,K-1))- edxs*(Hzs(I,J,K)-Hzs(I-1,J,K));
end

end
end
end    

for I=2:Nx;
for J=2:Ny;    
for K=1:Nz;

if (Material(I,J,K)==2)
Ezs(I,J,K) = CexeMetal*Ezs(I,J,K) + edxMetal*(Hys(I,J,K)-Hys(I-1,J,K))- edyMetal*(Hxs(I,J,K)-Hxs(I,J-1,K)) + CeicMetal*Ezi(I,J,K) + CeipMetal*Ezianterior(I,J,K);
else
Ezs(I,J,K) = Ezs(I,J,K) + edxs*(Hys(I,J,K)-Hys(I-1,J,K))- edys*(Hxs(I,J,K)-Hxs(I,J-1,K));
end

end
end
end



%% Scattered Electric Field Update in PECs

for I=1:Nx;%Ex calMetallated in the volume and in the six faces
for J=1:Ny;    
for K=1:Nz+1;                
         if ((Material(I,J,K)==4 && Material(I,J+1,K)==0)) %  Reversa                                          
         Exs(I,J,K) = -Exi(I,J,K);
         Exs(I,J+1,K) = -Exi(I,J+1,K);
         else
         end
         if (Material(I,J,K)==4 && Material(I,J+1,K)==4) % vol                                          
         Exs(I,J,K) = -Exi(I,J,K);
         else
         end                  
         if (Material(I,J,K)==0 && Material(I,J+1,K)==4) %  Directa                                          
         Exs(I,J,K) = -Exi(I,J,K);
         Exs(I,J+1,K) = -Exi(I,J+1,K);
         else
         end        
end
end
end
for I=1:Nx;
for J=1:Ny+1;    
for K=1:Nz;  
         if (Material(I,J,K)==4 && Material(I,J,K+1)==0) %  Reversa                                          
         Exs(I,J,K) = -Exi(I,J,K);
         Exs(I,J,K+1) = -Exi(I,J,K+1);
         else
         end
         if (Material(I,J,K)==4 && Material(I,J,K+1)==4) %  vol                                          
         Exs(I,J,K) = -Exi(I,J,K);
         else
         end        
         if (Material(I,J,K)==0 && Material(I,J,K+1)==4) %  Directa                                          
         Exs(I,J,K) = -Exi(I,J,K);
         Exs(I,J,K+1) = -Exi(I,J,K+1);
         else
         end        
end
end
end




for I=1:Nx;%Ey calMetallated in the volume and in the six faces
for J=1:Ny;    
for K=1:Nz+1;                        
         if (Material(I,J,K)==4 && Material(I+1,J,K)==0) %  Reversa                                          
         Eys(I,J,K) = -Eyi(I,J,K);
         Eys(I+1,J,K) = -Eyi(I+1,J,K);
         else
         end
         if (Material(I,J,K)==4 && Material(I+1,J,K)==4) %  vol                                          
         Eys(I,J,K) = -Eyi(I,J,K);
         else
         end                  
         if (Material(I,J,K)==0 && Material(I+1,J,K)==4) %  Directa                                          
         Eys(I,J,K) = -Eyi(I,J,K);
         Eys(I+1,J,K) = -Eyi(I+1,J,K);
         else
         end          
end
end
end
for I=1:Nx+1;
for J=1:Ny;    
for K=1:Nz;          
         if (Material(I,J,K)==4 && Material(I,J,K+1)==0) %  Reversa                                          
         Eys(I,J,K) = -Eyi(I,J,K);
         Eys(I,J,K+1) = -Eyi(I,J,K+1);
         else
         end
         if (Material(I,J,K)==4 && Material(I,J,K+1)==4) %  vol                                          
         Eys(I,J,K) = -Eyi(I,J,K);
         else
         end                  
         if (Material(I,J,K)==0 && Material(I,J,K+1)==4) %  Directa                                          
         Eys(I,J,K) = -Eyi(I,J,K);
         Eys(I,J,K+1) = -Eyi(I,J,K+1);
         else
         end          
end
end
end


for I=1:Nx+1;%Ez calMetallated in the volume and in the six faces
for J=1:Ny;    
for K=1:Nz;  
         if (Material(I,J,K)==4 && Material(I,J+1,K)==0) %  Reversa                                          
         Ezs(I,J,K) = -Ezi(I,J,K);
         Ezs(I,J+1,K) = -Ezi(I,J+1,K);
         else
         end
         if (Material(I,J,K)==4 && Material(I,J+1,K)==4) %  vol                                          
         Ezs(I,J,K) = -Ezi(I,J,K);
         else
         end                  
         if (Material(I,J,K)==0 && Material(I,J+1,K)==4) %  Directa                                          
         Ezs(I,J,K) = -Ezi(I,J,K);
         Ezs(I,J+1,K) = -Ezi(I,J+1,K);
         else
         end        
end
end
end
for I=1:Nx;
for J=1:Ny+1;    
for K=1:Nz;    
         if (Material(I,J,K)==4 && Material(I+1,J,K)==0) %  Reversa                                          
         Ezs(I,J,K) = -Ezi(I,J,K);
         Ezs(I+1,J,K) = -Ezi(I+1,J,K);
         else
         end
         if (Material(I,J,K)==4 && Material(I+1,J,K)==4) %  vol                                          
         Ezs(I,J,K) = -Ezi(I,J,K);
         else
         end                  
         if (Material(I,J,K)==0 && Material(I+1,J,K)==4) %  Directa                                          
         Ezs(I,J,K) = -Ezi(I,J,K);
         Ezs(I+1,J,K) = -Ezi(I+1,J,K);
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

Ezt1(t) = Eztotal(round(Nx/2),round(2900/2),12); 
Ezt2(t) = Eztotal(round(Nx/2),round(3100/2),12);

Jinicial = 1;
Jfinal = round(Ny/CeldaY-1);
for J=(Jinicial):(Jfinal);
    Eztlongitudinal(t,J)=Eztotal(round(Nx/2),J*CeldaY,12);
end;


%Este punto es importante para ver el campo transmitido
%Cuando el sistema está corriendo en verdad, debe colocarse:
%Eztotal(round(Nx/2),3000,12) ******


%Ezt(t) = Eztotal(round(Nx/2),10,round(Nz/2));
%Ezgraficable(:,:)=Eztotal(round(Nx/2),:,:);

%% Point Source
%Ezs(round(Nx/2)+1,round(Ny/2)+1,round(Nz/2)) =...
%sin(2*pi*f0*dt*t)./(1+exp(-.2*(t-60)));
%% Plane Source

Ezi(:,1,:) =...
sin(2*pi*f0*dt*(t-600))./((1+exp(-.025*(t-300)))*(1+exp(.025*(t-900))));
Ezi0(t) = Ezi(round(Nx/2),1500,12);%Este punto es importante para ver el campo incidente *******
%Cuando el sistema está corriendo en verdad, debe colocarse: Ezi0(t)=Ezi(round(Nx/2),3000,12)
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
grafico3 = Ezt1;
grafico4 = Ezt2;
subplot(5,1,1); plot(grafico0);figure(gcf);
%title('Income')
subplot(5,1,2); plot(grafico3);figure(gcf);
%title('Income')
subplot(5,1,3); plot(grafico4);figure(gcf);
%title('Outgo')

%plot(Ezs(round(Nx/2),:,round(Nz/2))); figure(gcf);  ****
%Evacionoacoplado(t)=(Ezi(1+round(Nx/2),round(275*1.25),round(Nz/2)));  ***
%Evacioacoplado(t)=(Eztotal(1+round(Nx/2),round(202*1.25),round(Nz/2)));
%Edielacoplado(t)=(Eztotal(1+round(Nx/2),round(2*Ny/3-50),round(Nz/2)));
%Edielnoacoplado(t)=(Ezi(1+round(Nx/2),round(2*Ny/3-50),round(Nz/2)));


iteracion = t

if(rem(t,40)==0);
save('alfa0410_2450MHz');%********************************************************************************************
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

Eztaux = Ezt1;
subplot(5,1,2); plot(Eztaux);figure(gcf);
for I=538:nt;
   Ezt1(I)=Eztaux(I-537); 
end    


%Jinicial = 1;
%Jfinal = round(Ny/Celday-1);
%for J=(Jnicial):(Jfinal);
%    Eztlongitudinal(t,J)=Eztotal(round(Nx/2),J*Celday,12);
%end;

tic
Ezmagfasorlong(Jfinal,20000)=0;

for J=Jinicial:Jfinal;

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
clear Ezt1;

Inicio = 1;
reduccion = 1;
ReducSpam = 200;
Mmax = round((nt-Inicio)/reduccion);
Ezt1=Eztlongitudinal(:,J);
Datacruda=Ezt1(Inicio:end);
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
f= Fs*(0:(L/2))/(L*ReducSpam); % vector f, which includes the considered bandwidth in which it must be also included the column that corresponds to f0att in turn.
%plot(f,P1,'green');
K=0;
while K < numel(P1);
K=K+1;
Ezmagfasorlong(J,K)=P1(K); % Ezmagfasorlong is the magnitude of the harmonic 
%field in the bandwidth given by the vector f. Please, select the column of 
%Ezmagfasorlong that, according to the data of the vector f, corresponds to 
%f0att. Such column corresponds to the field standing wave pattern along the 
%propagation direction which helps to estimate the attenuation constant in 
%the metamaterial observing the field decrescent exponential behaviour in it. 
%Nevertheless, the metamaterial is located in the volume Ny > y > (2/3)*Ny.
% The distance taken between the given samples in the examinated column is
% an edge of a metamaterial cell, it means the value of "a".
end
end
toc

save('alfa0410_2450MHz');

