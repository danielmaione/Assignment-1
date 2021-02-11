%Daniel Maione
%101076393
%The code will be generated below using the recomendation of Professor Smy
%who recomends in the assignment to use arrays of positions and 
%velocities (Px, Py, Vx, Vy)

set(0,'DefaultFigureWindowStyle','docked')
set(0,'defaultaxesfontsize',20)
set(0,'defaultaxesfontname','Times New Roman')
set(0,'DefaultLineLineWidth',2)

clear all;
close all;

m0 = 9.10938356e-31;
m = 0.26*m0;
T = 300;
k = 1.38064852e-23;
h = 100e-9;
l = 200e-9;
it=2000;
vth = sqrt(2*k*T/m)
Constanttime=0.2e-12;
showmovie=0;

%the next step is to find the mean free path with the equation below
mfp = vth*Constanttime;


tstep=1e-15;

%now probability of scatter can be calculated
P=1-exp(-tstep/Constanttime);

%next is to declare vaiables for electron and number of electron
nelectron=10000;
plotelectrons=12;

%declare a matrix for the electrons that will determine the x and y
%positions and velocities
position=zeros(nelectron,4);



for i=1:nelectron
    
%initial velocity is set up with maxwell boltzman function
nx = randn*(vth/sqrt(2));
ny = randn*(vth/sqrt(2));

%initial positions
px = l*rand; % Sets the initial x positions as a vector of randomly selected numbers over the length of the simulation
py= h*rand;

%fill matrix with values
position(i,:)=[px py nx ny];
end


%have to make a plot with temperatures and trajectories at each iteration, use matrix
Temp=zeros(it,1);
Traj=zeros(it,plotelectrons*2); %has to be *2 since it needs to contain x and y vectors

%now that everything is initalized, iterate through matrixs to get new
%velocities and positions
for i=1:it
    %to get the new distance of the electrons, must add the previous to the
    %velocity times the time the elctrons moved (ex m/s*s=m, new distance)
    position(:,1:2)=position(:,1:2)+tstep.*position(:,3:4);
   
    %first check to see if electrons are out of boundaries and then
    %reassign them based upon if crosses y reflects (mirror), if crosses x goes
    %through and continues (pac-man)
    overbpx=position(:,1)>l;
    position(overbpx,1)=position(overbpx,1)-l; %if greater meets on other side
    
    underbpx=position(:,1)<0;
    position(underbpx,1)=position(underbpx,1)+l; %if less meets on the other side
    
    overbpy=position(:,2)>h;
    position(overbpy,2)=h*2-position(overbpy,2);%reflects from top
    position(overbpy,4)=-position(overbpy,4);%reverses velocity so electron goes the other way
    
    underbpy=position(:,2)<0;
    position(underbpy,2)=-position(underbpy,2);%reflects from bottom    
    position(underbpy,4)=-position(underbpy,4);%reverses velocity so electron goes the other way
        
    %Update the trajectory matrix for each velocity matrix
    for v=1:plotelectrons
        Traj(i,(v*2):((v*2)+1))=position(v,1:2);
    end
    
    
    
    
    %Now that positions and velocitities are corrected, trajectores can
    %plotted
    if showmovie==1
    figure(1);
    trajcolors=hsv(plotelectrons);
    plot(position(1:plotelectrons,1),position(1:plotelectrons,2),'o')
    title('Position and Trajectories of Electrons')
    xlabel('Position in X (m)')
    ylabel('Position in Y (m)')
    grid on
    axis([0 l 0 h]) 
    pause(0.05)
    end
    
    %average velocity
    avgvelocity=(sum(position(:,3).^2)+sum(position(:,4).^2))/nelectron;
    
    %have to update temperaturE matrix
    %Temperature is found using the equation 1/2*mn*vavg/k   
    
    Temp(i) = (avgvelocity*0.5*m)/k;
    
    
    
    
end
    
%plotting temperatures
%Variables used for plotting
    time=0:it-1;
    minT=min(Temp);
    maxT=max(Temp);
    
    figure(2)
    hold off
    plot(time*tstep,Temp);
    axis([0 tstep*it minT-1 maxT+1]);
    title ('Temperature of Semiconductor')
    xlabel ('Time (s)')
    ylabel ('Temperature (K)')
    grid on
    
    %crete plot for trajectories
    figure(1);
    
    title('Electron Trajectories');
    xlabel('Distance in X (m)');
    ylabel('Distance in Y (m)');
    axis([0 l 0 h]);
    hold on;
    for i=1:plotelectrons
        plot(Traj(:,i*2), Traj(:,i*2+1), '.');
    end
    grid on
    
    
%Question 2

%first Assign a random velocity to each of the particles at the start. To do this
%use a Maxwell-Boltzmann distribution for each velocity component.
%From the notes on Monte Carlo modeling it was given that Velocity for 
%each direction is distribution is Gaussian with std of sqrt(kT/m) This can
%be modeled using the makedist function in matlab

maxboltz=makedist('Normal','mu',0,'sigma',sqrt(k*T/m));

%now that the distribution is found for velocities, the same process can be
%followed as question one to fill in the position matrix

for i=1:nelectron
    
%initial velocity is set up with maxwell boltzman function
nx = random(maxboltz);
ny = random(maxboltz);

%initial positions
px = l*rand; % Sets the initial x positions as a vector of randomly selected numbers over the length of the simulation
py= h*rand;

%fill matrix with values
position(i,:)=[px py nx ny];
end

%the next step given in the assignment is to find the histogram, to do this
%the same steps as the previous question is followed just with implementing
%probability of scattering

for i=1:it
    %to get the new distance of the electrons, must add the previous to the
    %velocity times the time the elctrons moved (ex m/s*s=m, new distance)
    position(:,1:2)=position(:,1:2)+tstep.*position(:,3:4);
   
    %first check to see if electrons are out of boundaries and then
    %reassign them based upon if crosses y reflects (mirror), if crosses x goes
    %through and continues (pac-man)
    overbpx=position(:,1)>l;
    position(overbpx,1)=position(overbpx,1)-l; %if greater meets on other side
    
    underbpx=position(:,1)<0;
    position(underbpx,1)=position(underbpx,1)+l; %if less meets on the other side
    
    overbpy=position(:,2)>h;
    position(overbpy,2)=h*2-position(overbpy,2);%reflects from top
    position(overbpy,4)=-position(overbpy,4);%reverses velocity so electron goes the other way
    
    underbpy=position(:,2)<0;
    position(underbpy,2)=-position(underbpy,2);%reflects from bottom    
    position(underbpy,4)=-position(underbpy,4);%reverses velocity so electron goes the other way
        
    for sc=1:nelectron
    if P>rand
        %Now new velocities need to be assigned
        nx = random(maxboltz);
        ny = random(maxboltz);
        position(sc,3)=nx;
        position(sc,4)=ny;
    end
    end
    
    %Update the trajectory matrix for each velocity matrix
    for v=1:plotelectrons
        Traj(i,(v*2):((v*2)+1))=position(v,1:2);
    end

 %Now that positions and velocitities are corrected, trajectores can
    %plotted
    if showmovie==1
    figure(3);
    trajcolors=hsv(plotelectrons);
    plot(position(1:plotelectrons,1),position(1:plotelectrons,2),'o')
    title('Position and Trajectories of Electrons')
    xlabel('Position in X (m)')
    ylabel('Position in Y (m)')
    grid on
    axis([0 l 0 h]) 
    
    
    %average velocity
    avgvelocity=(sum(position(:,3).^2)+sum(position(:,4).^2))/nelectron;
    
    %have to update temperaturE matrix
    %Temperature is found using the equation 1/2*mn*vavg/k   
    
    Temp(i) = (avgvelocity*0.5*m)/k;
    
    %now that everything is calculated histogram can be made
    %a histogram needs to be made of the x and y velocities but the
    %finction in matlab only allows the use of one variable for velocity so
    %an average velocity is calculated
    figure(5)
    hist(avgvelocity);
    title('Histogram Of MFP')
    xlabel('Electron Speed (m/s)')
    ylabel('Total Number of Electrons')
    pause(0.05)
    end
end

    
	figure(3)
    
    title('Electron Trajectories');
    xlabel('Distance in X (m)');
    ylabel('Distance in Y (m)');
    axis([0 l 0 h]);
    hold on;
    for i=1:plotelectrons
        plot(Traj(:,i*2), Traj(:,i*2+1), '.');
    end
    grid on
    
    %plotting temperatures
%Variables used for plotting
    time=0:it-1;
    minT=min(Temp);
    maxT=max(Temp);
    
    figure(4)
    hold off
    plot(time*tstep,Temp);
    axis([0 tstep*it minT-1 maxT+1]);
    title ('Temperature of Semiconductor')
    xlabel ('Time (s)')
    ylabel ('Temperature (K)')
    grid on
    
    avgvelocity=sqrt((position(:,3).^2)+sum(position(:,4).^2));
    figure(5)
    hist(avgvelocity);
    title('Histogram Of MFP')
    xlabel('Electron Speed (m/s)')
    ylabel('Total Number of Electrons')

    %Question 3
    %First start by defining the boundaries of the boxes
    
    %Variables used in Boundaries
    minx=80e-9;
    maxx=120e-9;
    miny=0;
    bmidy=40e-9;
    tmidy=60e-9;
    maxy=100e-9;
    
    %generate intial position of electrons like before
    for i=1:nelectron
    
        %initial velocity is set up with maxwell boltzman function
        nx = random(maxboltz);
        ny = random(maxboltz);

        %initial positions
        px = l*rand; % Sets the initial x positions as a vector of randomly selected numbers over the length of the simulation
        py= h*rand;

        %fill matrix with values
        position(i,:)=[px py nx ny];
        
        %loop to determine if it is in the box and to 
        while((px>minx && px<maxx && py>miny && py<bmidy)||(px>minx && px<maxx && py>tmidy && py<maxy))
            px=l*rand;
            py=h*rand;
            position(i,1:2)=[px py];
        end
    end
    
    %Have to make it so all boundaries capable to be either specular or 
    %diffusive 
    %Choose through selecting a value of 1 for diffusive and 0 for specular
    %Specular is so they bounce off at the same angle they came in at and
    %diffusive is so they bounce off at random angles, this is only
    %required in the y direction and when it hits the boxes
    SpecorDiffuse=0;
    %repeat assigning intital positions as before
    
    for i=1:it
    %to get the new distance of the electrons, must add the previous to the
    %velocity times the time the elctrons moved (ex m/s*s=m, new distance)
    position(:,1:2)=position(:,1:2)+tstep.*position(:,3:4);
   
    %first check to see if electrons are out of boundaries and then
    %reassign them based upon if crosses y reflects (mirror), if crosses x goes
    %through and continues (pac-man)
    overbpx=position(:,1)>l;
    position(overbpx,1)=position(overbpx,1)-l; %if greater meets on other side
    
    underbpx=position(:,1)<0;
    position(underbpx,1)=position(underbpx,1)+l; %if less meets on the other side
    
    overbpy=position(:,2)>h;
    if SpecorDiffuse==0
        position(overbpy,2)=h*2-position(overbpy,2);%reflects from top
        position(overbpy,4)=-position(overbpy,4);%reverses velocity so electron goes the other way
    
    else
        %The formula for calculating a new angle for the velocity and then
        %the new component in matrix "position" is Vx=v*cos(angle) and
        %Vy=V*sin(angle), this is taken from MC slides
        angle=2*pi*rand(nelectron,1);
        vavgDiff=sqrt((position(overbpy,3).^2)+(position(overby,4).^2));
        position(overbpy,3)=vavgDiff.*cos(angle);
        position(overbpy,4)=vavDiff.*sin(angle);
    end
        
    
    underbpy=position(:,2)<0;
     if SpecorDiffuse==0
        position(underbpy,2)=-position(underbpy,2);%reflects from bottom    
        position(underbpy,4)=-position(underbpy,4);%reverses velocity so electron goes the other way
    
    else
        %The formula for calculating a new angle for the velocity and then
        %the new component in matrix "position" is Vx=v*cos(angle) and
        %Vy=V*sin(angle), this is taken from MC slides
        angle=2*pi*rand(nelectron,1);
        vavgDiff=sqrt((position(underbpy,3).^2)+(position(underby,4).^2));
        position(underbpy,3)=vavgDiff.*cos(angle);
        position(underbpy,4)=-abs(vavDiff.*sin(angle));
     end
     %now check if the particles are in the box
     
     for sc=1:nelectron
        if P>rand
            %Now new velocities need to be assigned
            nx = random(maxboltz);
            ny = random(maxboltz);
            position(sc,3)=nx;
            position(sc,4)=ny;
        end
     end
   
     
     for b=1:it
        %loop for first box
        if ((position(b,1)>minx)&&(position(b,1)<maxx)&&(position(b,2)>miny)&&(position(b,2)<bmidy))            
            if SpecorDiffuse==0
                %only need to change the velocity since the points did not
                %start in the box
                position(b,3)=-position(b,3); 
                position(b,4)=-position(b,4);  
            else
                %The formula for calculating a new angle for the velocity and then
                %the new component in matrix "position" is Vx=v*cos(angle) and
                %Vy=V*sin(angle), this is taken from MC slides
                %have to determine if the velocity for x and y componenets
                %were positive or negative to give proper velocity
                signx=sign(position(b,3)); 
                signy=sign(position(b,4));
                angle=2*pi*rand(nelectron,1);
                vavgDiff=sqrt((position(b,3).^2)+(position(b,4).^2));
                position(b,3)=signx*abs(vavgDiff.*cos(angle));
                position(b,4)=signy*abs(vavDiff.*sin(angle));
            end
        elseif ((position(b,1)>minx)&&(position(b,1)<maxx)&&(position(b,2)>tmidy)&&(position(b,2)<maxy))            
            if SpecorDiffuse==0
                %only need to change the velocity since the points did not
                %start in the box
                position(b,3)=-position(b,3); 
                position(b,4)=-position(b,4);  
            else
                %The formula for calculating a new angle for the velocity and then
                %the new component in matrix "position" is Vx=v*cos(angle) and
                %Vy=V*sin(angle), this is taken from MC slides
                signx=sign(position(b,3)); 
                signy=sign(position(b,4));
                angle=2*pi*rand(nelectron,1);
                vavgDiff=sqrt((position(b,3).^2)+(position(b,4).^2));
                position(b,3)=signx*abs(vavgDiff.*cos(angle));
                position(b,4)=signy*abs(vavDiff.*sin(angle));
            end
        end
     end
    
            
       
     

    %Update the trajectory matrix for each velocity matrix
    for v=1:plotelectrons
        Traj(i,(v*2):((v*2)+1))=position(v,1:2);
    end
       
     %Now that positions and velocitities are corrected, trajectores can
    %plotted
    if showmovie==1
    figure(6)
    trajcolors=hsv(plotelectrons);
    plot(position(1:plotelectrons,1),position(1:plotelectrons,2),'o')
    hold on
    title('Position and Trajectories of Electrons')
    xlabel('Position in X (m)')
    ylabel('Position in Y (m)')
    grid on
    axis([0 l 0 h]) 
    rectCnt=100e-9;
    rectDelta=20e-9;
    ylimone([0,40e-9])
    ylimtwo([60e-9,100e-9])
    
    rectX=rectCnt+rectDelta*[-1,1];
    rectYone=ylimone;
    rectYtwo=ylimtwo;
    %pch1=patch(sp(1),rectX([1,2,2,1]),rectYone([1 1 2 2]),'r');
    %pch2=patch(sp(1),rectX([1,2,2,1]),rectYtwo([1 1 2 2]),'r');
    
    
    
    
    %average velocity
    avgvelocity=(sum(position(:,3).^2)+sum(position(:,4).^2))/nelectron;
    
    %have to update temperaturE matrix
    %Temperature is found using the equation 1/2*mn*vavg/k   
    
    Temp(i) = (avgvelocity*0.5*m)/k;
    
    %now that everything is calculated histogram can be made
    %a histogram needs to be made of the x and y velocities but the
    %finction in matlab only allows the use of one variable for velocity so
    %an average velocity is calculated
    figure(8)
    hist(avgvelocity);
    title('Histogram Of MFP')
    xlabel('Electron Speed (m/s)')
    ylabel('Total Number of Electrons')
    pause(0.05)
    end
    end

    figure(6)
subplot(2,1,1);
    
    title('Electron Trajectories');
    xlabel('Distance in X (m)');
    ylabel('Distance in Y (m)');
    axis([0 l 0 h]);
    hold on;
    for i=1:plotelectrons
        plot(Traj(:,i*2), Traj(:,i*2+1), '.');
    end
    %pch1=patch(sp(1),rectX([1,2,2,1]),rectYone([1 1 2 2]),'r');
    %pch2=patch(sp(1),rectX([1,2,2,1]),rectYtwo([1 1 2 2]),'r');
    grid on
    
    %plotting temperatures
%Variables used for plotting
    time=0:it-1;
    minT=min(Temp);
    maxT=max(Temp);
    
    figure(7)
    hold off
    plot(time*tstep,Temp);
    axis([0 tstep*it minT-1 maxT+1]);
    title ('Temperature of Semiconductor')
    xlabel ('Time (s)')
    ylabel ('Temperature (K)')
    grid on
    
    vgvelocity=sqrt((position(:,3).^2)+sum(position(:,4).^2));
    figure(8)
    hist(avgvelocity);
    title('Histogram Of MFP')
    xlabel('Electron Speed (m/s)')
    ylabel('Total Number of Electrons')
    
    
    
    

    

    
    
