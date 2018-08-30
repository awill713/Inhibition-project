%IMPORTANT CHANGES

%converted I to Ie, designating excitatory input (assoc. w/ Vna = 0). This
%applies to I in script, in runmodel function, and in soma, spine, and
%dendrite functions

%added Ii, in the same formate as Ie, designate inhibitory input. It is
%going to be associated with Vcl = -70.

%Updated the runmodel, soma, spine, and dendrite functions to  take Ii as
%an argument. Accordingly included Ii as an argument when those functions
%are called within the runmodel function. It goes right behind Ie

%Updated the dv terms in each function include the term: ...+ Ii * (Vcl-v)

%%Added Vcl = -70 to each soma, spine, and dendrite function


% clear all
% close all

% shafts = 49;

global tau
tau = .1;%ms
T=100;%length of sim in ms
n=round(T/tau);
v=-60*ones(n,45); u=0*v; S = u;
Ie=u;
Ii=Ie;

%Experimental setup. I{1,x} is experiment X excitatory input, I{2,x} is
%experiment X inhibitory input
I_setup = cell(2,8);
% I_setup{1,1} = zeros(n,62); I_setup{1,1}(150:200,56)=1; I_setup{2,1} = zeros(n,62); I_setup{2,1}(150:200,16)=1;
I_setup{1,1} = zeros(n,totalCompartments); I_setup{1,1}(150:200,1) = 0; I_setup{2,1} = zeros(n,totalCompartments);
I_setup{2,1}(175:225,1) = 1;

Ie = I_setup{1,1};
Ii = I_setup{2,1};

[v S epsp] = runmodel(Ie,Ii,n,connectome,conductanceMat,totalCompartments,shafts);
figure;plot(v(:,1));


function[v,S,epsp] = runmodel(Ie,Ii,n,con,cond,tc,sha)
v=-60*ones(n,tc); u=0*v; S = u;%init vars

%cell in each column is connected to cell in row with value = 1 (ie connect_m = [ 0 1 0; 1 0 1; 0 1 0] 
%cell 1 connections == connect_m(:,1) = cell 2 only, cell 2 connections =
% %connect_(:,2) == cell 1 and cell 3, etc
C=zeros(62,62);
C(1,2) = 1;
C(2,[3 4]) = 1;

for x = 5:22
    C(x-2,x)=1;
end
for x = 23:42
    C(x-20,x) = 1;
end

for x = 43:62
    C(x-40,x) = 1;
end

C=C+triu(C,1)';

% connect_m = C;
% connect_g=C.*25;
% connect_g(1:22,1:22) = C(1:22,1:22).*150;

connect_m = con;
connect_g = cond;

for j = 1:size(v,1)-1
    state = connect_m.*v(j,:)'; %gets the current values of each neuron and puts them in to the connectivity matrix
    Vn = cell(1,size(state,2));
    Vg=Vn;
    for ii = 1:size(state,2)
        Vn{1,ii} = state(state(:,ii)~=0,ii);%stores these values where Vn(1)= connections for cell1, does this before changing any values(ie calculated for all cells at exactly timepointv(j)
        Vg{1,ii}= connect_g(state(:,ii)~=0,ii);
        
    end
    
    
    [v(j+1,1), u(j+1,1), S(j,1)] = simple_som(v(j,1),u(j,1),Ie(j,1),Ii(j,1),Vn{1},Vg{1}); % soma index = 1
        if S(j,1)>0
            v(j,1) = S(j,1);
        end
        
    for ii = 2 : sha %all the shafts
        [v(j+1,ii), u(j+1,ii), S(j,ii)] = simple_shaft(v(j,ii),u(j,ii),Ie(j,ii),Ii(j,ii),Vn{ii},Vg{ii});
        if S(j,ii)>0
            v(j,ii) = S(j,ii);
        end
    end
    
    for ii = (sha+1) : tc %all the dendrites
        [v(j+1,ii), u(j+1,ii), S(j,ii)] = simple_den(v(j,ii),u(j,ii),Ie(j,ii),Ii(j,ii),Vn{ii},Vg{ii});
        if S(j,ii)>0
            v(j,ii) = S(j,ii);
        end
    end
    
end

epsp = max(v(:,1))+60;
end
function[v1,u1,S] = simple_som(v,u,Ie,Ii,Vn,Vg)
C=100;%capacitance
vr=-60;%resting pot
vt =-50;%threshold
k=0.7;%k and b determined by rheobase and input resistance
a=0.01;%recovery time constant for u(potassium current)
b=5;%b>0 is a resonant neuron with a sag current, b<0 is amplifying
c=-55;
d=500;
%g=100;%nS
vpeak=50;
Vna=0;
Vcl = -70;
global tau;
    S=0; 
    if v>=vpeak
        S=vpeak;
        v = c;
        u=u+d;
    end

    dv = tau*(k*(v-vr)*(v-vt)-u+Ie*(Vna-v)+ Ii*(Vcl-v) +(sum(Vg(:,1).*(Vn(:,1)-v))))/C;
    du = tau*a*(b*(v-vr)-u);
    v1=v+dv;
    u1=u+du;
    
    
end

function[v1,u1,S] = simple_den(v,u,Ie,Ii,Vn,Vg)
C=50;
vr=-60;
vt =-50;
k=3;
a=0.01;
b=5;
c=-55;
d=500;
% g=100;
vpeak=50;
Vna=0;
Vcl = -70;
global tau;

    S=0;
    if v>=vpeak
        S=vpeak;
        v = c;
        u=u+d;
    end
    
    dv = tau*(k*(v-vr)*(v-vt)-u+Ie*(Vna-v)+ Ii*(Vcl-v) +(sum(Vg(:,1).*(Vn(:,1)-v))))/C;
    du = tau*a*(b*(v-vr)-u);
    v1=v+dv;
    u1=u+du;
    
    
end

function[v1,u1,S] = simple_shaft(v,u,Ie,Ii,Vn,Vg)
C=50;
vr=-60;
vt =-50;
k=0.7;
a=0.01;
b=5;
c=-55;
d=500;
g=100;
vpeak=50;
Vna=0;
Vcl = -70;
global tau;
    
    S=0;
    if v>=vpeak
        S=vpeak;
        v = c;
        u=u+d;
    end
    
    dv = tau*((vr-v)+Ie*(Vna-v) + Ii*(Vcl-v) +(sum(Vg(:,1).*(Vn(:,1)-v))))/C;
    du = tau*a*(b*(v-vr)-u);
    v1=v+dv;
    u1=u+du;
    
    
end




%%%Aaron's previous stuff

%%%Increasing excitatory input to a single dendritic shaft
% total = 50;
% epspList = zeros(1,total);
% expected = zeros(1,total);
% for x = 1:total
%     Ie = zeros(n,45);
%     Ie(150:200,17)=x;
%     Ii = zeros(n,45);
%     [v,S,epsp] = runmodel(Ie,Ii,n);
%     epspList(x) = epsp-60;
%     if x==1
%         expected(x) = epsp-60;
%     else
%         expected(x) = (expected(1)+60)*x-60;
%     end
%     x
% end
% 
% figure;
% plot(1:total,epspList,'.');
% hold on;
% plot(1:total,expected,'*');

%%%Increasing excitatory input to a single dendritic spine
% total = 15;
% epspList = zeros(1,total);
% expected = zeros(1,total);
% for x = 1:total
%     Ie = zeros(n,45);
%     Ie(150:1000,45)=x;
%     Ii = zeros(n,45);
%     [v,S,epsp] = runmodel(Ie,Ii,n);
%     epspList(x) = epsp-60;
%     if x==1
%         expected(x) = epsp-60;
%     else
%         expected(x) = (expected(1)+60)*x-60;
%     end
%     x
% end
% 
% figure;
% plot(1:total,epspList,'.');
% hold on;
% plot(1:total,expected,'*');
% 
% figure;
% for comp = 1:45
%     subplot(9,5,comp);
%     plot(v(:,comp));
%     ylim([-80 20]);
%     xlim([100 500]);
% end

%%%Increasing inhibitory input to a single dendritic shaft
% total = 50;
% voltageTraces = zeros(50,n,45);
% for x = 1:50
%     Ie = zeros(n,45);
%     Ii = zeros(n,45);
%     Ii(150:200,17) = x;
%     [v S epsp] = runmodel(Ie,Ii,n);
%     voltageTraces(x,:,:) = v;
%     x
% end
% figure;
% for compartment = 1:45
%     subplot(9,5,compartment);
%     for intensity = 1:total
%         plot(voltageTraces(intensity,:,compartment));
%         hold on;
%     end
%     xlim([100 500]);
% end
% 
% ipspMinima = min(voltageTraces(:,:,1)');
% figure;
% plot(1:total,ipspMinima);
% hold on;
% plot(1:total,(ipspMinima(1)+60)*(1:total)-60,'--');
% title('Increasing inhibitory inputs onto single dendritic shaft');

%%%Increasing inhibitory input to a single dendritic spine
% total = 50;
% voltageTraces = zeros(50,n,45);
% for x = 1:50
%     Ie = zeros(n,45);
%     Ii = zeros(n,45);
%     Ii(150:200,45) = x;
%     [v S epsp] = runmodel(Ie,Ii,n);
%     voltageTraces(x,:,:) = v;
%     x
% end
% figure;
% for compartment = 1:45
%     subplot(9,5,compartment);
%     for intensity = 1:total
%         plot(voltageTraces(intensity,:,compartment));
%         hold on;
%     end
%     xlim([100 500]);
% end
% 
% ipspMinima = min(voltageTraces(:,:,1)');
% figure;
% plot(1:total,ipspMinima);
% hold on;
% plot(1:total,(ipspMinima(1)+60)*(1:total)-60);
% title('Increasing inhibitory inputs onto dendritic spine');

%Increasing inhibitory input to multiple dendritic shafts along the same
%branch
% Ii_setup = cell(1,6);
% Ii_setup{1,1} = zeros(n,45); Ii_setup{1,1}(150:200,17)=1; Ii_setup{1,1}(150:200,16)=1;
% Ii_setup{1,2} = zeros(n,45); Ii_setup{1,2}(150:200,17)=1; Ii_setup{1,2}(150:200,15)=1;
% Ii_setup{1,3} = zeros(n,45); Ii_setup{1,3}(150:200,17)=1; Ii_setup{1,3}(150:200,11)=1;
% Ii_setup{1,4} = zeros(n,45); Ii_setup{1,4}(150:200,17)=1; Ii_setup{1,4}(150:200,3)=1;
% Ii_setup{1,5} = zeros(n,45); Ii_setup{1,5}(150:200,17)=1; Ii_setup{1,5}(150:200,2)=1;
% Ii_setup{1,6} = zeros(n,45); Ii_setup{1,6}(150:200,17)=1; Ii_setup{1,6}(150:200,1)=1;
% 
% figure
% for experiment = 1:6
%     total = 50;
%     voltageTraces = zeros(50,n,45);
%     for x = 1:50
%         Ie = zeros(n,45);
%         Ii = Ii_setup{1,experiment}*x;
%         
%         [v S epsp] = runmodel(Ie,Ii,n);
%         voltageTraces(x,:,:) = v;
%         [experiment x]
%     end
%     
%     ipspMinima = min(voltageTraces(:,:,1)');
%     plot(1:total,ipspMinima);
%     hold on;
% end

%Other
% Ii_setup = cell(1,6);
% Ii_setup{1,1} = zeros(n,45); Ii_setup{1,1}(150:200,15)=1;
% Ii_setup{1,2} = zeros(n,45); Ii_setup{1,2}(150:200,11)=1;
% Ii_setup{1,3} = zeros(n,45); Ii_setup{1,3}(150:200,15)=1; Ii_setup{1,3}(150:200,11)=1;
% Ii_setup{1,4} = zeros(n,45); Ii_setup{1,4}(150:200,15)=1; Ii_setup{1,4}(150:200,4)=1;
% 
% 
% figure
% for experiment = 1:4
%     total = 50;
%     voltageTraces = zeros(50,n,45);
%     for x = 1:50
%         Ie = zeros(n,45);
%         Ii = Ii_setup{1,experiment}*x;
%         
%         [v S epsp] = runmodel(Ie,Ii,n);
%         voltageTraces(x,:,:) = v;
%         [experiment x]
%     end
%     
%     ipspMinima = min(voltageTraces(:,:,1)');
%     plot(1:total,ipspMinima);
%     hold on;
% end

%Increasing inhibitory input to multiple dendritic shafts along different
%branches
% Ii_setup = cell(1,6);
% Ii_setup{1,1} = zeros(n,45); Ii_setup{1,1}(150:200,17)=1; Ii_setup{1,1}(150:200,6)=1;
% Ii_setup{1,2} = zeros(n,45); Ii_setup{1,2}(150:200,17)=1; Ii_setup{1,2}(150:200,5)=1;
% Ii_setup{1,3} = zeros(n,45); Ii_setup{1,3}(150:200,17)=1; Ii_setup{1,3}(150:200,4)=1;
% Ii_setup{1,4} = zeros(n,45); Ii_setup{1,4}(150:200,17)=1; Ii_setup{1,4}(150:200,3)=1;
% Ii_setup{1,5} = zeros(n,45); Ii_setup{1,5}(150:200,17)=1; Ii_setup{1,5}(150:200,2)=1;
% Ii_setup{1,6} = zeros(n,45); Ii_setup{1,6}(150:200,17)=1; Ii_setup{1,6}(150:200,1)=1;
% 
% figure
% for experiment = 1:6
%     total = 50;
%     voltageTraces = zeros(50,n,45);
%     for x = 1:50
%         Ie = zeros(n,45);
%         Ii = Ii_setup{1,experiment}*x;
%         
%         [v S epsp] = runmodel(Ie,Ii,n);
%         voltageTraces(x,:,:) = v;
%         [experiment x]
%     end
%     
%     ipspMinima = min(voltageTraces(:,:,1)');
%     plot(1:total,ipspMinima);
%     hold on;
% end

%Increasing inhibitory input to dendritic shaft on same branch as dendritic
%spine
% Ii_setup = cell(1,6);
% Ii_setup{1,1} = zeros(n,45); Ii_setup{1,1}(150:200,45)=1; Ii_setup{1,1}(150:200,17)=1;
% Ii_setup{1,2} = zeros(n,45); Ii_setup{1,2}(150:200,45)=1; Ii_setup{1,2}(150:200,15)=1;
% Ii_setup{1,3} = zeros(n,45); Ii_setup{1,3}(150:200,45)=1; Ii_setup{1,3}(150:200,11)=1;
% Ii_setup{1,4} = zeros(n,45); Ii_setup{1,4}(150:200,45)=1; Ii_setup{1,4}(150:200,3)=1;
% Ii_setup{1,5} = zeros(n,45); Ii_setup{1,5}(150:200,45)=1; Ii_setup{1,5}(150:200,2)=1;
% Ii_setup{1,6} = zeros(n,45); Ii_setup{1,6}(150:200,45)=1; Ii_setup{1,6}(150:200,1)=1;
% 
% figure
% for experiment = 1:6
%     total = 50;
%     voltageTraces = zeros(50,n,45);
%     for x = 1:50
%         Ie = zeros(n,45);
%         Ii = Ii_setup{1,experiment}*x;
%         
%         [v S epsp] = runmodel(Ie,Ii,n);
%         voltageTraces(x,:,:) = v;
%         [experiment x]
%     end
%     
%     ipspMinima = min(voltageTraces(:,:,1)');
%     plot(1:total,ipspMinima);
%     hold on;
% end

%Increasing inhibitory input to dendritic shaft on same branch as
%excitatory input onto spine
% Ii_setup = cell(1,6);
% Ii_setup{1,1} = zeros(n,45); Ii_setup{1,1}(150:200,17)=1;
% Ii_setup{1,2} = zeros(n,45); Ii_setup{1,2}(150:200,15)=1;
% Ii_setup{1,3} = zeros(n,45); Ii_setup{1,3}(150:200,11)=1;
% Ii_setup{1,4} = zeros(n,45); Ii_setup{1,4}(150:200,3)=1;
% Ii_setup{1,5} = zeros(n,45); Ii_setup{1,5}(150:200,2)=1;
% Ii_setup{1,6} = zeros(n,45); Ii_setup{1,6}(150:200,1)=1;
% 
% voltagesToPlot = cell(1,10);
% intensitiesOfInterest = [5,10,15,20,25,30,35,40,45,50];
% 
% for experiment = 1:6
%     total = 50;
%     voltageTraces = zeros(50,n,45);
%     for x = 1:50
%         Ie = zeros(n,45);
%         Ie(150:200,45) = x;
%         Ii = Ii_setup{1,experiment}*x;
%         
%         [v S epsp] = runmodel(Ie,Ii,n);
%         voltageTraces(x,:,:) = v;
%         [experiment x]
%         
%         if sum(ismember(intensitiesOfInterest,x))~=0
%             properIndex = find(intensitiesOfInterest==x);
%             voltagesToPlot{properIndex} = v(:,1);
%         end
%     end
%     
%     figure(1)
%     ipspMinima = min(voltageTraces(:,:,1)');
%     plot(1:total,ipspMinima);
%     hold on;
%     title('Minimum');
%     
%     figure(2)
%     epspMaxima = max(voltageTraces(:,:,1)');
%     plot(1:total,epspMaxima);
%     hold on;
%     title('Maximum');
%     
%     figure(3)
%     plot(1:total, epspMaxima-ipspMinima); 
%     hold on;
%     title('Amplitude');
%     
%     figure(4)
%     for int = 1:10
%         subplot(2,5,int)
%         hold on;
%         plot(voltagesToPlot{int});
%         titleString = ['Intensity ' num2str(intensitiesOfInterest(int))];
%         title(titleString);
%     end
%     hold on;
% 
%     figure(5)
%     plot(voltagesToPlot{1})
%     hold on;
%     title('Intensity 5');
%     
% end

%Increasing inhibitory input to dendritic shaft at various times relative
%to the excitatory input
% Ii_setup = cell(1,5);
% Ii_setup{1,1} = zeros(n,45); Ii_setup{1,1}(100:150,17)=1;
% Ii_setup{1,2} = zeros(n,45); Ii_setup{1,2}(125:175,17)=1;
% Ii_setup{1,3} = zeros(n,45); Ii_setup{1,3}(150:200,17)=1;
% Ii_setup{1,4} = zeros(n,45); Ii_setup{1,4}(175:225,17)=1;
% Ii_setup{1,5} = zeros(n,45); Ii_setup{1,5}(200:250,17)=1;
% 
% voltagesToPlot = cell(1,10);
% intensitiesOfInterest = [5,10,15,20,25,30,35,40,45,50];
% 
% for experiment = 1:5
%     total = 50;
%     voltageTraces = zeros(50,n,45);
%     for x = 1:50
%         Ie = zeros(n,45);
%         Ie(150:200,45) = x;
%         Ii = Ii_setup{1,experiment}*x;
%         
%         [v S epsp] = runmodel(Ie,Ii,n);
%         voltageTraces(x,:,:) = v;
%         [experiment x]
%         
%         if sum(ismember(intensitiesOfInterest,x))~=0
%             properIndex = find(intensitiesOfInterest==x);
%             voltagesToPlot{properIndex} = v(:,1);
%         end
%     end
%     
%     figure(1)
%     ipspMinima = min(voltageTraces(:,:,1)');
%     plot(1:total,ipspMinima);
%     hold on;
%     title('Minimum');
%     
%     figure(2)
%     epspMaxima = max(voltageTraces(:,:,1)');
%     plot(1:total,epspMaxima);
%     hold on;
%     title('Maximum');
%     
%     figure(3)
%     plot(1:total, epspMaxima-ipspMinima); 
%     hold on;
%     title('Amplitude');
%     
%     figure(4)
%     for int = 1:10
%         subplot(2,5,int)
%         hold on;
%         plot(voltagesToPlot{int});
%         titleString = ['Intensity ' num2str(intensitiesOfInterest(int))];
%         title(titleString);
%     end
%     hold on;
%     
%     figure(5)
%     plot(voltagesToPlot{1})
%     hold on;
%     title('Intensity 5');
%     legend('Exp 1','Exp 2','Exp 3','Exp 4','Exp 5');
%     
% end
