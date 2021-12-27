%Mass and Energy Balance II Project (Mateo/Nick/Ricky)
%All rates are per second
clear;clc
p=1;o=0;i=1; %loop conditions
n2=0; %Guess Recycle mole rate
T3=250; %Temperature of Reactor (250-400)

%Allocating memory for plots 
Bplot=zeros(1,151); %Butanol
Rplot=zeros(1,151); %Recycle rates
Yoplot=zeros(1,151); %Yield 
QhexPlot=zeros(1,151); %Heat exchanger Q
QsepPlot=zeros(1,151); %Seperator Q
QrctPlot=zeros(1,151); %Reactor Q

%Heat of formations (kJ/mol)
Hf_W=-241.83; %Water
Hf_B=-275.1; %Butanol
Hf_E=-235.31; %Ethanol
Hf_H=-316.2; %Hexanol

%Heat of Reactions
H_r1= Hf_B + Hf_W - 2*Hf_E; %Reaction 1
H_r2= Hf_H + Hf_W - (Hf_E+Hf_B); %Reaction 2
while p>o
%Heat Capacities (kJ/mol)
Cp_W=@(T)33.46e-3+T*.688e-5+T.^2*0.7604e-8; %Water
Cp_E=@(T)61.34e-3+T*15.72e-5-T.^2*8.749e-8; %Ethanol
Cp_B=@(T)(0.2702*T+29)/1000; %Butanol
Cp_H=@(T)0.1636*log(T)-0.7815; %Hexanol
 
    
%Feed Variables   
n1 = 100; %Fresh Mole feed
E1= n1*0.96; %Ethanol in Feed
W1= n1*0.04; %Water in Feed
T1= 30; %Feed Temperature

%Recycle Guess Variables
E2=0.96*n2; %Ethanol in recycle Feed
W2=0.04*n2; %Water in recycle Feed

%After Heat Exchanger
n3=n1+n2; %Combined molar feeds
E3= 0.96*n3; %Ethanol after HEX
W3=0.04*n3; %Water after HEX
Xcon=-2.294+T3*1.438e-2-T3^2*1.592e-5; %Single Pass Conversion
Ysp= 2.201 -T3*6.203e-3 + T3^2*3.610e-6; %Single Pass yield

%Matrix to solve systems of Eqn inv(A)*B= X where X=[E4;B4;W4;H4;eR1;eR2]
A=[1 0 0 0 2 1;0 1 0 0 -1 1;0 0 1 0 -1 -1;0 0 0 1 0 -1;1 0 0 0 0 0;...
    Ysp 2 0 0 0 0];
B=[E3;0;W3;0;(-E3)*(Xcon)+E3;(Ysp)*E3];
Solution=(A)\B;

%After Reactor
E4=Solution(1); %Ethanol
B4=Solution(2); %n-butanol 
W4=Solution(3); %Water 
H4=Solution(4); %Hexanol
eR1=Solution(5); %Extent of reaction 1
eR2=Solution(6); %Extent of reaction 2

%Recycle Stream 
E5=E4; %All leftover ethanol
n5=E5/(0.96); %Reycle Stream rate
W5=n5*(0.04); %Water needed to be 4% ratio

%After Seperator
B6=B4; %All n-butanol leaves
H6=H4; %All hexanol leaves
W6=(W4-W5); %Waste water

if abs(E5-E2) <.1 %Only runs if recycle stream is correct
    
    H_rct=(E4*integral(Cp_E,25,T3)+B4*integral(Cp_B,25,T3)...
        +W4*integral(Cp_W,25,T3)+H4*integral(Cp_H,25,T3)...
        )-(E3*integral(Cp_E,25,T3)+W3*integral(Cp_W,25,T3)); %Calculating delta H in reactor
    
    Qrct=H_rct+eR1*H_r1+eR2*H_r2; %Q for Reactor
    
    Yield=200*(B6/E1); %Overall yield at each temperature
    
    Qhex=(E3*integral(Cp_E,25,T3)+W3*integral(Cp_W,25,T3))-...
        (E2*integral(Cp_E,25,90)+W2*integral(Cp_W,25,90)+...
        E1*integral(Cp_E,25,T1)+W1*integral(Cp_W,25,T1)); %Heat duty in HEX
    
    Qsep=(B6*integral(Cp_B,25,T1)+W6*integral(Cp_W,25,T1)+...
        H6*integral(Cp_H,25,T1)+W5*integral(Cp_W,25,90)...
        +E5*integral(Cp_E,25,90))-(E4*integral(Cp_E,25,T3)+...
        B4*integral(Cp_B,25,T3)+W4*integral(Cp_W,25,T3)...
        +H4*integral(Cp_H,25,T3)); %Heat duty in Seperator
    
    
    %Storing Plot points
    Yoplot(i)=Yield; %Yield
    Bplot(i)=B6; %Butanol Production
    Rplot(i)=n2; %Recycle Rates
    QhexPlot(i)=Qhex; %Q for HEX
    QsepPlot(i)=Qsep; %Q for SEP
    QrctPlot(i)=Qrct; % for Reactor    
    
    T3=T3+1; %Increases Reactor temp by 1 per loop
    i=i+1; %Sets values through array
    
    if T3==274 %Finding Values at optimal Reactor Temperature
        Opt_Qhex=Qhex;
        
        Opt_Qrct=Qrct;
        
        Opt_Qsep=Qsep;
    
        Opt_B=B6;
        Opt_yield=200*(B6/E1);
    end
    
    if T3==401 %Stops loop at 400 C
        break
    end
    n2=0;
end
n2=n2+.1;
end

%Graphs
Tplot=linspace(250,400,151);

subplot(3,2,1);
plot(Tplot,Bplot,'y','LineWidth',2);xlabel('Temperature (°C)');...
    ylabel('Butanol Production (mol/s)');title('Butanol Production vs Temp')
subplot(3,2,2)
plot(Tplot,Rplot,'m','LineWidth',2);xlabel('Temperature (°C)');...
    ylabel('Recycle Stream (mol/s)');title('Recycle Rate vs Temp')
subplot(3,2,3)
plot(Tplot,Yoplot,'c','LineWidth',2);xlabel('Temperature (°C)');...
    ylabel('Overall Yield (%)');title('Overall Yield vs Temp')
subplot(3,2,4)
plot(Tplot,QhexPlot,'r','LineWidth',2);xlabel('Temperature (°C)');...
    ylabel('Heat Duty of HEX (kJ/s)');title('Heat Duty HEX vs Temp')
subplot(3,2,5)
plot(Tplot,QsepPlot,'g','LineWidth',2);xlabel('Temperature (°C)');...
    ylabel('Heat Duty of SEP (kJ/s)');title('Heat Duty SEP vs Temp')
subplot(3,2,6)
plot(Tplot,QrctPlot,'b','LineWidth',2);xlabel('Temperature (°C)');...
    ylabel('Heat Duty of Reactor (kJ/s)');title('Heat Duty RCT vs Temp')


%Optimization of reactor
fn=@(T)-(sqrt(-2.294+T*1.438e-2-T.^2*1.592e-5)+2.201 -T*6.203e-3 +...
    T.^2*3.610e-6); %Optimal reactor function as negative to find max
bestTemp=fminbnd(fn,250,400);

fprintf('The optimal reactor temperature is %.2f° C \n',bestTemp)
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
fprintf('Butanol Production: %.2f moles|| Yield: %.2f%%\n',Opt_B,Opt_yield)
fprintf('Seperator Q: %.2f kJ ||Heat Exchanger Q: %.2f kJ\n',Opt_Qsep,Opt_Qhex)
fprintf('Reactor Q: %.2f kJ \n',Opt_Qrct)
