close all 
clear 
L=25;%cm Length 
vi=0.188;% interstitial velocity cm/s
vs=0.1326;
eps=0.704;
vs=vi*eps;%superficial velocity 
cf=15;%g/ml
Nz=400;%#of GP
delta_z=L/Nz;
zbd=linspace(0,25,Nz+1);
zs=zbd(1:Nz)+delta_z/2;
tspan=[0 400];
H=3.49;%Hengry constant 
cn0=zeros(1,2*Nz);

v=vi;

k=18.3;% mass transfer constant 
Da=1.31*10^(-3);%cm^2/s


[t,cn] =ode45(@(t,cn) odefcn(t,cn,Nz,delta_z,v,cf,eps,k,Da,H), tspan,cn0);
ns=[];
cs=[];


for i = 1 : 2*Nz
    if rem(i,2)==1
        cs=[cs,cn(:,i)];
    else    
        ns=[ns,cn(:,i)];
    end 
end 

figure (1)

plot(t,cs(:,end)./cf)
axis([275 375 0 1])
ylabel('c/c_f')
xlabel('time(s)')
% figure (2)
% 
% plot(t,ns(:,end))
% 
% figure (3)
% 
% plot(zs,cs(end,:))
% 
% figure (4)
% 
% plot(zs,ns(end,:))



function dcndt=odefcn(t,cn,Nz,delta_z,v,cf,eps,k,Da,H)
    dcndt=zeros(2*Nz,1);
    for i = 1 : 2*Nz
        if rem(i,2)==1% odd eqn(25)
            if i == 1%boundary condition
                ci=cn(i);
                ciplus1=cn(i+2);
                ciminus1=cf;
                ni=cn(i+1);

                Fjminushalf=v*ciminus1-v*(ci-cf);
                Fjplushalf=v*ci-Da*(ciplus1-ci)/delta_z;
                dcndt(i)=(1/delta_z)*(Fjminushalf-Fjplushalf)-((1-eps)/eps)*(k*(ci-ni/H));
            elseif i == 2*Nz-1
                ci=cn(i);
                ciminus1=cn(i-2);
                ni=cn(i+1);

                Fjminushalf=v*ciminus1-Da*(ci-ciminus1)/delta_z;
                Fjplushalf=v*ci;
                dcndt(i)=(1/delta_z)*(Fjminushalf-Fjplushalf)-((1-eps)/eps)*(k*(ci-ni/H));
             
            else
                ci=cn(i);
                ciplus1=cn(i+2);
                ciminus1=cn(i-2);
                ni=cn(i+1);

                Fjminushalf=v*ciminus1-Da*(ci-ciminus1)/delta_z;
                Fjplushalf=v*ci-Da*(ciplus1-ci)/delta_z;
                dcndt(i)=(1/delta_z)*(Fjminushalf-Fjplushalf)-((1-eps)/eps)*(k*(ci-ni/H));
            end
        else %even LDF model 
            ci=cn(i-1);
            ni=cn(i);
            dcndt(i)=k*(ci-ni/H);
        end 
    end 
end

