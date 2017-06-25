E=[0.356 0.511 0.662]*10^6*1.602*10^-19; % in J
A0=[397000 374000 371000]; % in Bq
th=[10.54 2.603 30.17]; %in a
thInS=th*365*24*3600; % in s
t1=650494368;
t2=t1+12*3600;
m=75;
D=zeros(1,3);
ints=zeros(1,3);

funBa=@(t) exp(-t/thInS(1)*log(2));
funNa=@(t) exp(-t/thInS(2)*log(2));
funCs=@(t) exp(-t/thInS(3)*log(2));
ints(1)=integral(funBa,t1,t2);
ints(2)=integral(funNa,t1,t2);
ints(3)=integral(funCs,t1,t2);

for i=1:3
    D(i)=E(i)*A0(i)/m*ints(i);
end

sum(D)
