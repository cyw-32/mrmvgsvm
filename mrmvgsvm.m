function [Z,P] = svc2(A1,A2,B1,B2,det,epis1,epis2)

%初始化Z,P
[n1,d1]=size(A1);e=ones(n1,1);
MA1=[A1 e];
G1= [A1 e]'*[A1 e];
[~,d2]=size(A2);
MA2=[A2 e];
G2= [A2 e]'*[A2 e];

[n2,~]=size(B1);e=ones(n2,1);
MB1=[B1 e];
H1= [B1 e]'*[B1 e];
[~,d2]=size(B2);
MB2=[B2 e];
H2= [B2 e]'*[B2 e];

K1_1=blkdiag(G1,G2);
K1_2=fliplr(blkdiag(fliplr(MA1'*MA2),fliplr(MA2'*MA1)));
K1=(1+epis1)*K1_1+(-epis1)*K1_2;
K1=K1+det*eye(size(K1));

T1=blkdiag(H1,H2);
[eigVector,eigValue]=eig(K1,T1);%A是向量，B是特征值,第一列向量对应其最小特征值
eigValue=diag(eigValue);
[eigValue,index1]=min(eigValue);
Z=eigVector(:,index1(1,1)); %view 1&2 plane one
%--------------------------------------------------------------------------
K2_1=blkdiag(H1,H2);
K2_2=fliplr(blkdiag(fliplr(MB1'*MB2),fliplr(MB2'*MB1)));
K2=(1+epis1)*K2_1+(-epis1)*K2_2;
K2=K2+det*eye(size(K2));

T2=blkdiag(G1,G2);  
[eigVector2,eigValue2]=eig(K2,T2);
eigValue2=diag(eigValue2);
[eigValue2,index2]=min(eigValue2);
P=eigVector2(:,index2(1,1));
%%%%%%%%%%%%%%%%%%%%%%%%initilzation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_t=50;iter=1;iter2=1;
objList11=[];
value1=inf;value2=inf;
Q=ones(n1,n1);Q2=ones(n2,n2);
d11=d1+1;d21=d2+1;
E1=diag(ones(d11,1));
E2=diag(ones(d21,1));
I=blkdiag(E1,E2);
v=sqrt(sum(Q.*Q,2));
D=diag(1./(v));
v2=sqrt(sum(Q2.*Q2,2));
D2=diag(1./(v2));

while ( iter <= max_t && value1>0.00001 )%
w1=Z(1:size(A1,2),1);%一视角   
bias11=Z(size(A1,2)+1,1);
w2=Z(size(A1,2)+2:end-1,1);
bias21=Z(end,1);

%更新Q
Q11=2*epis1*(A2*w2)*(A2*w2)'+epis2*D;%%
Q12=2*epis1*(A2*w2)*(A1*w1)';
Q  = Q11\Q12;
v = sqrt(sum(Q.*Q,2));
D = diag(1./(v));
%更新Z
MQ1=A1*E1(1:d1,:);MQ2=Q'*A2*E2(1:d2,:);
K1_2=[MQ1,-MQ2]'*[MQ1,-MQ2];
K1=K1_1+epis1*K1_2+det*I;

[eigVector,eigValue]=eig(K1,T1);%A是向量，B是特征值,第一列向量对应其最小特征值
eigValue=diag(eigValue);
[eigValue,index1]=min(eigValue);
Z=eigVector(:,index1(1,1)); 
ww1=Z(1:size(MA1,2));
ww2=Z(d1+2:end);
obj(iter) =(Z'*K1*Z+epis2*sum(v))/(Z'*T1*Z);
    if iter >1
        value1=abs(obj(iter)-obj(iter-1)); 
        objList11=[objList11;value1];
    end
    iter=iter+1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for iter2 = 1:max_t
while (iter2 <= max_t && value2>0.00001)
u1=P(1:size(B1,2),1);%二视角
bias12=P(size(B1,2)+1,1);
u2=P(size(B1,2)+2:end-1,1);
bias22=P(end,1);

%更新Q2
Q21=2*epis1*(B2*u2)*(B2*u2)'+epis2*D2;
Q22=2*epis1*(B2*u2)*(B1*u1)';
Q2  = Q21\Q22;

v2 = sqrt(sum(Q2.*Q2,2));
D2 = diag(1./(v2));
%更新P
MQ21=B1*E1(1:d1,:);MQ22=Q2'*B2*E2(1:d2,:);
% K2_2=fliplr(blkdiag(fliplr(MQ21'*MQ22),fliplr(MQ2'*MQ22)));
% K2=(1+epis1)*K2_1+(-epis1)*K2_2;
% K2=K2+det*eye(size(K2));

K2_2=[MQ21,-MQ22]'*[MQ21,-MQ22];
K2=K2_1+epis1*K2_2+det*I;

[eigVector2,eigValue2]=eig(K2,T2);
eigValue2=diag(eigValue2);
[eigValue2,index2]=min(eigValue2);
P=eigVector2(:,index2(1,1));
uu1=P(1:size(MB1,2));
uu2=P(d1+2:end);
obj(iter2) =(P'*K2*P+epis2*sum(v2))/(P'*T2*P);
    if iter2 >1
       value2=abs(obj(iter2)-obj(iter2-1)); 
    end
    iter2=iter2+1;
end

