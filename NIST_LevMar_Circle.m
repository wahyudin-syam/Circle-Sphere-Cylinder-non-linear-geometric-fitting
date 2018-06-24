%%%%implementation of NIST Levenberg-Marquardt Algorithm for fittingNon-Linear Least Square%%%%%%%%
%CIRCLE FITTING
%INPUT: p0
%OUTPUT: p
%where: p and p0 are center x,y,and radius r
%by: Wahyudin P. Syam
%email: gebe_top@yahoo.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function p=NIST_LevMar_Circle(p0,M)
[n m]=size(M);
lamda=0.0001;
counter1=1;
stop1=0;
while (~stop1)
    p=p0;    
    lamda=lamda-0.04; 
    %NO NORMALIZATION FOR SPHERE & 2D-CIRCLE
    %Build distance d vector
    for i=1:n
        d(i)=sqrt((M(i,1)-p0(1))^2+(M(i,2)-p0(2))^2)-p0(3);
    end
    %Build F0 matrix
    for i=1:n
        F0(i,1)=-(M(i,1)-p0(1))/sqrt((M(i,1)-p0(1))^2+(M(i,2)-p0(2))^2);
        F0(i,2)=-(M(i,2)-p0(2))/sqrt((M(i,1)-p0(1))^2+(M(i,2)-p0(2))^2);        
        F0(i,3)=-1;
    end    
    size(d);
    size(F0);
    %calculate matrix U mxm = F0'F0 where m=4 (x,y,z,r)
    U=F0'*F0;
    %calculate vector v mx1=F0'*d where m=4 (x,y,z,r)
    v=F0'*d';
    %calculate J0=sum (d^2)
    j0=sum(d.^2);    
    
    stop2=0;
    counter2=1;
    while (~stop2)
        lamda=lamda+10;
        [n1 m1]=size(U);
        H=U+lamda*(eye(m1)+diag(diag(U)));
        %solving Hx=-v by using Cholesky decomposition
        R=chol(H);
        x=inv(R'*R)*-v;
        %x=inv(H)*-v;
        %%%size(x)
        %%%size(p0')
        pnew=p0'+x;
        for i=1:n
            d(i)=sqrt((M(i,1)-pnew(1))^2+(M(i,2)-pnew(2))^2)-pnew(3);
        end
        Jnew=sum(d.^2); 
        if((pnew(1)==p0(1)) &&(pnew(2)==p0(2))&&(pnew(3)==p0(3)) ) %if converge Pnew=P0 --> stop            
            p0=pnew; 
            p=p0;
            return;
        end
        
        if(Jnew<j0)
            stop2=1;
        end
        counter2=counter2+1;
        if(counter2>10000)
            stop2=1;
        end
    end %END OF INNER WHILE LOOP
    
    if(Jnew<j0)
        p0=pnew';  %p0 original format is a row vector      
    end
    counter1=counter1+1;
    if(counter1>1000)
        stop1=1;
    end
end %END OF OUTER WHILE LOOP

%p=p0;
end %END OF FUNCTION
