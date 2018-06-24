%%%%implementation of NIST Levenberg-Marquardt Algorithm for fittingNon-Linear Least Square%%%%%%%%
%CYLINDER FITTING
%INPUT: p0
%OUTPUT: p
%Where, p and p0 are center x,y,z, a1,a2,a3 and radius r. a1,a2 and a3 are the normal direction of the cylinder axis
%by: Wahyudin P. Syam
%email: gebe_top@yahoo.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function p=NIST_LevMar_Cylinder2(p0,M) %with normalization, but a is subtituted by A elements
[n m]=size(M);
lamda=0.0001;
firstStep=1;
counter1=1;
stop1=0;
while (~stop1)
    p=p0;    
    lamda=lamda-0.04; 
    %NORMALIZATION --> for >2th iteration    
    if(firstStep>1) %implement normalization
%         %normalization x <-- point on axis xlosest to the origin        
%         X0=origin; %column vector
%         X1=[p0(1) p0(2) p0(3)]'; %column vector
%         t=1;
%         A=[p0(4) p0(5) p0(6)]'; %column vector: direction vector
%         X2=X1+t*A; %column vector and set t=1;
%         t=-((X1-X0)'*(X2-X1))./sqrt(sum((X2-X1).^2));
%         Xclosest=X1+t*A;
%         p0(1)=Xclosest(1);
%         p0(2)=Xclosest(2);
%         p0(3)=Xclosest(3);
        %normalization A=A/|A|
%         lengthA=sqrt(p0(4)^2+p0(5)^2+p0(6)^2);
%         p0(4)=p0(4)/lengthA; %A
%         p0(5)=p0(5)/lengthA; %B
%         p0(6)=p0(6)/lengthA; %C
        
    else
        firstStep=2;
    end
    
    %calculate vector a=A/|A|
    lengthA=sqrt(p0(4)^2+p0(5)^2+p0(6)^2);
    a(1)=p0(4);
    a(2)=p0(5);
    a(3)=p0(6);
    
    %Calculate distance vector d(i) and matrix F0
    for i=1:n
        %calculate gi=a.(xi-x) --> scalar
        g(i)=a(1)*(M(i,1)-p0(1))+a(2)*(M(i,2)-p0(2))+a(3)*(M(i,3)-p0(3));
        
        %calculate fi=sqrt(u^2+v^2+w^2) --> scalar
        uSquare=(a(3)*(M(i,2)-p0(2))-a(2)*(M(i,3)-p0(3)))^2;
        vSquare=(a(1)*(M(i,3)-p0(3))-a(3)*(M(i,1)-p0(1)))^2;
        wSquare=(a(2)*(M(i,1)-p0(1))-a(1)*(M(i,2)-p0(2)))^2;
        f(i)=sqrt(uSquare+vSquare+wSquare);
        
        %calcualte distance d
        d(i)=f(i)-p0(7);
   
        %calculate matrix F
        if(f(i)~=0)
            F0(i,1)=(a(1)*g(i)-(M(i,1)-p0(1)))/f(i);
            F0(i,2)=(a(2)*g(i)-(M(i,2)-p0(2)))/f(i);
            F0(i,3)=(a(3)*g(i)-(M(i,3)-p0(3)))/f(i);
            F0(i,4)=g(i)*(a(1)*g(i)-(M(i,1)-p0(1)))/f(i);;
            F0(i,5)=g(i)*(a(2)*g(i)-(M(i,2)-p0(2)))/f(i);
            F0(i,6)=g(i)*(a(3)*g(i)-(M(i,3)-p0(3)))/f(i);
        else %f(i)=0
            F0(i,1)=sqrt(1-a(1)^2);
            F0(i,2)=sqrt(1-a(2)^2);
            F0(i,3)=sqrt(1-a(3)^2);
            F0(i,4)=g(i)*sqrt(1-a(1)^2);
            F0(i,5)=g(i)*sqrt(1-a(2)^2);
            F0(i,6)=g(i)*sqrt(1-a(3)^2);
        end
        F0(i,7)=-1;
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
            %calculate gi=a.(xi-x) --> scalar
            g(i)=a(1)*(M(i,1)-pnew(1))+a(2)*(M(i,2)-pnew(2))+a(3)*(M(i,3)-pnew(3));
        
            %calculate fi=sqrt(u^2+v^2+w^2) --> scalar
            uSquare=(a(3)*(M(i,2)-pnew(2))-a(2)*(M(i,3)-pnew(3)))^2;
            vSquare=(a(1)*(M(i,3)-pnew(3))-a(3)*(M(i,1)-pnew(1)))^2;
            wSquare=(a(2)*(M(i,1)-pnew(1))-a(1)*(M(i,2)-pnew(2)))^2;
            f(i)=sqrt(uSquare+vSquare+wSquare);
        
            %calcualte distance d
            d(i)=f(i)-pnew(7);
        end
        Jnew=sum(d.^2); 
        if((pnew(1)==p0(1)) &&(pnew(2)==p0(2))&&(pnew(3)==p0(3))&&(pnew(4)==p0(4))&&(pnew(5)==p0(5)) &&(pnew(6)==p0(6) && (pnew(7)==p0(7))) ) %if converge Pnew=P0 --> stop
            %p0 from normalized pnew
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
        origin=[p0(1) p0(2) p0(3)]'; %column vector
        p0=pnew';  %put p0 original format back to a "row vector"
    end
    counter1=counter1+1;
    if(counter1>10000)
        stop1=1;
    end
    
end %END OF OUTER WHILE LOOP

%p=p0;
end %END OF FUNCTION
