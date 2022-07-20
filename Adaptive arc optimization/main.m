%%%%%Path=[x1 y1;
%%%%%%%%%% x2 y2
%%%%%%%%%% ....];



Pathe=Path;
    Num=size(Path,1)+1;
    NumM=size(MAX,1);
    we=1;
    A1=Pathe(we,:)+.5;
    
    while we<Num-2
        A2=Pathe(we+1,:)+.5;
        A3=Pathe(we+2,:)+.5;
        A12=A1-A2;
        A23=A2-A3;
        QA=Pathe(we,:)-A2;
        
        if QA(1)*A23(2)-QA(2)*A23(1)~=0  
            d1=sqrt((A1(1)-A2(1))^2+(A1(2)-A2(2))^2);
            d2=sqrt((A3(1)-A2(1))^2+(A3(2)-A2(2))^2);
            if d1<=d2
                p=A1;
                q=A3;
            else
                p=A3;
                q=A1;
            end                            
            io=1;         
            A2q=q-A2;
            A2p=p-A2;
            p=0.5*A2p+A2;
            
            while io==1
                A2q=q-A2;
                A2p=p-A2;
                d3=sqrt(A2q(1)^2+A2q(2)^2);
                d4=sqrt(A2p(1)^2+A2p(2)^2);
                A2qq=A2q/d3*d4;
                A2he=A2qq+A2p;
                qq=A2qq+A2;      
                syms X0 Y0          
                eq1=(X0-p(1))*(p(1)-A2(1))+(Y0-p(2))*(p(2)-A2(2));
                eq2=(X0-qq(1))*(qq(1)-A2(1))+(Y0-qq(2))*(qq(2)-A2(2));
                [x0,y0]=solve(eq1,eq2,X0,Y0);
                x0=eval(x0);y0=eval(y0);       
                R=sqrt((x0-p(1))^2+(y0-p(2))^2); 
                xx=(min(p(1),qq(1)):0.001:max(p(1),qq(1)));
                if A2(2)>y0    
                    yy=sqrt(R^2-(xx-x0).^2)+y0;
                else
                    yy=-sqrt(R^2-(xx-x0).^2)+y0;
                end
                io=0;
                %
                
                for m=1:NumM  
                    for n=1:NumM
                        if  MAX(m,n)==1
                            bb=[m,n;m+1,n;m+1,n+1;m,n+1];
                            if ((min(xx)<=m&&max(xx)>=m)||(min(xx)<=m+1&&max(xx)>=m+1))&&((min(yy)<=n&&max(yy)>=n)||(min(yy)<=n+1&&max(yy)>=n+1))
                                if (bb(1,1)-x0)^2+(bb(1,2)-y0)^2<=R^2&&(bb(2,1)-x0)^2+(bb(2,2)-y0)^2<=R^2&&...
                                        (bb(3,1)-x0)^2+(bb(3,2)-y0)^2<=R^2&&(bb(4,1)-x0)^2+(bb(4,2)-y0)^2<=R^2
                                elseif  (bb(1,1)-x0)^2+(bb(1,2)-y0)^2>=R^2&&(bb(2,1)-x0)^2+(bb(2,2)-y0)^2>=R^2&&...
                                        (bb(3,1)-x0)^2+(bb(3,2)-y0)^2>=R^2&&(bb(4,1)-x0)^2+(bb(4,2)-y0)^2>=R^2
                                else
                                    io=1;
                                    break
                                end
                            end
                        end
                        if  io==1
                            break
                        end
                    end
                    if  io==1
                        break
                    end
                end
                if io==1
                    p=0.05*A2p+A2;
                end
            end
            AA=[xx(end),yy(end)];
            AA2=AA-A2;
            AA3=A3-A2;
             ec{we}=[AA;A1;A2;A3];
             
            if AA(1)<=max(A2(1),A3(1))+10e-3&&AA(1)>=min(A2(1),A3(1))-10e-3&&AA(2)<=max(A2(2),A3(2))+10e-3&&AA(2)>=min(A2(2),A3(2))-10e-3
                A1=AA;
                Pathee{we}=[xx',yy'];
            else
                A1=[xx(1),yy(1)];
                Pathee{we}=flipud([xx',yy']);
            end
            ecc{we}=[A1];
        else
            A1=A2;
        end
        we=we+1;
    end
    for i=1:size(Pathee,2)
        PLO=Pathee{i};
           plot(PLO(:,1)',PLO(:,2)','Color',[255/255 0/255 0/255],'linewidth',2.5) 
        hold on
    end
    