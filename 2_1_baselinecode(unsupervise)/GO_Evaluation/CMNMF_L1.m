function [L,W,H1,H2] = CMNMF_L1(K,MaxIter,Inner_MaxIter,V1,V2,W,H1,H2,M,alpha,gama,lamat1,lamat2)
 L=zeros(1+2*Inner_MaxIter,MaxIter);
 e=ones(K,K);
for i=1:MaxIter
   
    L(1,i)=objective(alpha,gama,lamat1,lamat2,V1,V2,W,H1,H2,M);
    s=2;
    for t=1:Inner_MaxIter
      W=W.*(V1*H1'+alpha*V2*H2')./(W*(H1*H1')+alpha*W*(H2*H2')+lamat1*W*e+eps);
      L(s,i)=objective(alpha,gama,lamat1,lamat2,V1,V2,W,H1,H2,M);

        if(L(s,i)>L((s-1),i))
          break;
        end
        s=s+1;
    end
    
    for t=1:Inner_MaxIter
      H1=H1.*(W'*V1+0.5*gama*H2*M')./(W'*W*H1+lamat2*e*H1+eps);
      H2=H2.*(alpha*W'*V2+0.5*gama*H1*M)./(alpha*(W'*W*H2)+lamat2*e*H2+eps);
      L(s,i)=objective(alpha,gama,lamat1,lamat2,V1,V2,W,H1,H2,M);
        if(L(s,i)>L((s-1),i))
            break;
        end
       s=s+1;
    end

end
end

function [O]= objective(alpha,gama,lamat1,lamat2,V1,V2,W,H1,H2,M)
    
    O1=sum(sum((V1-W*H1).^2))+alpha*sum(sum((V2-W*H2).^2));
    O2=-gama*trace(H1*M*H2')+lamat1*sum(sum(W,2).^2)+lamat2*(sum(sum(H1).^2)+sum(sum(H2).^2));
    O=O1+O2;
    
end
