function [Zindex] = Zfunc2(Zindex,P,prob)
    N=1000;
    P = cumsum(P,2); % comsum by the columns
    
    for i = 2:N
        
        if (prob(i) < P(Zindex(i),1))            
            Zindex(i) = 1;
            
        elseif (prob(i) < P(Zindex(i),2))
            Zindex(i) = 2;  
            
        elseif (prob(i) < P(Zindex(i),3))
            Zindex(i) = 3;  
            
        elseif (prob(i) < P(Zindex(i),4))
            Zindex(i) = 4;  
            
        else
            Zindex(i) = 5;  
        end
    end
end
