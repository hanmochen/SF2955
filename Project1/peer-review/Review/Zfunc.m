function [Zindex] = Zfunc(Zindex,P,prob)

    P = cumsum(P,2); % comsum by the columns
    
    for i = 2:length(prob)
        
        if (prob(i) < P(Zindex,1))            
            Zindex = 1;
            
        elseif (prob(i) < P(Zindex,2))
            Zindex = 2;  
            
        elseif (prob(i) < P(Zindex,3))
            Zindex = 3;  
            
        elseif (prob(i) < P(Zindex,4))
            Zindex = 4;  
            
        else
            Zindex = 5;  
        end
    end
end
