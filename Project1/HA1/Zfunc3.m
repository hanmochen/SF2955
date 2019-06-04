function [Zindex] = Zfunc3(Zindex,P,prob)

    P = cumsum(P,2); % comsum by the columns

    if (prob < P(Zindex,1))            
        Zindex = 1;

    elseif (prob < P(Zindex,2))
        Zindex = 2; 

    elseif (prob < P(Zindex,3))
        Zindex = 3; 
        
    elseif (prob < P(Zindex,4))
        Zindex = 4;

    else
        Zindex = 5; 
    end

end
