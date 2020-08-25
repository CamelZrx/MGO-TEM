function [z1]=normalize_up_to_t(eo)
% EAP is the set of positve local minimums
% local minimums, eo is the new local minimum, MP is the mean vector.
global EAP MP CMP Nt Np fmin fmax

Nt_new = fmax-2*min(0,fmin)+1;                       % the new normalization factor

if(Nt_new>=Nt)                                   % new normalization is found is found
    eap = EAP*Np/Nt_new;                    
    MP = MP*Nt/Nt_new; 
    Nt = Nt_new;
    m = mean([eap (eo-2*min(0,fmin)+1)/Nt]);       % new entry in the mean vector
    MP = [MP m];    
else
    if(Nt~=0)
        eap = EAP*Np/Nt;                         % normalization up to t
        m = mean([eap (eo-2*min(0,fmin)+1)/Nt]);   % new entry in the mean vector
        MP = [MP m];
    end
end

if (length(MP)>1)
    CMP = [CMP mean(MP)];
end

if(isempty(EAP))
    if(EAP>eo)
        Np = EAP;                               % initialize the Np
        EAP = 1;            
    else
        Np = fmax-2*min(0,fmin)+1;
        EAP = EAP/Np;                           % initialize the EAP
    end
end

if(Nt_new>=Nt)
    EAP = EAP*Np/Nt_new;                       % normalization up to t
    Np = Nt_new;    
else      
end

z1 = Np;
end