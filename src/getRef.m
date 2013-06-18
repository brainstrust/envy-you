function [out] = getRef(t,name)


 %% OSTBY Input   
        global lengthpulse startpulse
        t0 = startpulse;    %10;
        lengtht1 = 10;
        t1 = t0+lengtht1;    %20;
        t2 = t0+lengthpulse;    %30;
        t3 = t1+lengthpulse;    %40;


 if strcmp('ft',name)
        
        
     
        % Compute gamma distribution as done in Ostby:
        alpha =  2;
        beta  =  5;
        f_input = 2.5;
        
        deltaT = t1-t0;
        gab= factorial(alpha + beta -1);
        ga = factorial(alpha -1);
        gb = factorial(beta -1);

            if (t>=t0 && t<=t1)
                out = f_input*gab/(ga*gb) .* (1-(t-t0)./deltaT).^(beta-1).*((t-t0)./deltaT).^(alpha-1); 
            elseif (t>t2 && t<t3)
                out = -f_input; 
            else
                out = 0;
            end
            
 elseif strcmp('fluxft',name)
     %out = createPulse(t,10,40,0,1,0.0005,0.0005);
     out = createPulse(t,t0,t1+lengthpulse,0,1,0.0005,0.0005);



%% Hannah's Input

 pulse_start = 100;
 pulse_end   = 200;
 
 elseif strcmp('rho',name)
     out = createPulse(t,pulse_start,pulse_end,0.1,0.7,1,1); %[-] fraction between zero and one
 elseif strcmp('J_K_s',name)
     out = createPulse(t,pulse_start,pulse_end,1,8,1); % [microMs-1]
 elseif strcmp('Glu',name)
     out = createPulse(t,100,250,0,2,3,3); 
 elseif strcmp('wss',name)
     out = createPulse(t,450,800,0,4,3,3); 
 elseif strcmp('fthannah',name)
     if t>=100 && t<160
        out =0.35*(1+tanh((t-101)/3));        
     else
        out =0; 
     end
 else
     fprintf('Error reference input unknown set to zero\n')
     out = 0;
 end
end