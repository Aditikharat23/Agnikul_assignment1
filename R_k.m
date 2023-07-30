classdef R_k
    methods(Static)
         function dstatedt = dmp(t,state)
             m   = 0.3;
             b  = 0.77;
             k=0.11;
            x= state(1);
            x_dot  = state(2);
            x_ddot= (-b*x_dot - k*x)/m;
            dstatedt = [x_dot,x_ddot];
         end
         function out= sol_dmp(t)
             m   = 0.3;
             b  = 0.77;
             k=0.11;
             root  = roots([m b k]);
             C1 = -2735301378778141/2038317670839040950;
             C2 = 4350165479555896/203831767083904095;
             

             out =( C1*exp(root(1)*t))+ (C2*exp(root(2)*t));

         end
         function [t_track,state_track] = conststep(state,h,f,tend)
            
             t_track(1)=0;
            state_track(1,:)= state;
            t = 0;
            iterations = floor(tend/h) +1;
            for i  = 1:iterations

                k1 = h*f(t,state);
                k2 = h*f(t+0.5*h,state + 0.5*k1);
                k3 = h*f(t+0.5*h,state + 0.5*k2);
                k4 = h*f(t+h,state+ k3);


                state = state + (k1+2*k2+2*k3+k4)/6 ;

                state_track(i+1,:) = state;
                t = t+h;
                t_track(end+1) = t;
            end
        end
function [t_track,state_track,h_track] = stephalving45(state,tol,f,tend,h_ini)
            t_track(1)=0;
            state_track(1,:)= state;
            h_track(1) = h_ini;
            h  = h_ini;
            t = 0;i = 1;
            while(t<=tend)
                state1 = R_k.order4(t,state,h,@R_k.dmp);
                state2  =R_k.order4(t,state,0.5*h,@R_k.dmp);
                state3 = R_k.order4(t+0.5*h,state2,0.5*h,@R_k.dmp);
                delta = abs(state3-state1);

                if (delta(1) >tol) %% Error is large, h must decrease
                    h = h/2;
                else
                    
                    t = t+h;
                    t_track(end+1) = t;
                    h_track(end+1)=h;
                    i=i+1;
                    state = state1;
                    state_track(i,:) = state;
                    if delta(1)<tol/100 %%If the error is considerably smaller
                        h = 2*h;
                        if h>2
                            h = 2;
                        end
                        
                    end
                end
                end

end
function out_state = order4(t,state,h,f)
             k1 = h*f(t,state);
                k2 = h*f(t+0.5*h,state + 0.5*k1);
                k3 = h*f(t+0.5*h,state + 0.5*k2);
                k4 = h*f(t+h,state+ k3);


                out_state = state + (k1+2*k2+2*k3+k4)/6 ;
        end
    end
end
