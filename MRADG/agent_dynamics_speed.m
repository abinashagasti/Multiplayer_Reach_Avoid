function dydt = agent_dynamics_speed(t,y)
% Agent dynamics for the reach-avoid game
    
    ve=y(7);vp=y(8);B=y(9);
    alpha=ve/vp;
    xp=y(4:6,1);xe=y(1:3,1);
    [xc,rc]=apollonius_parameters(xp,xe,alpha,0); % find the apollonius 
    % circle center and radius
    Rc=norm(xc);
    Rp=norm(xp);
    Re=norm(xe);
    
    if B>=0
        if abs(norm(xp-xe))>1e-2
            gradV=(1/(1-alpha^2))*[xc/Rc-(alpha^2/(1-alpha^2))*((xe-xp)/rc);
                alpha^2*(-xc/Rc+(1/(1-alpha^2))*((xe-xp)/rc))];%finding the 
            % gradient of the Value function
            rhoe=norm(gradV(1:3,1));
            rhop=norm(gradV(4:6,1));    
            dydt=[-(ve/rhoe)*gradV(1:3,1);(vp/rhop)*gradV(4:6,1);0;0;0];
%             dydt=[-(ve/Re)*xe;(vp/rhop)*gradV(4:6,1);0;0;0];
        else
            dydt=zeros(9,1);
        end
    elseif B<0
        if abs(norm(xe))>1e-3
            gradV=[-(1/alpha)*(xe/Re);xp/Rp];
            rhoe=norm(gradV(1:3,1));
            rhop=norm(gradV(4:6,1)); 
%             dydt=[(alpha/rhoe)*gradV(1:3,1);-(1/rhop)*gradV(4:6,1);0;0;0];
            dydt=[(ve/rhoe)*gradV(1:3,1);-(vp/rhop)*gradV(4:6,1);0;0;0];
        else
            dydt=zeros(9,1);
        end
    end

end