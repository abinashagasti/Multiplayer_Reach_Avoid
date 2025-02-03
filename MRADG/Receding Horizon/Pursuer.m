classdef Pursuer < handle
% Pursuer class
    properties
        index % index of pursuer agent
        position % 3x1 vector representing current position
        speed % scalar representing maximum speed
        velocity % 3x1 vector representing current velocity
        status % indicator denoting status of pursuer
        % 0 - idle pursuer
        % -1 - chasing pursuer
        % 1 - pursuer after neutralisation of assigned evader
        evader % index of the evader assigned to be pursued
        % 0 if no evader is assigned
    end

    methods

        function p = Pursuer(initPos, speed, index)
            % Constructor, assigns initial values
            p.index = index;
            p.position = initPos;
            p.speed = speed;
            p.status = 0;
            p.evader = 0;
        end

        function updatePos(p,position)
            sz_input = size(position);
            sz_output = size(p.position);
            if sz_input~=sz_output
                if sz_input==fliplr(sz_output)
                    p.position = position';
                else
                    error('Wrong position dimensions.')
                end
            else
                p.position = position;
            end
        end

        function velocity = return_velocity(p,evader,B,tolerance)
            alpha=evader.speed/p.speed;
            xp = p.position;
            xe = evader.position;
            xc = (xe-alpha^2*xp)/(1-alpha^2);
            Rc = norm(xc);
            rc = (alpha/(1-alpha^2))*norm(xp-xe);
            if B>=0
                if abs(norm(p.position-evader.position))>tolerance
                    gradV = (alpha^2/(1-alpha^2))*(-xc/Rc+(1/(1-alpha^2))*((xe-xp)/rc));
                    % finding the gradient of the Value function
                    rhop=norm(gradV);
                    velocity = (p.speed/rhop)*gradV;
                else
                    velocity=zeros(3,1);
                end
            elseif B<0
                if abs(norm(xe))>1e-3
                    gradV=-xp/norm(xp);
                    rhop=norm(gradV); 
                    velocity = (p.speed/rhop)*gradV;
                else
                    velocity=zeros(3,1);
                end
            end
        end

    end

end