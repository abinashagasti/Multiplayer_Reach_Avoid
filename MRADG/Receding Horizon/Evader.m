classdef Evader < handle

    properties (SetAccess = public)
        index % index of evader agent
        position % 3x1 vector representing current position
        speed % scalar representing maximum speed
        velocity % 3x1 vector representing current velocity
        status % indicator denoting status of evader
        % 0 - free evader
        % 1 - captured evader
        % -1 - reached target
        pursuer % index of the pursuer assigned to pursuer
        % 0 if no pursuer is assigned
    end

    methods

        function e = Evader(initPos, speed, index)
            % Constructor, assigns initial values
            e.index = index;
            e.position = initPos;
            e.speed = speed;
            e.status = 0;
            e.pursuer = 0;
        end

        function updatePos(e,position)
            sz_input = size(position);
            sz_output = size(e.position);
            if sz_input~=sz_output
                if sz_input==fliplr(sz_output)
                    e.position = position';
                else
                    error('Wrong position dimensions.')
                end
            else
                e.position = position;
            end
        end

        function velocity = return_velocity(e,pursuer,B,tolerance)
            alpha=e.speed/pursuer.speed;
            xp = pursuer.position;
            xe = e.position;
            xc = (xe-alpha^2*xp)/(1-alpha^2);
            Rc = norm(xc);
            rc = (alpha/(1-alpha^2))*norm(xp-xe);
            if B>=0
                if abs(norm(e.position-pursuer.position))>tolerance
                    gradV = (1/(1-alpha^2))*(xc/Rc-(alpha^2/(1-alpha^2))*((xe-xp)/rc));
                    % finding the gradient of the Value function
                    rhoe=norm(gradV);
                    velocity = (-e.speed/rhoe)*gradV;
                else
                    velocity=zeros(3,1);
                end
            elseif B<0
                if abs(norm(xe))>1e-3
                    gradV=(1/alpha)*(xe/norm(xe));
                    rhoe=norm(gradV); 
                    velocity = (-e.speed/rhoe)*gradV;
                else
                    velocity=zeros(3,1);
                end
            end
        end

    end

end