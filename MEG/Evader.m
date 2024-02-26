classdef Evader < handle
    properties
        position % 2x1 vector representing current position
        speed % scalar representing maximum speed
    end
    methods
        function e = Evader(initPos, speed)
            % Constructor, assigns initial values
            e.position = initPos;
            e.speed = speed;
        end
        
        function updatePos(e, position)
            % Updates position based on random direction
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

        function position = getPos(e)
            position = e.position;
        end

        function velocity = heading_velocity(e, pursuer_position, target_position)
            xc = (e.position(1)+pursuer_position(1))/2;
            yc = (e.position(2)+pursuer_position(2))/2;
            m = (e.position(2)-pursuer_position(2))/(e.position(1)-pursuer_position(1));

            x_intercept = (m*(yc - target_position(2))+xc+m^2*target_position(1))/(1+m^2);
            y_intercept = target_position(2) + m*(x_intercept - target_position(1));
            
            velocity = [x_intercept - e.position(1), y_intercept - e.position(2)];
            velocity = (e.speed/norm(velocity,2))*(velocity);
        end

    end
end
