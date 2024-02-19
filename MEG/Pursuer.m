classdef Pursuer < handle
    properties (SetAccess = public)
        position % 2x1 vector representing current position
        speed % scalar representing maximum speed
    end

    methods

        function p = Pursuer(initPos, speed)
            % Constructor, assigns initial values
            p.position = initPos;
            p.speed = speed;
        end
        
        function updatePos(p, position)
            % Updates position and velocity based on evader position
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

        function position = getPos(p)
            position = p.position;
        end

        function theta = heading_direction(p,evader_positions)
            theta=0;
        end

    end
end
