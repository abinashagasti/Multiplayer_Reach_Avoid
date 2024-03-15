classdef Pursuer < handle
    properties (SetAccess = public)
        position % 2x1 vector representing current position
        speed % scalar representing maximum speed
        motor1 % motor variables 
        motor2
        motor3
        pos_vrpn
        ori_vrpn
        wheel_radius
        wheel_centre_radius
    end

    methods

        function p = Pursuer(initPos, speed)
            % Constructor, assigns initial values
            p.position = initPos;
            p.speed = speed;
            p.wheel_radius = 0.025;
            p.wheel_centre_radius = 0.15;
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

        function [psi_star,theta_star, m, xm, ym] = optimal_headings_Ei(p, evader_position, target_position)
            m = (evader_position(2)-p.position(2))/(evader_position(1)-p.position(1));
            xm = (evader_position(1)+p.position(1))/2;
            ym = (evader_position(2)+p.position(2))/2;
            xt = target_position(1); yt = target_position(2);
            xI = (m*(ym-yt)+xm+m^2*xt)/(1+m^2);
            yI = (m*(xm-xt)+yt+m^2*ym)/(1+m^2);

            psi_star = atan2(yI-evader_position(2),xI-evader_position(1));
            theta_star = atan2(yI-p.position(2),xI-p.position(1));
        end

        function d_i = objective_Ei(p, evader_position, target_position, r, theta)
            [psi_star, ~, m, xm, ym] = p.optimal_headings_Ei(evader_position, target_position);
            xt = target_position(1); yt = target_position(2);
            Tx = abs(m*yt+xt-m*ym-xm)/sqrt(1+m^2);
            Ty = abs(yt-m*xt-evader_position(2)+m*evader_position(1))/sqrt(1+m^2);
            k = norm(p.position-evader_position)/2;
            delta = atan2(p.position(2)-evader_position(2),p.position(1)-evader_position(1));
            d_i = Tx + (r/(2*k))*sqrt(Ty^2+k^2)*(-1-cos(theta+psi_star-2*delta));
        end

        function [theta_min, theta_max, min_evader, max_evader] = concave_domain(p, evader_positions, target_position)
            shape = size(evader_positions); % Shape is 2*n
            n = shape(2);
            theta_keys = 1:n;
            theta_values = zeros(1,n);
            for i=1:n
                [~,theta_star] = p.optimal_headings_Ei(evader_positions(:,i), target_position);
                theta_values(i) = theta_star;
            end
            [theta_values, sort_order] = sort(theta_values);
            % if theta_values(n)-theta_values(1)>pi
            %    error("The max range of heading angles is more than pi. Cannot handle such situations.")
            % end
            theta_keys = theta_keys(sort_order);
            [theta_largest, max_evader] = max(theta_values);
            [theta_smallest, min_evader] = min(theta_values);
            max_evader = theta_keys(max_evader);
            min_evader = theta_keys(min_evader);
            theta_min = max(theta_smallest, theta_largest-pi/2);
            theta_max = min(theta_largest, theta_smallest+pi/2);
        end

        function cost = objective_fun(p, evader_positions, target_position, r, theta)
            shape = size(evader_positions);
            n = shape(2);
            cost = 0;
            for i=1:n
                cost = cost + 1/(p.objective_Ei(evader_positions(:,i),target_position,r,theta));
            end
        end

        function theta = heading_direction(p, evader_positions, target_position, r)
            cost = @(theta) p.objective_fun(evader_positions, target_position, r, theta);
            [theta_min, theta_max, ~, ~] = p.concave_domain(evader_positions, target_position);
            options = optimoptions("fmincon",...
                    "Algorithm","interior-point",...
                    "EnableFeasibilityMode",true,...
                    "SubproblemAlgorithm","cg", "Display","none");
            theta = fmincon(cost,0.5*(theta_min+theta_max),[],[],[],[],theta_min,theta_max,[],options);
        end

        function [velocity, theta] = heading_velocity(p, evader_positions, target_position, r)
            theta = p.heading_direction(evader_positions, target_position, r);
            velocity = p.speed*[cos(theta);sin(theta)];
        end

        function init_start_mtr(p,myev3)
            p.motor1 = motor(myev3,'A');
            p.motor2 = motor(myev3,'B');
            p.motor3 = motor(myev3,'C');
            start(p.motor1);
            start(p.motor2);
            start(p.motor3);
        end

        function stop_mtr(p)
            stop(p.motor1);
            stop(p.motor2);
            stop(p.motor3);
        end

        function set_mtr_speed(p,speed)
            p.motor1.Speed = speed(1);
            p.motor1.Speed = speed(2);
            p.motor1.Speed = speed(3);
        end

        function callback(p,message)
            disp("poda")
            p.pos_vrpn =[message.Pose.Position.X message.Pose.Position.Y message.Pose.Position.Z];
            p.ori_vrpn = [message.Pose.Orientation.X message.Pose.Orientation.Y message.Pose.Orientation.Z message.Pose.Orientation.W]; 
        end
    end
end
