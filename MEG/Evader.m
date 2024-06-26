classdef Evader < handle
    properties
        position % 2x1 vector representing current position
        speed % scalar representing maximum speed
        index % integer indexing of each evader
        name % evader name string
        motor1 % motor variables 
        motor2
        motor3
        pos_vrpn
        ori_vrpn
        wheel_radius
        wheel_centre_radius
    end
    methods
        function e = Evader(initPos, speed, index)
            % Constructor, assigns initial values
            e.position = initPos;
            e.speed = speed;
            e.index = index;
            e.name = "evader"+int2str(index);
            e.wheel_radius = 0.025;
            e.wheel_centre_radius = 0.15;
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

        function [velocity, psi] = heading_velocity(e, pursuer_position, target_position, win, alpha)
            if win && alpha==1
                xc = (e.position(1)+pursuer_position(1))/2;
                yc = (e.position(2)+pursuer_position(2))/2;
                m = (e.position(2)-pursuer_position(2))/(e.position(1)-pursuer_position(1));
    
                x_intercept = (m*(yc - target_position(2))+xc+m^2*target_position(1))/(1+m^2);
                y_intercept = target_position(2) + m*(x_intercept - target_position(1));
                
                velocity = [x_intercept - e.position(1), y_intercept - e.position(2)]';
            elseif win && alpha<1
                xc = (e.position(1)-alpha^2*pursuer_position(1))/(1-alpha^2);
                yc = (e.position(2)-alpha^2*pursuer_position(2))/(1-alpha^2);
                rc = (alpha/(1-alpha^2))*norm(pursuer_position-e.position);
                Rc = norm([xc;yc]-target_position);

                x_intercept = target_position(1)+(1-(rc/Rc))*xc;
                y_intercept = target_position(2)+(1-(rc/Rc))*yc;

                velocity = [x_intercept - e.position(1), y_intercept - e.position(2)]';
            else
                sz_target = size(target_position);
                sz_pos = size(e.position);
                if sz_target~=sz_pos
                    if sz_target==fliplr(sz_pos)
                        target_position = target_position';
                    else
                        disp('Wrong target position dimensions.')
                    end
                end
                velocity = target_position - e.position;
            end
            velocity = (e.speed/norm(velocity,2))*(velocity);
            psi = atan2(velocity(2),velocity(1));
        end

        function init_start_mtr(e,myev3)
            e.motor1 = motor(myev3,'A');
            e.motor2 = motor(myev3,'B');
            e.motor3 = motor(myev3,'C');
            start(e.motor1);
            start(e.motor2);
            start(e.motor3);
        end

        function stop_mtr(e)
            stop(e.motor1);
            stop(e.motor2);
            stop(e.motor3);
        end

        function set_mtr_speed(e,speed)
            e.motor1.Speed = speed(1);
            e.motor1.Speed = speed(2);
            e.motor1.Speed = speed(3);
        end

        function callback(e,~,message)
            e.pos_vrpn =[message.Pose.Position.X message.Pose.Position.Y message.Pose.Position.Z];
            e.ori_vrpn = [message.Pose.Orientation.X message.Pose.Orientation.Y message.Pose.Orientation.Z message.Pose.Orientation.W]; 
        end

    end
end
