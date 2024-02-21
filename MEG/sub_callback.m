function sub_callback(~, message)
    %exampleHelperROSPoseCallback Subscriber callback function for pose data    
    %   exampleHelperROSPoseCallback(~,MESSAGE) returns no arguments - it instead sets 
    %   global variables to the values of position and orientation that are
    %   received in the ROS message MESSAGE.
    %   
    %   See also ROSPublishAndSubscribeExample.
    
    %   Copyright 2014-2015 The MathWorks, Inc.
    
    % Declare global variables to store position and orientation
    global pos_vrpn
    global orient_vrpn
    
    % Extract position and orientation from the ROS message and assign the
    % data to the global variables.
    pos_vrpn = [message.Pose.Position.X message.Pose.Position.Y message.Pose.Position.Z];
    orient_vrpn = [message.Pose.Orientation.X message.Pose.Orientation.Y message.Pose.Orientation.Z message.Pose.Orientation.W];
end

