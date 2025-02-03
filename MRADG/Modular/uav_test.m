% Create UAV scenario
scene = uavScenario("UpdateRate", 10, "StopTime", 20);

% Define UAV platform (Quadcopter)
uav = uavPlatform("MyDrone", scene, "ReferenceFrame", "NED", ...
    "InitialPosition", [0, 0, -10]); % Start at (0,0,-10)

% Define a simple motion model
traj = waypointTrajectory([0, 0, -10; 50, 50, -20], [0, 10]); % A to B

% Attach trajectory to UAV
uav.Trajectory = traj;

% Run the scenario
while scene.CurrentTime < scene.StopTime
    scene.update();
    pose = uav.read(); % Get UAV position
    disp(pose.Position); % Print position in real-time
    pause(0.1);
end