function [xc,rc] = apollonius_parameters(xp,xe,alpha,draw)
% This function takes in the starting position of the pursuer and the
% evader and returns the center and radius of the Apollonius circle

if alpha<1
    xc=(xe-alpha^2*xp)/(1-alpha^2);
    rc=(alpha/(1-alpha^2))*norm(xp-xe);
end
if alpha>1
    xc=(alpha^2*xp-xe)/(alpha^2-1);
    rc=(alpha/(alpha^2-1))*norm(xp-xe);
end

if draw==1
    [X,Y,Z]=sphere;
    X=rc*X;Y=rc*Y;Z=rc*Z;
    surf(X+xc(1),Y+xc(2),Z+xc(3),FaceColor="flat",FaceAlpha=0.1,EdgeAlpha=0.1)
    axis equal
end

end