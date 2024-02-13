function v = Value(xp,xe,alpha)
% Takes in the position a pursuer and an evader and their speed ratios and
% provides the Value function of the 1v1 differential game

if alpha==1
    v=0.5*((norm(xe)^2-norm(xp)^2)/norm(xp-xe));
else if alpha<1
    [xc,rc]=apollonius_parameters(xp,xe,alpha,0);
    v=norm(xc)-rc;
else
    v=norm(xp)-norm(xe)/alpha;
    v=-v;
end

end