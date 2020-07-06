function [isensors,zhat]=F_sensor_convex_sub(U, p, maxiteration)

    [n,r]=size(U);
    UU=[];
    UU(:,:,1)=U(1:n,1:r); 
    [zhat, ~, ~] = F_sensor_convex_approxnt_vec(UU, p, maxiteration);
    isensors = find(zhat);

end