function [time_convex, H, sensors, NT_TOL_cal, iter]=F_sensor_convex(U, p, maxiteration)

    [n,~]=size(U);
    tic;
    [sensors, ~, NT_TOL_cal, iter]=F_sensor_convex_sub(U, p, maxiteration);
    time_convex=toc;
    [H]=F_calc_sensormatrix(p, n, sensors);
    
end
