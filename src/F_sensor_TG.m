function [time, H, sensors]=F_sensor_TG(U, p)
    % objective function: trace[inv(CCT)], trace[inv(CTC)]

    [n,~]=size(U);
    tic;
    [sensors]=F_sensor_TG_calc_trace(U, p);
    time=toc;
    [H]=F_calc_sensormatrix(p, n, sensors);
    
end