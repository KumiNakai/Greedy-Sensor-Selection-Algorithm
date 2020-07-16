%% Main program
%% ///////////////////////////////////////////////////////////////////
% Comments:
% 	Collaborator: Yuji Saito, Keigo Yamada, Taku Nonomura
%                 Kumi Nakai, Takayuki Nagata
% 	Last modified: 2020/7/16
% Nomenclature:
% - Scalars
%   n : Number of degrees of freedom of spatial POD modes (state dimension)
%   p : Number of sensors
%   r : Number of rank for truncated POD
%   m : Number of snaphot (temporal dimension)
% - Matrices
% 	X : Supervising data matrix
% 	Y : Observation matrix
% 	H : Sparse sensor location matrix
% 	U : Spatial POD modes matrix
% 	C : Measurement matrix
% 	Z : POD mode amplitude matrix
%% ===================================================================

clear; close all;
warning('off','all')

%% Selection of Problems ============================================
% num_problem=1; % //Randomized sensor problem//
  num_problem=2; % //NOAA-SST//
% !<NOAA-SST> It takes a long time to obtain the solution in the convex 
% !<NOAA-SST> approximation method and the convex method is commented out 
% !<NOAA_SST> as default setting for reduction of demo time.
%
%% Parameters =======================================================
r = 10;
% pmin = 1;
% pinc = 1;
% pmax = 1;
% ps   = pmin:pinc:pmax;
ps = [5 8 10];
num_ave = 2; % Number of iteration for averaging operation
CNT = 0; % Counter
maxiteration = 200; % Max iteration for convex approximation
% //Randomized sensor problem//
n = 2000;
% //NOAA-SST//
m = 52*10; % 10-years (52weeks/year)
num_video = 2; % maxmum: m

%% Preparation of output directories ================================
workdir   = ('../work');
videodir  = [workdir,'/video'];
sensordir = [workdir,'/sensor_location'];
mkdir(workdir);
mkdir(videodir);
mkdir(sensordir);

%% Randomized sensor problem ========================================
if num_problem == 1

    %% Sensor selection =============================================
    for p = ps
        CNT = CNT+1;
        text = [ num2str(p),' sensor selection started --->' ];
        disp(text);

        %% Average loop =============================================
        for w=1:1:num_ave

            %% Preprocess for Randomized problem ====================
            U = randn(n,r);

            %% Random selection -----------------------------------------
            [time_rand(CNT,w+1), H_rand, sensors_rand] = F_sensor_random(n,p);
            det_rand (CNT,w+1) = F_calc_det  (p,H_rand,U);
            tr_rand  (CNT,w+1) = F_calc_trace(p,H_rand,U);
            eig_rand (CNT,w+1) = F_calc_eigen(p,H_rand,U);

            %% Maximization of det - Convex--------------------------
            %!! This is very time consuming proceduce, We do not recommend to try this
            % [time_convex(CNT,w+1), H_convex, sensors_convex, ...
            %  NT_TOL_cal_convex(CNT,w+1), iter_convex(CNT,w+1)] ...
            %  = F_sensor_convex(U,p,maxiteration);
            % [det_convex (CNT,w+1)] = F_calc_det  (p,H_convex,U);
            % [tr_convex  (CNT,w+1)] = F_calc_trace(p,H_convex,U);
            % [eig_convex (CNT,w+1)] = F_calc_eigen(p,H_convex,U);
            %!! I recommend you use the following dummy values
            %!! if you do not need the solution in the convex approximation in NOAA-SST.        
            time_convex(CNT,w+1) = time_rand(CNT,w+1);
            det_convex (CNT,w+1) = det_rand (CNT,w+1);
            tr_convex  (CNT,w+1) = tr_rand  (CNT,w+1);
            eig_convex (CNT,w+1) = eig_rand (CNT,w+1);
            H_convex=H_rand;
            sensors_convex=sensors_rand;
            
            %% Maximization of row norm - Greedy based on QR --------
            [time_QR(CNT,w+1), H_QR, sensors_QR] = F_sensor_QR(U,p);
            det_QR (CNT,w+1) = F_calc_det  (p,H_QR,U);
            tr_QR  (CNT,w+1) = F_calc_trace(p,H_QR,U);
            eig_QR (CNT,w+1) = F_calc_eigen(p,H_QR,U);

            %% Maximization of det - Greedy -------------------------
            [time_DG(CNT,w+1), H_DG, sensors_DG] = F_sensor_DG(U,p);
            det_DG (CNT,w+1) = F_calc_det  (p,H_DG,U);
            tr_DG  (CNT,w+1) = F_calc_trace(p,H_DG,U);
            eig_DG (CNT,w+1) = F_calc_eigen(p,H_DG,U);

            %% QD (Hybrid of QR and DG) -----------------------------
            [time_QD(CNT,w+1), H_QD, sensors_QD] = F_sensor_QD(U,p);
            det_QD (CNT,w+1) = F_calc_det  (p,H_QD,U);
            tr_QD  (CNT,w+1) = F_calc_trace(p,H_QD,U);
            eig_QD (CNT,w+1) = F_calc_eigen(p,H_QD,U);

            %% Minimization of trace(inverse) -  Greedy -------------
            [time_TG(CNT,w+1), H_TG, sensors_TG] = F_sensor_TG(U,p);
            det_TG (CNT,w+1) = F_calc_det  (p,H_TG,U);
            tr_TG  (CNT,w+1) = F_calc_trace(p,H_TG,U);
            eig_TG (CNT,w+1) = F_calc_eigen(p,H_TG,U);

            %% Maximization of minimum eigen -  Greedy --------------
            [time_EG(CNT,w+1), H_EG, sensors_EG] = F_sensor_EG(U,p);
            det_EG (CNT,w+1) = F_calc_det  (p,H_EG,U);
            tr_EG  (CNT,w+1) = F_calc_trace(p,H_EG,U);
            eig_EG (CNT,w+1) = F_calc_eigen(p,H_EG,U);
        end
        
        %% Averaging ================================================
        [ time_rand, det_rand, tr_rand, eig_rand ]...
        = F_data_ave1( CNT, num_ave, time_rand, det_rand, tr_rand, eig_rand );
        [ time_convex, det_convex, tr_convex, eig_convex ]...
        = F_data_ave1( CNT, num_ave, time_convex, det_convex, tr_convex, eig_convex );
        [ time_QR, det_QR, tr_QR, eig_QR ]...
        = F_data_ave1( CNT, num_ave, time_QR, det_QR, tr_QR, eig_QR );
        [ time_DG, det_DG, tr_DG, eig_DG ]...
        = F_data_ave1( CNT, num_ave, time_DG, det_DG, tr_DG, eig_DG );
        [ time_QD, det_QD, tr_QD, eig_QD ]...
        = F_data_ave1( CNT, num_ave, time_QD, det_QD, tr_QD, eig_QD );
        [ time_TG, det_TG, tr_TG, eig_TG ]...
        = F_data_ave1( CNT, num_ave, time_TG, det_TG, tr_TG, eig_TG );
        [ time_EG, det_EG, tr_EG, eig_EG ]...
        = F_data_ave1( CNT, num_ave, time_EG, det_EG, tr_EG, eig_EG );
    
        NT_TOL_cal_convex(CNT,1)=mean(NT_TOL_cal_convex(CNT,2:w+1));
        iter_convex(CNT,w+1)=mean(iter_convex(CNT,2:w+1));
        
        %% Sensor location ==========================================
        sensor_memo = zeros(p,7);
        sensor_memo(1:p,1) = sensors_rand(1:p)';
        sensor_memo(1:p,2) = sensors_convex(1:p);
        sensor_memo(1:p,3) = sensors_QR(1:p)';
        sensor_memo(1:p,4) = sensors_DG(1:p);
        sensor_memo(1:p,5) = sensors_QD(1:p)';
        sensor_memo(1:p,6) = sensors_TG(1:p)';
        sensor_memo(1:p,7) = sensors_EG(1:p)';
        filename = [workdir, '/sensors_p_', num2str(p), '.mat'];
        save(filename,'sensor_memo');

        text = [ '---> ', num2str(p), ' sensor selection finished!' ];
        disp(text);
    end
end

%% NOAA-SST =========================================================
if num_problem == 2

    %% Preprocces for NOAA-SST ======================================
    text='Readinng/Arranging a NOAA-SST dataset';
    disp(text);
    [Lat, Lon, time, mask, sst]...
    = F_pre_read_NOAA_SST( ['sst.wkmean.1990-present.nc'], ['lsmask.nc'] );
    [Uorg, Sorg, Vorg, Xorg, meansst, n] = F_pre_SVD_NOAA_SST(m, time, mask, sst);
    F_map_original(num_video, Xorg, meansst, mask, time, videodir);
    [U, Error_ave_pod, Error_std_pod]...
    = F_pre_truncatedSVD(r, Xorg, Uorg, Sorg, Vorg, num_video, meansst, mask, time, m, videodir);
    Error_ave_pod = repmat( Error_ave_pod , size(ps,2) );
    text='Complete Reading/Arranging a NOAA-SST dataset!';
    disp(text);

    %% Sensor selection =============================================
    for p = ps
        CNT = CNT+1;
        text = [ num2str(p),' sensor selection started --->' ];
        disp(text);

        %% Random selection -----------------------------------------
        % Average loop
        for w=1:1:num_ave
            [time_rand(CNT,w+1), H_rand, sensors_rand] = F_sensor_random(n,p);
            det_rand(CNT,w+1) = F_calc_det  (p,H_rand,U);
            tr_rand (CNT,w+1) = F_calc_trace(p,H_rand,U);
            eig_rand(CNT,w+1) = F_calc_eigen(p,H_rand,U);
            [Zestimate_rand, Error_rand(CNT,w+1), Error_std_rand(CNT,w+1)] ...
            = F_calc_error(m, Xorg, U, H_rand);
        end
        % Averaging
        [ time_rand, det_rand, tr_rand, eig_rand, Error_rand, Error_std_rand ]...
        = F_data_ave2( CNT, num_ave, time_rand, det_rand, tr_rand, eig_rand, Error_rand, Error_std_rand );

        %% Maximization of det - Convex------------------------------
        %!! This is very time consuming proceduce, We do not recommend to try this
        % [time_convex(CNT,1), H_convex, sensors_convex, ...
        %  NT_TOL_cal_convex(CNT,1), iter_convex(CNT,1)] ...
        %  = F_sensor_convex(U,p,maxiteration);
        % [det_convex(CNT,1)] = F_calc_det(p,H_convex,U);
        % [tr_convex(CNT,1)]  = F_calc_trace(p,H_convex,U);
        % [eig_convex(CNT,1)] = F_calc_eigen(p,H_convex,U);
        %!! I recommend you use the following dummy values 
        %!! if you do not need the solution in the convex approximation in NOAA-SST.
        time_convex(CNT,1) = time_rand(CNT,1);
        det_convex (CNT,1) = det_rand (CNT,1);
        tr_convex  (CNT,1) = tr_rand  (CNT,1);
        eig_convex (CNT,1) = eig_rand (CNT,1);
        H_convex=H_rand;
        sensors_convex=sensors_rand;
        %!!
        [ Zestimate_convex, Error_convex(CNT,1), Error_std_convex(CNT,1) ] ...
        = F_calc_error(m, Xorg, U, H_convex);
        
        %% Maximization of row norm - Greedy based on QR ------------
        [time_QR(CNT,1), H_QR, sensors_QR] = F_sensor_QR(U,p);
        det_QR (CNT,1) = F_calc_det  (p,H_QR,U);
        tr_QR  (CNT,1) = F_calc_trace(p,H_QR,U);
        eig_QR (CNT,1) = F_calc_eigen(p,H_QR,U);
        sensors_QR=sensors_rand;
        [Zestimate_QR, Error_QR(CNT,1), Error_std_QR(CNT,1)] ...
        = F_calc_error(m, Xorg, U, H_QR);

        %% Maximization of det - Greedy -----------------------------
        [time_DG(CNT,1), H_DG, sensors_DG] = F_sensor_DG(U,p);
        det_DG (CNT,1) = F_calc_det  (p,H_DG,U);
        tr_DG  (CNT,1) = F_calc_trace(p,H_DG,U);
        eig_DG (CNT,1) = F_calc_eigen(p,H_DG,U);
        [Zestimate_DG, Error_DG(CNT,1), Error_std_DG(CNT,1)] ...
        = F_calc_error(m, Xorg, U, H_DG);

        %% QD (Hybrid of QR and DG) ---------------------------------
        [time_QD(CNT,1), H_QD, sensors_QD] = F_sensor_QD(U,p);
        det_QD (CNT,1) = F_calc_det  (p,H_QD,U);
        tr_QD  (CNT,1) = F_calc_trace(p,H_QD,U);
        eig_QD (CNT,1) = F_calc_eigen(p,H_QD,U);
        [Zestimate_QD, Error_QD(CNT,1), Error_std_QD(CNT,1)] ...
        = F_calc_error(m, Xorg, U, H_QD);

        %% Minimization of trace(inverse) -  Greedy -----------------
        [time_TG(CNT,1), H_TG, sensors_TG] = F_sensor_TG(U,p);
        det_TG (CNT,1) = F_calc_det  (p,H_TG,U);
        tr_TG  (CNT,1) = F_calc_trace(p,H_TG,U);
        eig_TG (CNT,1) = F_calc_eigen(p,H_TG,U);
        [Zestimate_TG, Error_TG(CNT,1), Error_std_TG(CNT,1)] ...
        = F_calc_error(m, Xorg, U, H_TG);

        %% Maximization of minimum eigen -  Greedy ------------------
        [time_EG(CNT,1), H_EG, sensors_EG] = F_sensor_EG(U,p);
        det_EG (CNT,1) = F_calc_det  (p,H_EG,U);
        tr_EG  (CNT,1) = F_calc_trace(p,H_EG,U);
        eig_EG (CNT,1) = F_calc_eigen(p,H_EG,U);
        [Zestimate_EG, Error_EG(CNT,1), Error_std_EG(CNT,1)] ...
        = F_calc_error(m, Xorg, U, H_EG);
        
        %% Sensor location ==========================================
        sensor_memo = zeros(p,7);
        sensor_memo(1:p,1) = sensors_rand(1:p)';
        sensor_memo(1:p,2) = sensors_convex(1:p);
        sensor_memo(1:p,3) = sensors_QR(1:p)';
        sensor_memo(1:p,4) = sensors_DG(1:p);
        sensor_memo(1:p,5) = sensors_QD(1:p)';
        sensor_memo(1:p,6) = sensors_TG(1:p)';
        sensor_memo(1:p,7) = sensors_EG(1:p)';
        filename = [workdir, '/sensors_p_', num2str(p), '.mat'];
        save(filename,'sensor_memo');

        %% Video ====================================================
        name='rand';
        F_map_reconst(r, num_video, Xorg, meansst, U, mask, time, p, ...
                        sensors_rand, Zestimate_rand, name, videodir, sensordir)
        name='convex';
        F_map_reconst(r, num_video, Xorg, meansst, U, mask, time, p, ...
                        sensors_convex, Zestimate_convex, name, videodir, sensordir)
        name='QR';
        F_map_reconst(r, num_video, Xorg, meansst, U, mask, time, p, ...
                        sensors_QR, Zestimate_QR, name, videodir, sensordir)
        name='DG';
        F_map_reconst(r, num_video, Xorg, meansst, U, mask, time, p, ...
                        sensors_DG, Zestimate_DG, name, videodir, sensordir)
        name='QD';
        F_map_reconst(r, num_video, Xorg, meansst, U, mask, time, p, ...
                        sensors_QD, Zestimate_QD, name, videodir, sensordir)
        name='TG';
        F_map_reconst(r, num_video, Xorg, meansst, U, mask, time, p, ...
                        sensors_TG, Zestimate_TG, name, videodir, sensordir)
        name='EG';
        F_map_reconst(r, num_video, Xorg, meansst, U, mask, time, p, ...
                        sensors_EG, Zestimate_EG, name, videodir, sensordir)
        
        text = [ '---> ', num2str(p), ' sensor selection finished!' ];
        disp(text);
    end
end

%% Data organization ================================================
% Arrange
[time_all] = F_data_arrange1( ps,   CNT, time_rand, time_convex, time_QR,...
                              time_DG,   time_QD,   time_TG,     time_EG );
[det_all]  = F_data_arrange1( ps,   CNT, det_rand,  det_convex,  det_QR, ...
                              det_DG,    det_QD,    det_TG,      det_EG  );
[tr_all]   = F_data_arrange1( ps,   CNT, tr_rand,   tr_convex,   tr_QR,...
                              tr_DG,     tr_QD,     tr_TG,       tr_EG   );
[eig_all]  = F_data_arrange1( ps,   CNT, eig_rand,  eig_convex,  eig_QR,...
                              eig_DG,    eig_QD,    eig_TG,      eig_EG  );
if num_problem == 2
    [Error] = F_data_arrange2( ps,   CNT, ...
                               Error_rand,   Error_std_rand,  ...
                               Error_convex, Error_std_convex,... 
                               Error_QR,     Error_std_QR,    ...
                               Error_DG,     Error_std_DG,    ...
                               Error_QD,     Error_std_QD,    ...
                               Error_TG,     Error_std_TG,    ...
                               Error_EG,     Error_std_EG,    ...
                               Error_ave_pod );
end
% Normalize
[Normalized_det] = F_data_normalize( ps, CNT, det_rand, det_convex, det_QR, ...
                                     det_DG,  det_QD,   det_TG,     det_EG );
[Normalized_tr]  = F_data_normalize( ps, CNT, tr_rand,  tr_convex,  tr_QR,  ...
                                     tr_DG,   tr_QD,    tr_TG,      tr_EG  );
[Normalized_eig] = F_data_normalize( ps, CNT, eig_rand, eig_convex, eig_QR, ...
                                     eig_DG,  eig_QD,   eig_TG,     eig_EG );

%% Save =============================================================
cd(workdir)
save('time.mat','time_all');
save('det.mat','det_all');
save('trace.mat','tr_all');
save('eigen.mat','eig_all');
save('Normalized_det.mat','Normalized_det');
save('Normalized_trace.mat','Normalized_tr');
save('Normalized_eigen.mat','Normalized_eig');
save('time_rand.mat','time_rand');
save('det_rand.mat','det_rand');
save('trace_rand.mat','tr_rand');
save('eigen_rand.mat','eig_rand');
if num_problem == 1
    save('time_convex.mat','time_convex');
    save('time_QR.mat','time_QR');
    save('time_DG.mat','time_DG');
    save('time_QD.mat','time_QD');
    save('time_TG.mat','time_TG');
    save('time_EG.mat','time_EG');
    save('det_convex.mat','det_convex');
    save('det_QR.mat','det_QR');
    save('det_DG.mat','det_DG');
    save('det_QD.mat','det_QD');
    save('det_TG.mat','det_TG');
    save('det_EG.mat','det_EG');
end
if num_problem == 2
    save('Error.mat','Error');
    save('Error_rand.mat','Error_rand');
end
log_convex(1:CNT,1)=ps';
log_convex(1:CNT,2)=NT_TOL_cal_convex(1:CNT,1);
log_convex(1:CNT,3)=ter_convex(1:CNT,1);
save('log_convex.mat','log_convex');

warning('on','all')
disp('Congratulations!');
cd ../src
%% ///////////////////////////////////////////////////////////////////
%% Main program end