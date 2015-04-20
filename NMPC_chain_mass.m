clear all
close all
clc

%% SIMULATION OPTIONS

EXPORT  = 1;                  % export code for ACADO simulator and solver

COMPILE = 1;                  % compile exported code for ACADO simulator and solver

NRUNS   = 5;                  % run closed-loop simulation NRUNS times and store minimum timings (to minimize OS interference)

ACADOSOLVER = 'qpDUNES';      % 'qpDUNES', 'qpOASES' or 'FORCES' (i.e., FORCES Pro) 

QPCONDENSINGSTEPS = 10;       % number of stages for block condensing (for qpDUNES and FORCES Pro only)

VISUAL = 1;                   % set to 1 to visualize chain of masses (only for the first out of the NRUNS simulations)

WALL = -0.05;                 % position of wall (re-export if changed)

Ts = 0.1;                     % sampling time [s]

N  = 50;                      % number of shooting intervals

NMASS = 4;                    % number of masses in chain (data available from 3 to 6 masses)

INITMODE = 1;                 % 1: initialize lin. at initial condition
                              % 2: initialize lin. at reference

%% Initialization

if QPCONDENSINGSTEPS < 1 || mod(N,QPCONDENSINGSTEPS)~= 0
    error('Invalid block size for given horizon length N.')
end

if strcmp(ACADOSOLVER,'qpDUNES') == 0 && strcmp(ACADOSOLVER,'qpOASES') == 0 && strcmp(ACADOSOLVER,'FORCES') == 0
    error('Invalid solver name.')
end

M  = NMASS - 2;     % number of intermediate masses
NX = (2*M + 1)*3;   % differential states
NU = 3;             % control inputs

DifferentialState xEnd(3,1);    % 3-dimensional position of end point (M+1)
DifferentialState x(M*3,1);     % 3-dimensional position of masses 1, 2, ..., M
DifferentialState v(M*3,1);     % 3-dimensional velocity of masses 1, 2, ..., M
Control u(3,1);

x0 = [0; 0; 0];
L  = 0.033;
D  = 1.0;
m  = 0.03;

%% Differential Equation

A1 = zeros(3,3);
B1 = eye(3,3);

% Compute the spring forces:
g = acado.Expression([0; 0; -9.81]);
f = is(repmat(g, M, 1));
Force = [];
for i = 1:M+1
    if i == 1
        dist = is(x((i-1)*3+1:i*3) - x0);
    elseif( i <= M )
        dist = is(x((i-1)*3+1:i*3) - x((i-2)*3+1:(i-1)*3));
    else
        dist = is(xEnd - x((M-1)*3+1:end));
    end
    
    scale = D/m*(1-L/norm(dist));
    F = is(scale*dist);
    
    Force = [Force; F];
    
    % mass on the right
    if i < M+1
        f((i-1)*3+1:i*3) = f((i-1)*3+1:i*3) - F;
    end
    % mass on the left
    if i > 1
        f((i-2)*3+1:(i-1)*3) = f((i-2)*3+1:(i-1)*3) + F;
    end
end

ode = [ dot(x) == v; ...
        dot(v) == f ];

%% SIMexport

acadoSet('problemname', 'sim');

sim = acado.SIMexport( Ts );
sim.setLinearInput(A1,B1);
sim.setModel(ode);
sim.set( 'INTEGRATOR_TYPE',        'INT_IRK_GL2' );
sim.set( 'NUM_INTEGRATOR_STEPS',        2        );

if EXPORT
    sim.exportCode( 'export_SIM' );
end
if COMPILE
    cd export_SIM
    make_acado_integrator('../integrate_chain')
    cd ..
end

%% MPCexport

acadoSet('problemname', 'mpc');

ocp = acado.OCP( 0.0, N*Ts, N );

rf = [xEnd; x; v; u];
S  = acado.BMatrix(eye(NX+NU));

ocp.minimizeLSQ( S, rf );

rfN = [xEnd; x; v];
SN  = acado.BMatrix(eye(NX));

ocp.minimizeLSQEndTerm( SN, rfN );

ocp.subjectTo( WALL <= [x([2:3:end]); xEnd(2)] <= 10 ); % constraint on y-position TODO upper bound 10?
ocp.subjectTo( -1 <= u <= 1 );                          % box constraints on controls

ocp.setLinearInput(A1,B1);
ocp.setModel(ode);

mpc = acado.OCPexport( ocp );

mpc.set( 'HESSIAN_APPROXIMATION',       'GAUSS_NEWTON'       );
mpc.set( 'DISCRETIZATION_TYPE',         'MULTIPLE_SHOOTING'  );

if strcmp(ACADOSOLVER,'qpDUNES')
    
    mpc.set( 'QP_SOLVER',               'QP_QPDUNES'         );
    mpc.set( 'SPARSE_QP_SOLUTION',      'BLOCK_CONDENSING_N2');
    mpc.set( 'CONDENSING_BLOCK_SIZE',    QPCONDENSINGSTEPS   );
    
elseif strcmp(ACADOSOLVER,'FORCES')
    
    if QPCONDENSINGSTEPS == 1
        mpc.set( 'SPARSE_QP_SOLUTION',  'SPARSE_SOLVER');
    else
        mpc.set( 'SPARSE_QP_SOLUTION',  'BLOCK_CONDENSING_N2');
        mpc.set( 'CONDENSING_BLOCK_SIZE', QPCONDENSINGSTEPS  );
    end
    
    mpc.set( 'QP_SOLVER',               'QP_FORCES'          );
    
elseif strcmp(ACADOSOLVER,'qpOASES')
    
    mpc.set( 'SPARSE_QP_SOLUTION',      'FULL_CONDENSING_N2' );
    mpc.set( 'QP_SOLVER',               'QP_QPOASES'    	 );
    
end

mpc.set( 'INTEGRATOR_TYPE',             'INT_IRK_GL2'        );
mpc.set( 'NUM_INTEGRATOR_STEPS',        2*N                  );

mpc.set( 'MAX_NUM_QP_ITERATIONS', 2000 );
mpc.set( 'PRINTLEVEL', 'LOW' );

if EXPORT
    mpc.exportCode( 'export_MPC' );
end

if COMPILE
    cd export_MPC
    make_acado_solver('../acado_MPCstep')
    cd ..
end

%% SIMULATIONS

minACADOtLog = [];

for iRUNS = 1:NRUNS
    
    eval(['X0  = textread(' '''' 'chain_mass_model_dist_M' num2str(NMASS) '.txt' '''', ', ''''' ');']);
    eval(['ref = textread(' '''' 'chain_mass_model_eq_M' num2str(NMASS) '.txt' '''', ', ''''' ');']);
    
    X0 = [X0(end-3+1:end); X0(1:end-3)];
    ref= [ref(end-3+1:end); ref(1:end-3)];
    
    Xref = repmat(ref.',N+1,1);
    
    Uref = zeros(N,NU);
    
    if INITMODE == 1
        input.x = repmat(X0.',N+1,1);
    elseif INITMODE == 2
        input.x = repmat(ref.',N+1,1);
    else
        error('wrong initialization flag INITMODE')
    end
    
    input.u = Uref;
    
    input.y = [Xref(1:N,:) Uref];
    input.yN = ref;
       
    display('------------------------------------------------------------------')
    display('               Simulation Loop'                                    )
    display('------------------------------------------------------------------')
    
    iter = 0; 
    time = 0;
    Tf   = 4;
        
    controls_MPC = [];
    state_sim    = X0.';
    
    ACADOtLog    = [];  % log timings of solvers
    ACADOoutputs = {};  % log all ACADO outputs
    ACADOnIter   = [];  % log iterations (if available)
    
    if VISUAL && iRUNS == 1
        visualize;
    end
    
    % weights
    input.W  = blkdiag(15*eye(3), 10*eye(length(x)), eye(length(v)), 0.05*eye(3));
    input.WN = blkdiag(15*eye(3), 10*eye(length(x)), eye(length(v)));
    
    if NRUNS > 1
        % Re-initialize ACADO solver in case of multiple runs
        input.x0      = state_sim';
        input.control = 1;
        output        = acado_MPCstep(input);
        input.control = 0;
    end
    
    while time(end) < Tf
        
        % Solve NMPC OCP with ACADO
        input.x0 = state_sim(end,:).';
        output   = acado_MPCstep(input);
        
        if output.info.status ~= 1 &&  output.info.status ~= 0
            keyboard
        end
        
        ACADOnIter = [ACADOnIter output.info.nIterations];
        ACADOoutputs{end+1} = output;      
        ACADOtLog = [ACADOtLog; 1000*output.info.cpuTime];
        
        % Save the MPC step
        controls_MPC = [controls_MPC; output.u(1,:)];
        
        % Shift trajectories
        input.x = [output.x(2:end,:); output.x(end,:)];
        input.u = [output.u(2:end,:); output.u(end,:)];
                
        % Simulate system
        sim_input.x = state_sim(end,:).';
        sim_input.u = output.u(1,:).';
        [states,outputs] = integrate_chain(sim_input);
        state_sim  = [state_sim; states.value'];
        
        iter = iter + 1;
        nextTime = iter*Ts;
        
        disp(['current time: ' num2str(nextTime) '   ' char(9) ' (RTI step: ' num2str(output.info.cpuTime*1e3) ' ms)'])
        
        time = [time nextTime];
        
        if VISUAL && iRUNS == 1
            visualize;
        end
        
    end
    
    % delete very first QP from measurements TODO
    ACADOtLog(1,:) = [];
    
    % delete last time instant
    time(end) = [];
    
    if iRUNS > 1
        minACADOtLog = min(ACADOtLog, minACADOtLog);
    else
        minACADOtLog = ACADOtLog;
    end
        
end


figure
plot(time(2:end),minACADOtLog,'lineWidth',2)
title('Timings of RTI scheme in closed loop','fontSize',16,'fontWeight','Normal')
xlabel('Time [s]','fontSize',16)
ylabel('CPU time [ms]','fontSize',16)
leg = legend(ACADOSOLVER);
leg.FontSize = 16;
set(gca,'fontSize',16)



