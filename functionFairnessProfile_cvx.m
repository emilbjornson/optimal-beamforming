function [finalInterval,WBestBeamforming,nbrOfEvaluations] = functionFairnessProfile_cvx(H,D,Qsqrt,q,delta,lowerPoint,upperPoint,specialMode,specialParam)
%Solves the fairness-profile optimization (FPO) problem in Example 2.8 in
%the book:
%
%Emil Björnson, Eduard Jorswieck, “Optimal Resource Allocation in
%Coordinated Multi-Cell Systems,” Foundations and Trends in Communications
%and Information Theory, vol. 9, no. 2-3, pp. 113-381, 2013.
%
%This is a generalization of max-min fairness optimization. The FPO
%problem searches on a line between LOWERPOINT and UPPERPOINT and finds
%the intersection between the line and the Pareto boundary of the rate
%region. The search is based on bisection and solving a power minimization
%under QoS requirements; for example, as in (2.30).
%
%The FPO problem is equivalent to finding (g_1,...,) that solves
%
%maximize min_k ( g_k - lowerPoint(k) ) / ( upperPoint(k) - lowerPoint(k) )
%
%subject to     g_k >= lowerPoint(k) for all users k,
%               Power constraints.
%
%This problem is quasi-convex, meaning that it can be solved as a short
%sequence of convex optimization problem. The computational complexity is
%therefore polynomial in the number of users, antennas, and power
%constraints. The implementation can, at least, handle 30 users, 50
%antennas, and 50 power constraints.
%
%This is version 1.1. (Last edited: 2014-03-26)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%original article listed above.
%
%The implementation utilizes and requires CVX: http://cvxr.com/
%
%INPUT:
%H           = Kr x Kt*Nt matrix with row index for receiver and column
%              index transmit antennas
%D           = Kt*Nt x Kt*Nt x Kr diagonal matrix. Element (j,j,k) is one if
%              j:th transmit antenna can transmit to user k and zero otherwise
%Qsqrt       = N x N x L matrix with matrix-square roots of the L weighting
%              matrices for the L power constraints
%q           = Limits of the L power constraints
%delta       = Accuracy of the final solution. The algorithm terminates when
%                 norm(upperPoint - lowerPoint) <= delta
%lowerPoint  = Start point of the line (must be inside the rate region)
%upperPoint  = End point of the line (must be outside of the rate region)
%specialMode = (Optional) Consider different non-idealities in the system
%              model. Normally we have specialMode==0. 
%              Transceiver impairments is given by specialMode==1.
%specialParam = (Optional) Additional parameters for each specialMode.
%               specialMode==1: 2 x 1 vector with EVM at the transmit and 
%               receive antennas, respectively. These are used along with
%               linear distortion functions.
%
%OUTPUT:
%finalInterval    = Kr x 2 matrix with lowerPoint and upperPoint at
%                   termination
%WBestBeamforming = Kt*Nt x Kr matrix with beamforming that achieves the 
%                   lower point in the final interval
%nbrOfEvaluations = Number of times that the convex subproblem (power
%                   minimization under QoS requirements) is solved

if nargin < 8
    specialMode = 0;
end

Kr = size(H,1); %Number of users
L = size(Qsqrt,3); %Number of power constraints


%Pre-allocation of matrix for storing optimal beamforming
WBestBeamforming = [];

%Count the number of feasibility problem solved
nbrOfEvaluations = 0;

%%Part 1: Solve the problem by bisection.

%Solve the problem by bisection - iterate until different between
%current lower and upper point
while norm(upperPoint - lowerPoint) > delta
    
    candidatePoint = (lowerPoint+upperPoint)/2; %Compute the midpoint at the line
    
    gammavar = 2.^(candidatePoint)-1; %Transform midpoint into SINR requirements
    
    %Check the feasibility at the midpoint by solving a feasibility
    %problem. Different feasibility problems are solved depending on the
    %mode. 
    if specialMode == 0 %Ideal case
        [feasible,Wcandidate] = functionFeasibilityProblem_cvx(H,D,Qsqrt,q,gammavar);
    elseif specialMode == 1 %Transceiver impairments
        [feasible,Wcandidate] = functionFeasibilityProblem_Impairments_cvx(H,D,Qsqrt,q,gammavar,specialParam);
    end
    
    %If the problem was feasible, then replace lowerPoint with
    %candidatePoint and store W as current best solution.
    if feasible
        lowerPoint = candidatePoint;
        WBestBeamforming = Wcandidate;
    else
        %If the problem was not feasible,then replace upperPoint with candidatePoint
        upperPoint = candidatePoint;
    end
    
    %Increase the number of function evaluations
    nbrOfEvaluations = nbrOfEvaluations+1;
end




%%Part 2: Prepare the achieved solution for output

%If the midpoints analyzed by the algorithm have never been feasible,
%then obtain a feasible beamforming solution using the lowerPoint. This
%happens when delta is too large or when the optimal point is very
%close to lowerPoint.
if isempty(WBestBeamforming)
    gammavar = 2.^(lowerPoint)-1;
    
    if specialMode == 0 %Ideal case
        [feasible,Wcandidate] = functionFeasibilityProblem_cvx(H,D,Qsqrt,q,gammavar);
    elseif specialMode == 1 %Transceiver impairments
        [feasible,Wcandidate] = functionFeasibilityProblem_Impairments_cvx(H,D,Qsqrt,q,gammavar,specialParam);
    end
    
    if feasible
        WBestBeamforming = Wcandidate;
    else
        %The algorithm requires that the start point is inside of the
        %rate region, which is not the case if we end up here.
        error('Fairness-profile optimization problem is infeasible');
    end
end

%Prepare for output, depending on scenario mode
if specialMode == 0 %Ideal case
    
    %Change scaling of achieved beamforming to satisfy at least one power
    %constraint with equality (based on Theorem 1.2).
    scaling = zeros(L,1);
    for l = 1:L
        scaling(l) = norm(Qsqrt(:,:,l)*WBestBeamforming,'fro').^2/q(l);
    end
    
    %Scale beamforming to satisfy at least one power constraint with equality
    WBestBeamforming = WBestBeamforming/sqrt(max(scaling));
    
    
    %Compute the rates that are actually achieved by WBestBeamforming
    channelGains = abs(H*WBestBeamforming).^2;
    signalGains = diag(channelGains);
    interferenceGains = sum(channelGains,2)-signalGains;
    rates = log2(1+signalGains./(1+interferenceGains));
    
    %Store the final interval between lower and upper point
    if sum(rates>lowerPoint) == Kr
        finalInterval = [rates upperPoint];
    else
        finalInterval = [lowerPoint upperPoint];
    end
    
elseif specialMode == 1 %Transceiver impairments
    
    %Store the final interval between lower and upper point
    finalInterval = [lowerPoint upperPoint];
    
end
