function [feasible,Wsolution] = functionFeasibilityProblem_cvx(H,D,Qsqrt,q,gammavar)
%Solves the feasibility problem with quality-of-service (QoS) constraints
%in (2.29). The implementation is based on the alternative formulation in
%(2.30): power minimization under QoS requirements.
%
%The references to equations refer to the following book:
%
%Emil Björnson, Eduard Jorswieck, “Optimal Resource Allocation in
%Coordinated Multi-Cell Systems,” Foundations and Trends in Communications
%and Information Theory, vol. 9, no. 2-3, pp. 113-381, 2013.
%
%The power minimization under QoS requirements is
%
%minimize   betavar
%subject to SINR_k >= gammavar(k) for all users k,
%           Power constraints scaled by betavar.
%
%If this optimization problem is feasible and betavar<=1, then the 
%feasibility problem with QoS constraints is also feasible.
%
%This optimization problem is convex. The computational complexity is
%therefore polynomial in the number of users, antennas, and power
%constraints. The implementation can, at least, handle 30 users, 50
%antennas, and 50 power constraints.
%
%This is version 1.2. (Last edited: 2015-08-24)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%original article listed above.
%
%The implementation utilizes and requires CVX: http://cvxr.com/
%
%INPUT:
%H          = Kr x Kt*Nt matrix with row index for receiver and column
%             index transmit antennas
%D          = Kt*Nt x Kt*Nt x Kr diagonal matrix. Element (j,j,k) is one if
%             j:th transmit antenna can transmit to user k and zero otherwise
%Qsqrt      = N x N x L matrix with matrix-square roots of the L weighting 
%             matrices for the L power constraints
%q          = Limits of the L power constraints
%gammavar   = Kr x 1 vector with SINR constraints for all users.
%
%OUTPUT:
%feasible  = This variable is feasible=true if the feasibility problem is
%            feasible. Otherwise we have feasible=false.
%Wsolution = Kt*Nt x Kr matrix with beamforming achieved by the power
%            minimization problem. This matrix is empty if this problem is
%            infeasible.


Kr = size(H,1); %Number of users
N = size(H,2); %Number of transmit antennas (in total)
L = size(Qsqrt,3); %Number of power constraints


%Solve the power minimization under QoS requirements problem using CVX
cvx_begin
cvx_quiet(true); % This suppresses screen output from the solver

variable W(N,Kr) complex;  %Variable for N x Kr beamforming matrix
variable betavar %Scaling parameter for power constraints

minimize betavar %Minimize the power indirectly by scaling power constraints

subject to

%SINR constraints (Kr constraints)
for k = 1:Kr
    
    %Channels of the signal intended for user i when it reaches user k
    hkD = zeros(Kr,N);
    for i = 1:Kr
        hkD(i,:) = H(k,:)*D(:,:,i);
    end
    
    imag(hkD(k,:)*W(:,k)) == 0; %Useful link is assumed to be real-valued
    
    %SOCP formulation for the SINR constraint of user k
    real(hkD(k,:)*W(:,k)) >= sqrt(gammavar(k))*norm([1 hkD(k,:)*W(:,[1:k-1 k+1:Kr])  ]);
end

%Power constraints (L constraints) scaled by the variable betavar
for l = 1:L
    norm(Qsqrt(:,:,l)*W,'fro') <= betavar*sqrt(q(l));
end

betavar >= 0; %Power constraints must be positive

cvx_end


%Analyze result and prepare the output variables.
if isempty(strfind(cvx_status,'Solved')) %Both power minimization problem and feasibility problem are infeasible.
    feasible = false;
    Wsolution = [];
elseif betavar>1 %Only power minimization problem is feasible.
    feasible = false;
    Wsolution = W;
else %Both power minimization problem and feasibility problem are feasible.
    feasible = true;
    Wsolution = W;
end
