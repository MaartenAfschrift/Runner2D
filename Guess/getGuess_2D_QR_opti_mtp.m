% This script provides an inital guess for the design variables. The guess is
% quasi-random (QR). We set constant values to the muscle variables, arm
% variables and most joint variables. We only ensure that the distance traveled
% is non-null. The model is moving forward at a constant speed and is standing
% on the ground. We use a pre-defined guess for the final time..
%
% Author: Antoine Falisse
% Date: 9/9/2019
% 
function guess = getGuess_2D_QR_opti_mtp(N,d,nq,NMuscle,scaling,v_tgt,jointi,varargin)

if isempty(varargin)
   PelvisTy = 0.9385;
else
   PelvisTy = varargin{1}; 
end
%% Final time
guess.tf = 0.6;

%% Qs
% The model is moving forward but with a standing position (Qs=0)
guess.Qs = zeros(N,nq.all);
guess.Qs(:,jointi.pelvis.tx) = linspace(0,guess.tf*v_tgt,N);
% The model is standing on the ground
guess.Qs(:,jointi.pelvis.ty) = PelvisTy;

%% Qdots
guess.Qdots = zeros(N,nq.all);
% The model is moving forward with a constant speed
guess.Qdots(:,jointi.pelvis.tx) = v_tgt;
% Qs and Qdots are intertwined
guess.QsQdots = zeros(N,2*nq.all);
guess.QsQdots(:,1:2:end) = guess.Qs;
guess.QsQdots(:,2:2:end) = guess.Qdots;

%% Qdotdots
guess.Qdotdots = zeros(N,nq.all);

%% Muscle variables
guess.a = 0.1*ones(N,NMuscle);
guess.vA = 0.01*ones(N,NMuscle);
guess.FTtilde = 0.1*ones(N,NMuscle);
guess.dFTtilde = 0.01*ones(N,NMuscle);

%% Back variables
guess.a_b = 0.1*ones(N,nq.trunk);
guess.e_b = 0.1*ones(N,nq.trunk);

%% Scaling
guess.QsQdots   = guess.QsQdots./repmat(scaling.QsQdots,N,1);
guess.Qdotdots  = guess.Qdotdots./repmat(scaling.Qdotdots,N,1);
guess.a         = (guess.a)./repmat(scaling.a,N,size(guess.a,2));
guess.FTtilde   = (guess.FTtilde)./repmat(scaling.FTtilde,N,1);
guess.vA        = (guess.vA)./repmat(scaling.vA,N,size(guess.vA,2));
guess.dFTtilde  = (guess.dFTtilde)./repmat(scaling.dFTtilde,N,...
    size(guess.dFTtilde,2));

%% Update guess information for use with opti implementation
% (also evaluate at collocation points

% add values for last mesh point (states)
guess.QsQdots   = [guess.QsQdots; guess.QsQdots(end,:)];
guess.a         = [guess.a; guess.a(end,:)];
guess.FTtilde   = [guess.FTtilde; guess.FTtilde(end,:)];
guess.a_b       = [guess.a_b; guess.a_b(end,1)];

% approximate values at the collocation points
gmesh = linspace(0,1,N+1);
gcoll = linspace(0,1,N*d);
guess.QsQdots_c = interp1(gmesh, guess.QsQdots,gcoll);
guess.a_c       = interp1(gmesh, guess.a,gcoll);
guess.FTtilde_c = interp1(gmesh, guess.FTtilde,gcoll);
guess.a_b_c = interp1(gmesh, guess.a_b,gcoll);


end
