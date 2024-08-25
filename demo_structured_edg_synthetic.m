% ------------------------------------------------------------------------------------
% File Name: demo_structured_edg_protein.m
%
% Description:
% This script demonstrates the AB_completion.m script. We test it on synthetic data,
% with idealized sampling model. Below, we provide a concise 
% description of the problem and highlight the key parameters involved.
%
% Problem Overview:
% We consider a scenario with a total of p = n + m nodes, categorized into two types:
% - "Anchor Nodes" (m nodes): These nodes have some known distances between them.
% Remark: Anchor nodes typically refers to nodes whose positions are known.
% In our case, we don't know any position information about the anchor nodes.
% Rather, we have partial distance information between these nodes.
% 
% - Mobile Nodes (n nodes): We have no distance information between these nodes.
%
% Assumptions:
% 1. No distance information is available between the mobile nodes.
% 2. The distance matrix for the anchor nodes is sampled according to the Bernoulli model, 
%    where (gamma) denotes the probability of selecting each entry.
% 3. From each mobile node, we observe the distances to alpha-1 anchor nodes
%    that are sampled uniformly at random.
% 4. There is a central anchor node with known distances to all other nodes.
%
% The squared distance matrix D has the following block structure:
%
%          D = [E   F
%               F^T G]
%
% Where:
% - E represents the squared distance matrix between anchor nodes.
% - F represents the squared distance matrix between anchor and mobile nodes.
% - G represents the squared distance matrix between mobile nodes.
%
% Objective:
% The goal is to estimate the configuration of the nodes given the parameters m, n, 
% gamma, and alpha. Additionally, the method of selecting the m anchor nodes plays 
% a critical role; in this script, we select m anchors uniformly spaced within the 
% range [1,p]
%
% For a detailed explanation of the algorithm and its theoretical foundations, please 
% refer to our publication:
% 
% Citation:
% Samuel Lichtenberg and Abiy Tasissa, “Localization from Structured Distance Matrices 
% via Low-Rank Matrix Recovery”, IEEE Transactions on Information Theory, 2024.
% arXiv link: https://arxiv.org/pdf/2311.18076
%
% ---------------------------------------------------------------------------------------
% Set up problem parameters
% r: embedding dimension
% tot_samples: number of samples in a column of F. This does not include
% the central node. tot_samples = alpha-1;
% gamma = percentage of sampling for E block
r = 2;
alpha = 10;
tot_samples = alpha-1;
gamma = 0.1;
% m: number of anchors
m = 50;
% n: number of mobile nodes
n = 500;
% Create synthetic data. 
P = randn(r,m+n);
p = length(P);
% Center points and set-up anchors
P(:,1:m) = P(:,1:m)-mean(P(:,1:m),2);
% Form the Gram matrix
X = P(:,1:m);
Y = P(:,m+1:end);
A = X'*X;
B = X'*Y;
C = Y'*Y;
G = [A B ; B' C];
% Squared distance matrix
dist = squareform(pdist(P'));
D = dist.*dist;
% Number of trials
num_trials = 1;
% Store the RMSE between the true point cloud and the reconstructed point cloud
RMSE_modified = zeros(num_trials,1);
for q = 1:num_trials
    q   
    % Construct the weight matrix for F
    Weight_F1=zeros(m-1,n); 
    for i = 1:n
        ridx= randperm(m-1);
        Weight_F1(ridx(1:tot_samples),i)=1;
    end
    Weight_F  = [Weight_F1; ones(1,n)];
    % Construct the weight matrix for E
    Weight_E1=rand(m-1,m-1);
    Weight_E1(Weight_E1>1-gamma)=1;
    Weight_E1(Weight_E1<1)=0;
    Weight_E1(Weight_E1>0)=1;
    for i=1:m-1
    Weight_E1(i,i)=1;
    for j=i+1:m-1
        Weight_E1(i,j)=Weight_E1(j,i);
    end
    end
    Weight_E = [Weight_E1 ones(m-1,1) ; ones(1,m)];
    % Choice of algorithm
    G_modified_nys = AB_completion(D,Weight_E,Weight_F,m,n);
    % Find the reconstructed point
    [V,Lam] = eigs(G_modified_nys,r,'lm');
    lam = diag(Lam);
    [lam,IJ] = sort(lam,'descend');
    V = V(:,IJ);
    Pt_modified = real(V(:,1:r)*diag(sqrt(lam(1:r))));
    % RMSE error for algorithm
    [RMSE_modified(q),nystrom_aligned,~] =  Compute_RMSE(P',Pt_modified);        
    fprintf("The RMSE error for modified Nystrom is %6f\n", RMSE_modified(q))
 
end