clc;
clear;
close all;

%% Problem Definition

global NFE;
NFE=0;                  %***Check***%


% CostFunction=@(x,y) CostFunction(x,y);      % Cost Function

nVar=6;             % Number of Decision Variables

VarSize=[1 nVar];   % Size of Decision Variables Matrix

VarMin=1;           % Lower Bound of Variables %%changing%%
VarMax=4;          % Upper Bound of Variables %%changing%%

% Number of Objective Functions
nObj=3;

% Inputting the Structural Parameters
InputParameters = [10,20,5]; % radius (A), pitch angle (degree),  diameter(A)

%% NSGA-II Parameters

MaxIt=100;      % Maximum Number of Iterations

nPop=10;        % Population Size    %***Check***%

pCrossover=0.7;                         % Crossover Percentage
nCrossover=2*round(pCrossover*nPop/2);  % Number of Parnets (Offsprings)

pMutation=0.4;                          % Mutation Percentage
nMutation=round(pMutation*nPop);        % Number of Mutants

mu=0.166;                    % Mutation Rate

%sigma=0.1*(VarMax-VarMin);  % Mutation Step Size


%% Initialization

empty_individual.Position=[];
empty_individual.Cost=[];
empty_individual.Rank=[];
empty_individual.DominationSet=[];
empty_individual.DominatedCount=[];
empty_individual.CrowdingDistance=[];
empty_individualt.Position=[];
empty_individualt.SimulationNo=[];
empty_individualt.radius=[];
empty_individualt.pitchangle=[];
empty_individualt.diameter=[];
empty_individualt.lengthccnt=[];
empty_individualt.pitchl=[];
empty_individualt.numberofturns=[];


pop=repmat(empty_individual,nPop,1);
poptotal = repmat(empty_individualt,1,1);

fileID1=fopen('good.txt','w');
fileID2=fopen('bad.txt','w');

for i=1:nPop                    %***Check***%
    
    while true                %checking the structures and eliminating the bad structures
        for k=1:nVar          %creating the positions
            x(k)=randi([VarMin VarMax]);
            x(5)=randi([1 2]);
        end
        if ismember(x,poptotal.Position)
            continue
        else
            [y,z,ytemp]=teststructure(x,InputParameters);
            if y==1
                Writing_Bad(x);
                continue
            end
        end
        break
    end
    fileID1=fopen('good.txt','a+');
    fprintf(fileID1,'%d\n',x);
    fclose all;
    pop(i).Position=x;
   pop(i).Cost = z;
   poptotal=[poptotal;ytemp];
end

% the initial population from last simulation
% popinput=load('popinput');

% for i=100:nPop                 %***Check***%
%     [pop(i).Cost,ytemp]=CostFunction(popinput.pop(i-99).Position,InputParameters);
%     poptotal=[poptotal;ytemp];
% end

% generating poptotal again
% poptotalinput=load('poptotalinput');
% 
% for i=100:nPop                 %***Check***%
%     poptotal(i).Position=poptotalinput.poptotal(i-99).Position;
%     poptotal(i).SimulationNo=poptotalinput.poptotal(i-99).SimulationNo;
%     poptotal(i).YoungModulus=poptotalinput.poptotal(i-99).YoungModulus;
%     poptotal(i).Toughness=poptotalinput.poptotal(i-99).Toughness;
%     poptotal(i).YieldStress=poptotalinput.poptotal(i-99).YieldStress;
%     poptotal(i).YieldStrain=poptotalinput.poptotal(i).YieldStrain;
%     poptotal(i).Stressend=poptotalinput.poptotal(i).Stressend;
%     poptotal(i).Strainend=poptotalinput.poptotal(i).Strainend;
% end


% Non-Dominated Sorting
[pop F]=NonDominatedSorting(pop);

% Calculate Crowding Distance
pop=CalcCrowdingDistance(pop,F);

% Sort Population
[pop F]=SortPopulation(pop);

%Array to hold Number of Function Evaluation
nfe=zeros(MaxIt,1);

%% NSGA-II Main Loop

for it=1:MaxIt
    
    % Crossover
    popc=repmat(empty_individual,nCrossover/2,2);
    for k=1:nCrossover/2
        while true
        ss=randperm(nPop,2);      %two different random numbers
        p1=pop(ss(1));
        p2=pop(ss(2));
        [popc(k,1).Position,popc(k,2).Position]=Crossover(p1.Position,p2.Position);
        if (teststructure(popc(k,1).Position)==1)
            Writing_Bad(popc(k,1).Position);
            while true
                ss=randperm(nPop,2);
                 p1=pop(ss(1));
                 p2=pop(ss(2));
                 [popc(k,1).Position,~]=Crossover(p1.Position,p2.Position);
                 if (teststructure(popc(k,1).Position)==1)
                     Writing_Bad(popc(k,1).Position);
                     continue
                 else
                     [popc(k,1).Cost,ytemp]=CostFunction(popc(k,1).Position,InputParameters);
                     poptotal=[poptotal;ytemp];
                     fileID1=fopen('good.txt','a+');
                     fprintf(fileID1,'\n');
                     fprintf(fileID1,'%d',popc(k,1).Position);
                     fclose all;
                 end
                 break
            end
        else
            [popc(k,1).Cost,ytemp]=CostFunction(popc(k,1).Position,InputParameters);
            poptotal=[poptotal;ytemp];
            fileID1=fopen('good.txt','a+');
            fprintf(fileID1,'\n');
            fprintf(fileID1,'%d',popc(k,1).Position);
            fclose all;
        end
        
        if (teststructure(popc(k,2).Position)==1)
            Writing_Bad(popc(k,2).Position);
            while true
                ss=randperm(nPop,2);
                 p1=pop(ss(1));
                 p2=pop(ss(2));
                 [~,popc(k,2).Position]=Crossover(p1.Position,p2.Position);
                 if (teststructure(popc(k,2).Position)==1)
                     Writing_Bad(popc(k,2).Position);
                     continue
                 else
                     [popc(k,2).Cost,ytemp]=CostFunction(popc(k,2).Position,InputParameters);
                     poptotal=[poptotal;ytemp];
                     fileID1=fopen('good.txt','a+');
                     fprintf(fileID1,'\n');
                     fprintf(fileID1,'%d',popc(k,2).Position);
                     fclose all;
                 end
                 break
            end
        else
            [popc(k,2).Cost,ytemp]=CostFunction(popc(k,2).Position,InputParameters);
            poptotal=[poptotal;ytemp];
            fileID1=fopen('good.txt','a+');
            fprintf(fileID1,'\n');
            fprintf(fileID1,'%d',popc(k,2).Position);
            fclose all;
        end
        save all;
        
    end
    popc=popc(:);
    
    % Mutation
    popm=repmat(empty_individual,nMutation,1);
    for k=1:nMutation
        while true
            i=randi([1 nPop]);
            p=pop(i);
        
            popm(k).Position=MutateAM2(p.Position,mu,VarMax,VarMin,nVar);
        
            if ((teststructure(popm(k).Position))==1)
                Writing_Bad(popm(k).Position);
                continue
            end
            break
        end
            [popm(k).Cost,ytemp]=CostFunction(popm(k).Position,InputParameters);
            poptotal=[poptotal;ytemp];
            fileID1=fopen('good.txt','a+');
            fprintf(fileID1,'\n');
            fprintf(fileID1,'%d',popm(k).Position);
            fclose all;
            save all;
    end
    F2=pop(F{1});	
    % Merge
    pop=[pop
         popc
         popm];
     
    % Non-Dominated Sorting
    [pop F]=NonDominatedSorting(pop);

    % Calculate Crowding Distance
    pop=CalcCrowdingDistance(pop,F);

    % Sort Population
    [pop F]=SortPopulation(pop);
    
    % Truncate
    pop=pop(1:nPop);
    
    % Non-Dominated Sorting
    [pop F]=NonDominatedSorting(pop);

    % Calculate Crowding Distance
    pop=CalcCrowdingDistance(pop,F);

    % Sort Population
    [pop F]=SortPopulation(pop);
    
    % Store F1
    F1=pop(F{1});
    
    % Store NFE
    nfe(it)=NFE;
    save all;
    
    % Show Iteration Information
    disp(['Iteration ' num2str(it) ': NFE = ' num2str(nfe(it))...
        ': Number of F1 Members = ' num2str(numel(F1))]);
        
    % Plot F1 Costs
    figure(1);
    PlotCosts(F1,F2);
    
end
%     fclose(fileID1);
%     fclose(fileID2);

%% Results