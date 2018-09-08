function [x,endPop,bPop,traceInfo] = fmaxga(bounds,evalFN,options,startPop,...
          termFN,termOps,selectFN,selectOps,xOverFNs,xOverOps,mutFNs,mutOps,varargin)
% [X,ENDPOP,BPOP,TRACEINFO] = FMINGA(BOUNDS,EVALFN) runs a Genetic
%  Algorithm (GA) to *MAXIMIZE* the function defined by 'evalFN' (usually a .m
%  file).
%  
%   x            - the best solution found during the course of the run
%   endPop       - the final population 
%   bPop         - a trace of the best population
%   traceInfo    - a matrix of best, mean and variance of the GA for each generation
%
%   bounds       - a matrix which contains the bounds of each variable, i.e.
%                  [var1_low var1_high; var2_low var2_high; ....]
%   evalFN       - the name of the evaluation .m function
%
% [X,ENDPOP,BPOP,TRACEINFO] = FMINGA(BOUNDS,EVALFN,OPTIONS,STARTPOP,
%   TERMFN,TERMOPS,SELECTFN,SELECTOPS,XOVERFNS,XOVEROPS,MUTFNS,MUTOPS,P1,P2,...)
%   allows you to specify your own genetic operators and to pass the
%   problem specific parameters P1,P2,... to the evaluation function. If
%   you want to use the default (given in parenthesis) for a parameter
%   pass the empty matrix [].
%
%   bounds       - a matrix which contains the bounds of each variable, i.e.
%                  [var1_high var1_low; var2_high var2_low; ....]
%
%   evalFN       - the name of the evaluation .m function
%
%   options      - [epsilon floatGA display prec] ([1e-6 1 0 1e-1])
%                    epsilon ... change required to consider two solutions different,
%		     floatGA ... 0 if you want to use the binary GA,
%		                 1 for the float GA;
%                    display ... is 1 to output progress 0 for quiet 
%
%                    prec   ... used to calculate how many bits should be used to
%                               represent one variable
%
%   startPop     - a matrix of initial solutions that is initialized
%                  with initializega.m if [] is given.
%
%   termFN       - name of the .m termination function (['maxGenTerm'])
%
%   termOps      - options string to be passed to the termination function
%                  ([100]).
%
%   selectFN     - name of the .m selection function (['normGeomSelect'])
%
%   selectOpts   - options string to be passed to select after
%                  select(pop,#,opts) ([0.08])
%
%   xOverFNS     - a string containing blank seperated names of Xover.m
%                  files (['arithXover heuristicXover simpleXover']) 
%
%   xOverOps     - A matrix of options to pass to Xover.m files with the
%                  first column being the number of that xOver to perform
%                  ([2 0;2 3;2 0])
%
%   mutFNs       - a string containing blank separated names of mutation.m 
%                  files (['boundaryMutation multiNonUnifMutation ...
%                           nonUnifMutation unifMutation'])
%
%   mutOps       - A matrix of options to pass to mutation.m files with the
%                  first column being the number of that mutation to perform
%                  ([4 0 0;6 100 3;4 100 3;4 0 0])

% Binary and Real-Valued Simulation Evolution for Matlab 
% Copyright (C) 1996 C.R. Houck, J.A. Joines, M.G. Kay 
%
% Modified by Thomas Natschlaeger 1999/01/05
%
% C.R. Houck, J.Joines, and M.Kay. A genetic algorithm for function
% optimization: A Matlab implementation. ACM Transactions on Mathmatical
% Software, Submitted 1996.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 1, or (at your option)
% any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details. A copy of the GNU 
% General Public License can be obtained from the 
% Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

%%$Log: ga.m,v $
%Revision 1.10  1996/02/02  15:03:00  jjoine
% Fixed the ordering of imput arguments in the comments to match
% the actual order in the ga function.
%
%Revision 1.9  1995/08/28  20:01:07  chouck
% Updated initialization parameters, updated mutation parameters to reflect
% b being the third option to the nonuniform mutations
%
%Revision 1.8  1995/08/10  12:59:49  jjoine
%Started Logfile to keep track of revisions
%
 

if nargin < 2
  error('Insufficient arguements. You must supply ''bounds'' and ''evalFN''!');
end
bounds = fliplr(bounds);

if nargin <  3, options  = []; end
if nargin <  4, startPop = []; end
if nargin <  5, termFN   = []; end
if nargin <  6, termOps  = []; end
if nargin <  7, selectFN = []; end
if nargin <  8, selectOps= []; end
if nargin <  9, xOverFNs = []; end
if nargin < 10, xOverOps = []; end
if nargin < 11, mutFNs   = []; end
if nargin < 12, mutOps   = []; end

%
% Set default parameters
%
if isempty(options)
  options = [1e-6 1 0 1e-1];
end

if isempty(startPop) % start Populations
  startPop=initializega(80,fliplr(bounds),evalFN,options,varargin{:});
end

if isempty(termFN)   % termination function
  termOps=[100];
  termFN='maxGenTerm';
end

if isempty(selectFN) % selection function
  selectFN=['normGeomSelect'];
  selectOps=[0.08];
end

if isempty(xOverFNs) % crossover functions
  if options(2)==1 
    % Float GA
    xOverFNs=['arithXover heuristicXover simpleXover'];
    xOverOps=[2 0;2 3;2 0];
  else 
    % Binary GA
    xOverFNs=['simpleXover'];
    xOverOps=[0.6];
  end
end

if isempty(mutFNs)   % muatation functions
  if options(2)==1
    % Float GA
    mutFNs=['boundaryMutation multiNonUnifMutation nonUnifMutation unifMutation'];
    mutOps=[4 0 0;6 termOps(1) 3;4 termOps(1) 3;4 0 0];
  else
    % Binary GA
    mutFNs=['binaryMutation'];
    mutOps=[0.05];
  end
end

if options(2)==0 % if binary version calc number of needed bits
  bits=calcbits(bounds,options(4));
end

%
% check evalFN
%
if any(evalFN<48)
  %
  % Not using a .m file: It is not possible to pass additional parameters!
  %
  if options(2)==1
    % Float GA
    e1str=['x=c1(1:end-1)''; c1(xZomeLength)=' evalFN ';'];
    e2str=['x=c2(1:end-1)''; c2(xZomeLength)=' evalFN ';'];  
  else
    % Binary GA
    e1str=['x=b2f(endPop(j,1:end-1)'',bounds,bits); endPop(j,xZomeLength)=' evalFN ';'];
  end
else
  %
  % If a .m file is used to define the evaluation function, then
  % we pass all additional parameters to evalFN;
  %
  if options(2)==1
    % Float GA
    e1str=['[c1 c1(xZomeLength)]=' evalFN '(c1,gen,varargin{:});'];  
    e2str=['[c2 c2(xZomeLength)]=' evalFN '(c2,gen,varargin{:});'];  
  else 
    % Binary GA
    e1str=['x=b2f(endPop(j,1:end-1),bounds,bits);[x v]=' evalFN '([x endPop(j,end)],gen,varargin{:}); endPop(j,:)=[f2b(x,bounds,bits) v];'];  
  end
end

xOverFNs=parse(xOverFNs);
mutFNs=parse(mutFNs);

xZomeLength  = size(startPop,2); 	%Length of the xzome=numVars+fittness
numVar       = xZomeLength-1; 		%Number of variables
popSize      = size(startPop,1); 	%Number of individuals in the pop
endPop       = zeros(popSize,xZomeLength); %A secondary population matrix
c1           = zeros(1,xZomeLength); 	%An individual
c2           = zeros(1,xZomeLength); 	%An individual
numXOvers    = size(xOverFNs,1); 	%Number of Crossover operators
numMuts      = size(mutFNs,1); 		%Number of Mutation operators
epsilon      = options(1);                 %Threshold for two fittness to differ
oval         = max(startPop(:,xZomeLength)); %Best value in start pop
bFoundIn     = 1; 			%Number of times best has changed
done         = 0;                       %Done with simulated evolution
gen          = 1; 			%Current Generation Number
collectTrace = (nargout>3); 		%Should we collect info every gen
floatGA      = options(2)==1;              %Probabilistic application of ops
display      = options(3);                 %Display progress 

while(~done)
  %Elitist Model
  [bval,bindx] = max(startPop(:,xZomeLength)); %Best of current pop
  best =  startPop(bindx,:);

  if collectTrace
    traceInfo(gen,1)=gen; 		          %current generation
    traceInfo(gen,2)=startPop(bindx,xZomeLength);       %Best fittness
    traceInfo(gen,3)=mean(startPop(:,xZomeLength));     %Avg fittness
    traceInfo(gen,4)=std(startPop(:,xZomeLength)); 
  end
  
  if ( (abs(bval - oval)>epsilon) | (gen==1)) %If we have a new best sol
    if display
      fprintf(1,'\n%d %f\n',gen,bval);          %Update the display
    end
    if floatGA
      bPop(bFoundIn,:)=[gen startPop(bindx,:)]; %Update bPop Matrix
    else
      bPop(bFoundIn,:)=[gen b2f(startPop(bindx,1:numVar),bounds,bits)...
	  startPop(bindx,xZomeLength)];
    end
    bFoundIn=bFoundIn+1;                      %Update number of changes
    oval=bval;                                %Update the best val
  else
    if display
      fprintf(1,'%d ',gen);	              %Otherwise just update num gen
    end
  end
  
  endPop = feval(selectFN,startPop,[gen selectOps]); %Select
  
  if floatGA %Running with the model where the parameters are numbers of ops
    for i=1:numXOvers,
      for j=1:xOverOps(i,1),
	a = round(rand*(popSize-1)+1); 	%Pick a parent
	b = round(rand*(popSize-1)+1); 	%Pick another parent
	xN=deblank(xOverFNs(i,:)); 	%Get the name of crossover function
	[c1 c2] = feval(xN,endPop(a,:),endPop(b,:),bounds,[gen xOverOps(i,:)]);
	
	if c1(1:numVar)==endPop(a,(1:numVar)) %Make sure we created a new 
	  c1(xZomeLength)=endPop(a,xZomeLength); %solution before evaluating
	elseif c1(1:numVar)==endPop(b,(1:numVar))
	  c1(xZomeLength)=endPop(b,xZomeLength);
	else 
	  %[c1(xZomeLength) c1] = feval(evalFN,c1,[gen evalOps]);
	  eval(e1str);
	end
	if c2(1:numVar)==endPop(a,(1:numVar))
	  c2(xZomeLength)=endPop(a,xZomeLength);
	elseif c2(1:numVar)==endPop(b,(1:numVar))
	  c2(xZomeLength)=endPop(b,xZomeLength);
	else 
	  %[c2(xZomeLength) c2] = feval(evalFN,c2,[gen evalOps]);
	  eval(e2str);
	end      
	
	endPop(a,:)=c1;
	endPop(b,:)=c2;
      end
    end
  
    for i=1:numMuts,
      for j=1:mutOps(i,1),
	a = round(rand*(popSize-1)+1);
	c1 = feval(deblank(mutFNs(i,:)),endPop(a,:),bounds,[gen mutOps(i,:)]);
	if c1(1:numVar)==endPop(a,(1:numVar)) 
	  c1(xZomeLength)=endPop(a,xZomeLength);
	else
	  %[c1(xZomeLength) c1] = feval(evalFN,c1,[gen evalOps]);
	  eval(e1str);
	end
	endPop(a,:)=c1;
      end
    end
    
  else %We are running a probabilistic model of genetic operators
    for i=1:numXOvers,
      xN=deblank(xOverFNs(i,:)); 	%Get the name of crossover function
      cp=find(rand(popSize,1)<xOverOps(i,1)==1);
      if rem(size(cp,1),2) cp=cp(1:(size(cp,1)-1)); end
      cp=reshape(cp,size(cp,1)/2,2);
      for j=1:size(cp,1)
	a=cp(j,1); b=cp(j,2); 
	[endPop(a,:) endPop(b,:)] = feval(xN,endPop(a,:),endPop(b,:),...
	  bounds,[gen xOverOps(i,:)]);
      end
    end
    for i=1:numMuts
      mN=deblank(mutFNs(i,:));
      for j=1:popSize
	endPop(j,:) = feval(mN,endPop(j,:),bounds,[gen mutOps(i,:)]);
	eval(e1str);
      end
    end
  end
  
  gen=gen+1;
  done=feval(termFN,[gen termOps],bPop,endPop); %See if the ga is done
  startPop=endPop; 			%Swap the populations
  
  [bval,bindx] = min(startPop(:,xZomeLength)); %Keep the best solution
  startPop(bindx,:) = best; 		%replace it with the worst
end

[bval,bindx] = max(startPop(:,xZomeLength));
if display 
  fprintf(1,'\n%d %f\n',gen,bval);	  
end

x=startPop(bindx,:);
if options(2)==0 %binary
  x=b2f(x,bounds,bits);
  bPop(bFoundIn,:)=[gen b2f(startPop(bindx,1:numVar),bounds,bits)...
      startPop(bindx,xZomeLength)];
else
  bPop(bFoundIn,:)=[gen startPop(bindx,:)];
end

if collectTrace
  traceInfo(gen,1)=gen; 	                  % current generation
  traceInfo(gen,2)=startPop(bindx,xZomeLength);   % Best fittness
  traceInfo(gen,3)=mean(startPop(:,xZomeLength)); % Avg fittness
  traceInfo(gen,4)=std(startPop(:,xZomeLength));  % standard deviation
end






