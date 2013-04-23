if~(exist('restart'))   % in case we want to start from last iteration
close all
clear all
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Link to wherever pmpack is
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  addpath ../pmpack/ ../pmpack/test

  global totalevals

  totalevals = 0;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Define the ACES input file
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %input_file = @sample_input_decomp;
  %input_file = @sample_input_discont2;
  %input_file = @sample_input_dimAdaptive; % 2D dim adaptive test case
  %input_file = @sample_input_KLE;
  %input_file = @sample_input_navierstokes;
  input_file = @sample_input_explication;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % INITIALIZATION
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  maxIter = 20 ;
  refineTol = 0.5;   % refine elemets within refineTol of max 
  Etol = 1e-10;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % A bunch of parameters and settings
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  dim = 2;          % dimension of param space
  h = 1/32;         % initial mesh size (smooth = 1/8, discont = 1/64,NS=0.05)
  order = zeros(1,dim) + 1 ;    % initial order for Q 
  eorder = @(ord) 2 * ord ;     % how does eorder depend on order 
  qOrder = @(ord) max(1,4*ord); % control quadrature order for error calcs

  %order = [4 2 1 0 0 0 0 0 0 0] ;        % initial order for Q 
  %order = zeros(1,dim) + 1 ;        % initial order for Q 
  %eorder = @(ord) ord+1 ; % how does eorder depend on order 

  %------------------------------------------------------------------------------
  % refinement percenetages used to control when to split elements in parameter
  % space if local contribution accounts for over alpha percentage
  % alpha = 0   - always split every element
  % alpha = 1   - only split if all the error is in one element
  % alpha = 0.5 - split if half of the error is contained in that element
  % alpha = 1.5 - never split (identical to pure p-refinement)
  %------------------------------------------------------------------------------
  hAlpha = 0.25;
  NAlpha = 0.25;

  %------------------------------------------------------------------------------
  % parameter settings
  %------------------------------------------------------------------------------
  param=struct();
  param.inputFile = input_file;
  param.dim = dim;
  %param.min = [0.01 1];        % NS params
  %param.max = [0.1 3];         % NS params
  param.min = 0*ones(1,dim);
  param.max = 1*ones(1,dim);
  param.hInitial = h;
  param.numVars = 2; % total number of variables (forward+adjoint)

  %------------------------------------------------------------------------------
  %------------------------------------------------------------------------------
  % Control parameters - dictate how refinement is done etc
  %------------------------------------------------------------------------------
  %------------------------------------------------------------------------------
  % initialSplits - how many times to split all elements
  %         number of elem =  2^(dim*initialSplits) elements
  %      0 - dont split (start with one element)
  param.initialSplits = 2;
  %---------------------

  %---------------------
  % simultRefin - overrides all other refinement settings to always perfrom
  % h-refinement in physical space and p-refinement in parameter space
  param.simultRefine = false;
  %---------------------

  %---------------------
  % pRefine - flag to control p-refinement
  param.pRefine = false;
  param.forcePrefine = false;
  pIncrement = 1;

  % anisotropic - preferentially choose dimensions to refine based on coefficients
  %               of the error estimate
  param.anisotropic = false;
  %---------------------

  %---------------------
  % hRefine - flag to control h-refinement
  param.hRefine = true;
  param.uniformHrefine = false;
  %---------------------

  %---------------------
  % refinePhysical - flag to turn off dual refinement and only perform parameter
  %                  space refinement
  param.refinePhysical = true;
  param.onlyRefinePhysical = false;
  %---------------------

  %---------------------
  % smartSplit - controls if we use the neibhboring elements error indicators to
  %              inform split
  param.smartSplit = false;
  %---------------------

  %---------------------
  % resuse - flag to resuse lower order poly evaluation instead of resolving pde
  param.reusePoly = false;
  %---------------------

  %---------------------
  % evalTrue - whether to compute truth and true error or not
  param.evalTrue = true;
  %---------------------


  % maxNumElems - initial elem array size
  param.maxNumElems = max(100,...
    1+sum(2.^(dim*linspace(1,param.initialSplits,param.initialSplits)))...
    );
  %---------------------


  %------------------------------------------------------------------------------
  % setup initial param domain
  %------------------------------------------------------------------------------
  s = parameter('Legendre',param.min(1),param.max(1),0,0);
  if param.dim > 1
      for k = 2:param.dim
          s = [s parameter('Legendre',param.min(k),param.max(k),0,0)];
      end
  end

  %------------------------------------------------------------------------------
  % set up initial ACES structure
  %------------------------------------------------------------------------------
  param.ACESref(1).ACES = input_file('setup',s,param.hInitial);
  %param.ACESref(1).ACES = input_file('refine',[],param.ACESref(1).ACES);
  param.ACESref(1).ACES.checkSolve = 0;


  %------------------------------------------------------------------------------
  % Initial setup of parameter element structure
  %------------------------------------------------------------------------------
  newElem = struct(...
      'active',                 0,...
      'update',                 0,...
      'ACESref',                0,...
      's',                      struct([]),...
      'weight',                 0.,...
      'order',                  0*ones(1,dim),...
      'eorder',                 0*ones(1,dim),...
      'qOrder',                 0*ones(1,dim),...
      'Upoly',                  struct('coefficients',[]),...
      'Qpoly',                  struct('coefficients',[]),...
      'Epoly',                  struct('coefficients',[]),...
      'SEpoly',                 struct('coefficients',[]),...
      'errIndpoly',             struct('coefficients',[]),...
      'Eest',                   0.,...
      'Eh',                     0.,...
      'EN',                     0.,...
      'Etrue',                  0.,...
      'dof',                    zeros(param.numVars,1),...
      'parent',                 1,...
      'children',               [],...
      'dependentChildren',      [],...
      'solvetime',              0.,...
      'neighbors',              []...
      );
  param.elem = repmat(newElem,1,param.maxNumElems);

  param.elem(1) = struct(...
      'active',                 1,...
      'update',                 1,...
      'ACESref',                1,...
      's',                      s,...
      'weight',                 prod(param.max-param.min),...
      'order',                  order,...
      'eorder',                 eorder(order),...
      'qOrder',                 qOrder(eorder(order)),...
      'Upoly',                  struct('coefficients',[]),...
      'Qpoly',                  struct('coefficients',[]),...
      'Epoly',                  struct('coefficients',[]),...
      'SEpoly',                 struct('coefficients',[]),...
      'errIndpoly',             struct('coefficients',[]),...
      'Eest',                   0.,...
      'Eh',                     0.,...
      'EN',                     0.,...
      'Etrue',                  0.,...
      'dof',                    zeros(param.numVars,1),...
      'parent',                 1,...
      'children',               [],...
      'dependentChildren',      [],...
      'solvetime',              0.,...
      'neighbors',              []...
      );
  param.numElemsUsed = 1;

  % initialize error estimates and other running totals
  EtrueGlobal = 0;
  EestGlobal  = 0;
  EhGlobal    = 0;
  ENGlobal    = 0;
  Eff         = [];
  Eff2        = [];
  dofs        = zeros(param.numVars,1);
  cost        = zeros(param.numVars,1);
  nSolves     = [];
  omegaDofs   = [];
  refineSpace = [];

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Plot stuff
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %[p,w]=gaussian_quadrature(param.elem(1).s,20,0);
  %[XX,YY]=ndgrid(p{:});

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % construct initial piecewise parameter space 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for nRefine =1:param.initialSplits;
      for adaptElemId = find([param.elem(:).active])
          
          %create children
          children = splitElem(param.elem(adaptElemId));
          
          % set parent and mark children for update
          [children.parent]   = deal(adaptElemId);
          [children.children] = deal([]);
          [children.update]   = deal(1);
          
          % remove parent from active elements
          param.elem(adaptElemId).update    = 0;
          param.elem(adaptElemId).active    = 0;
          
          % add new elements to list
          param.elem(adaptElemId).children  = [param.numElemsUsed+1:param.numElemsUsed+length(children)];
          param.elem(param.numElemsUsed+1:param.numElemsUsed+length(children)) = children;
          param.numElemsUsed = param.numElemsUsed+length(children);
      end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Initial construction of approximation
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  fprintf('------------------------------------------------------------------\n')
  fprintf('              ITER 1 (initial settings)      \n')
  fprintf('------------------------------------------------------------------\n')
  nIter = 1;


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % update all flagged elements
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  [param,cost(:,nIter)] = updateElements_reuse(param);


  activeIndices = find([param.elem.active]);
  % calculate global quantities
  if param.evalTrue
    EtrueGlobal(nIter)  = sum([param.elem(activeIndices).Etrue]);
  else
    EtrueGlobal(nIter)  = 0;
  end
  EestGlobal(nIter)   = sum([param.elem(activeIndices).Eest]);
  EhGlobal(nIter)     = sum([param.elem(activeIndices).Eh]);
  ENGlobal(nIter)     = sum([param.elem(activeIndices).EN]);
  dofs(:,nIter)       = sum([param.elem.dof],2);
  cost(:,nIter)       = cost(:,nIter);

  nSolves   = [nSolves totalevals];
  dummyParam = {param.elem(activeIndices).order}';
  omegaDofs = [omegaDofs sum(prod([cell2mat(dummyParam)+1]'))];

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % PRINTING RESULTS
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % update effectivity
  if param.evalTrue
    Eff   = [Eff EestGlobal/EtrueGlobal];
    Eff2  = [Eff2 (EhGlobal+ENGlobal)/EtrueGlobal];
  else
    Eff   = [Eff 0];
    Eff2  = [Eff2 0];
  end
  nIter=1;

  fprintf('\n-------------------------------------------------------------------------------------\n')
  fprintf('nIter   EtrueGlobal   EestGlobal   EhGlobal       ENGlobal     eff.     eff2.   \n')
  fprintf(...
      ' %2i     %1.4e    %1.4e    %1.4e    %1.4e    %1.4f     %1.4f\n',...
      nIter,EtrueGlobal,EestGlobal,EhGlobal,ENGlobal,Eff(nIter),Eff2(nIter))
  fprintf('\n-------------------------------------------------------------------------------------\n\n')
end
restart=true;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADAPTIVITY LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for nIter = 2:maxIter
    fprintf('------------------------------------------------------------------\n')
    fprintf('              ITER %1i \n',nIter)
    fprintf('------------------------------------------------------------------\n')
    
    activeIndices = find([param.elem.active]);
    fprintf('   %0.1i active elements  \n',length(activeIndices))
    
    if (param.uniformHrefine)
      splitIndices = [];
    else
      splitIndices = ...
       find([param.elem(activeIndices).Eest] >= refineTol * max([param.elem(activeIndices).Eest])); 
    end


     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % which space to refine ??
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     if (param.refinePhysical && (ENGlobal(nIter-1) <= EhGlobal(nIter-1)) ) ...
       || param.simultRefine || param.onlyRefinePhysical
         fprintf('--> performing physical refinement \n')
         refineSpace = [refineSpace,0];
         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % Check - are there any proposals for splitting
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
         
         if (~isempty(splitIndices) && param.hRefine) && ~param.simultRefine 
             fprintf('----> Elements marked for refinement (%0.1i)\n',length(splitIndices))
             fprintf('       ')
             fprintf('%0.1i ',activeIndices(splitIndices))
             fprintf('\n')
             
             
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             % check if splitting elements is viable
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             for adaptElemId = activeIndices(splitIndices)
                 
                 if (param.elem(adaptElemId).dof(1) == 0)
                     % if this element is a dependent child we simply compute its sol
                     param.elem(adaptElemId).update = 1;
                 else
                     % construct new ACES structure for refined mesh
                     param.ACESref(length(param.ACESref)+1).ACES =...
                       input_file('refine',[],param.ACESref(param.elem(adaptElemId).ACESref).ACES,param.simultRefine);
                         
                     
                     % construct chlidren
                     children = splitElem(param.elem(adaptElemId));
                     [children.parent]             = deal(adaptElemId);
                     [children.children]           = deal([]);
                     [children.dependentChildren]  = deal([]);
                     [children.dof]                = deal(zeros(param.numVars,1));

                     % calculate what contribution would be from each element
                     for eID = 1:length(children)
                         % get quadrature poitns and weights
                         [pe,we] = gaussian_quadrature(children(eID).s,children(eID).eorder*2,1);
                         SEpts   = evaluate_expansion(children(eID).SEpoly,pe);
                         QofIpts = evaluate_expansion(children(eID).Qpoly,pe);
                         TEpts   = evaluate_expansion(children(eID).Epoly,pe);

                         % calculate error over this element
                         % --  since quadrature rule is based on local element basis we need to
                         %     account for weight
                         children(eID).Eh     = children(eID).weight *abs(SEpts)*we;
                         children(eID).Eest   = children(eID).weight *abs(TEpts)*we;
                         children(eID).EN     = children(eID).weight *abs(TEpts-SEpts)*we;

                     end % FOR - eID = 1:nElem

                     if param.evalTrue && (param.dim == 2)
                       for eID = 1:length(children)
                         [pe,we] = gaussian_quadrature(children(eID).s,children(eID).eorder*2,1);
                         Qtrue   = @(s1,s2) ...
                             arrayfun(@(t1,t2)input_file('integrate',[t1,t2],param.ACESref(children(eID).ACESref).ACES),s1,s2);
                         QTruepts  = Qtrue(pe(:,1),pe(:,2))';
                         children(eID).Etrue  = children(eID).weight * abs(QofIpts-QTruepts)*we;
                       end
                     end
       

                     % make sure piecewise contributions roughly add to total elem
                     if (abs(sum([children.Eh]) - param.elem(adaptElemId).Eh ) > 1e-6 );
                         sum([children.Eh])
                         param.elem(adaptElemId).Eh
                         %assert(false)
                     end

                     if param.smartSplit
                       % dummy ACESchildren to keep track of error indicators
                       % only really used if we enable smart split
                       ACESchildren = struct('ACES',[]);
                       for eID = 1:length(children)
                         [pe,we] = gaussian_quadrature(children(eID).s,children(eID).eorder*2,1);
                         % calculate avg error indicator over this region
                         ACESchildren(eID).ACES = param.ACESref(children(eID).ACESref).ACES;
                         ACESchildren(eID).ACES.mesh.expression = evaluate_expansion(children(eID).errIndpoly,pe)*we;
                       end
                     end
                     
                     fprintf('------> Checking for splitting Element %0.1i\n',adaptElemId)
                     dominateElem = find([children.Eh] > hAlpha * param.elem(adaptElemId).Eh);
                     if dominateElem
                         % split element
                         fprintf('------> Splitting Element %0.1i\n',adaptElemId)
                         
                         % add children to parent.children array
                         param.elem(adaptElemId).children = ...
                           [param.numElemsUsed+1:param.numElemsUsed+length(children)];
                         
                         %if it is a dominate element want to update
                         for eID = dominateElem

                             children(eID).ACESref = length(param.ACESref);
                             children(eID).update = 1;

                             if param.smartSplit
                               % check if error indicators are similar and recombine if so
                               [Y,I] = sort(ACESchildren(eID).ACES.mesh.expression);
                               num2refine = round(0.1*length(I));
                               refine = I(end-num2refine:end);
                             
                                 for neighbor = [children(eID).neighbors]'
                                     [Yn,In] = sort(ACESchildren(neighbor).ACES.mesh.expression);
                                     refineN = In(end-num2refine:end);
                                     notInCommon = length(setxor(refine,refineN));
                                     % recombine elements
                                     if (num2refine <= notInCommon)
                                         children(neighbor).active = 0;
                                         mins = min([children(eID).s(:).l],[children(neighbor).s(:).l]);
                                         maxs = max([children(eID).s(:).r],[children(neighbor).s(:).r]);
                                         children(eID).s = parameter('Legendre',mins(1),maxs(1),0,0);
                                         if dim > 1
                                             for k = 2:dim
                                                 children(eID).s = [children(eID).s parameter('Legendre',mins(k),maxs(k),0,0)];
                                             end
                                         end
                                     end % if num2refine <= notInCommon
                                 end % for neihbor = [children(eID).neihbors]'
                             end % if smartSplit
                             
                         end % for eID = dominateElem
                         
                         % -- check for overlapping elements after combining neighbors
                         for eID = setxor(find([children.active]),dominateElem)
                             for dominate = dominateElem
                                 % if this element is covered by dominate element remove it
                                 if (...
                                         sum([children(eID).s(:).r]<=[children(dominate).s(:).r]) == dim ...
                                         && ...
                                         sum([children(eID).s(:).l]>=[children(dominate).s(:).l]) == dim ...
                                         )
                                     children(eID).active = 0;
                                 end
                             end
                         end
                         
                         % -- if remaining active children that dont dominate it is  a depenedent child
                         for eID = setxor(find([children.active]),dominateElem)
                             param.elem(adaptElemId).dependentChildren = ...
                                 [param.elem(adaptElemId).dependentChildren param.numElemsUsed+eID];
                         end
                         
                         
                         % append children to param.elem vector
                         %param.elem = [param.elem children];
                         param.elem(param.numElemsUsed+1:param.numElemsUsed+length(children)) = children;
                         param.numElemsUsed = param.numElemsUsed+length(children);
                         
                         %mark element as inactive
                         param.elem(adaptElemId).update=0;
                         param.elem(adaptElemId).active=0;
                         
                     else % IF - find(piecewiseEh > hAlpha * param.elem(adaptElemId))
                         fprintf('------> Keeping Element %0.1i\n',adaptElemId)
                         param.elem(adaptElemId).ACESref = length(param.ACESref);
                         param.elem(adaptElemId).update = 1;
                         
                     end % IF - find(piecewiseEh > hAlpha * param.elem(adaptElemId))
                 end % IF param.elem(adaptElemId).dof(1) == 0
             end % FOR - adaptElemId = splitIndices
             
             
         else % IF - no elements proposed for splitting
             %select all elements to refine
             fprintf('----> No Proposals for splitting (refining all)\n')
             
             for eID = find([param.elem.active])
                 % refine unless the element is a dependent child, then just solve
                 if (param.elem(eID).dof(1) ~= 0)
                     param.ACESref(length(param.ACESref)+1).ACES = ...
                       input_file('refine',[],param.ACESref(param.elem(eID).ACESref).ACES,param.simultRefine);
                     param.elem(eID).ACESref=length(param.ACESref);
                 end % IF param.elem(adaptElemId).dof(1) ~= 0
                 param.elem(eID).update = 1;
             end % FOR - eID = find([param.elem.active])
             
         end % IF - ~isempty(splitIndices)

         
     end
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % PARAMETER REFINEMENT
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     if ( (param.pRefine || param.hRefine) ...
         && (ENGlobal(nIter-1) > EhGlobal(nIter-1)) ) ...
         || param.simultRefine || param.forcePrefine
        fprintf('--> performing parameter refinement \n')
         refineSpace = [refineSpace,1];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Check - are there any proposals for splitting
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % we have different ways to decide this
        %splitIndices = find([param.elem(:).EN] > NAlpha * ENGlobal(nIter-1));
        %splitIndices = find([param.elem(:).EN] >= 0.9* max([param.elem(:).EN]));
        %    [Y,I] = sort([param.elem(:).EN]);
        %    num2refine = ceil(0.1*length(find([param.elem.active]));
        %    splitIndices = I(end-num2refine+1:end);
        %splitIndices = find([param.elem(activeIndices).EN] >= ENGlobal(nIter-1)/10);
        %splitIndices = find([param.elem(activeIndices).Eest] >= EestGlobal(nIter-1)/10);
        %splitIndices = find([param.elem(activeIndices).Eest] >= refineTol/length(activeIndices));
        %splitIndices = find([param.elem(activeIndices).EN] >= refineTol.*[param.elem(activeIndices).weight]);
        %splitIndices = find([param.elem(activeIndices).EN] >= refineTol.*[param.elem(activeIndices).weight]);

        % this is a check in case we want to keep inital order but no elements
        % are proposed for splitting => split all
        refineAll = false;
        splitIndices;
        if ~param.pRefine && isempty(splitIndices)
          fprintf('----> No Proposals for splitting (splitting all)\n')
          fprintf('         param.pRefine == false         \n')
          refineAll = true;
          splitIndices = find([param.elem(activeIndices).active]);
        end

    
        if ~isempty(splitIndices) && param.hRefine
            fprintf('----> Elements marked for refinement\n')
            fprintf('       ')
            fprintf('%0.1i ',activeIndices(splitIndices))
            fprintf('\n')
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % check if splitting elements is viable
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for adaptElemId = activeIndices(splitIndices)
                
                if (param.elem(adaptElemId).dof(1) == 0)
                  % if this element is a dependent child we simply compute its sol
                  param.elem(adaptElemId).update = 1;
                else
                    

                  % construct chlidren
                  children = splitElem(param.elem(adaptElemId));
                  [children.parent]             = deal(adaptElemId);
                  [children.children]           = deal([]);
                  [children.dependentChildren]  = deal([]);
                  [children.dof]                = deal(zeros(param.numVars,1));
                    
                  % calculate what contribution would be from each element
                  for eID = 1:length(children)
                      % get quadrature poitns and weights
                      [pe,we] = gaussian_quadrature(children(eID).s,children(eID).eorder*2,1);
                      SEpts   = evaluate_expansion(children(eID).SEpoly,pe);
                      TEpts   = evaluate_expansion(children(eID).Epoly,pe);

                      % calculate error over this element
                      % --  since quadrature rule is based on local element basis we need to
                      %     account for weight
                      children(eID).Eh     = children(eID).weight *abs(SEpts)*we;
                      children(eID).Eest   = children(eID).weight *abs(TEpts)*we;
                      children(eID).EN     = children(eID).weight *abs(TEpts-SEpts)*we;
                  end % FOR - eID = 1:nElem


                  if param.evalTrue && (param.dim == 2)
                    for eID = 1:length(children)
                         [pe,we] = gaussian_quadrature(children(eID).s,children(eID).eorder*2,1);
                      QofIpts = evaluate_expansion(children(eID).Qpoly,pe);
                      Qtrue   = @(s1,s2) ...
                          arrayfun(@(t1,t2)input_file('integrate',[t1,t2],param.ACESref(children(eID).ACESref).ACES),s1,s2);
                      QTruepts  = Qtrue(pe(:,1),pe(:,2))';
                      children(eID).Etrue  = children(eID).weight * abs(QofIpts-QTruepts).^2*we;
                    end
                  end


                    if (abs(sum([children.EN]) - param.elem(adaptElemId).EN ) > 1e-6 );
                        sum([children.EN])
                        param.elem(adaptElemId).EN
                        %assert(false)
                    end
                    
                    
                    fprintf('------> Checking for splitting Element %0.1i\n',adaptElemId)
                    if refineAll
                      dominateElem = find([children.weight]);
                    else
                      dominateElem = find([children.EN] > NAlpha * param.elem(adaptElemId).EN);
                    end

                    % a check again in case we want to keep initial order but
                    % still refine something and not already refining all
                    if ~param.pRefine && isempty(dominateElem) 
                      [~,dominateElem] = max([children.EN]);
                    end

                    if dominateElem
                        % split element
                        fprintf('------> Splitting Element %0.1i\n',adaptElemId)
                        
                       % add children to parent.children array
                       param.elem(adaptElemId).children = ...
                         [param.numElemsUsed+1:param.numElemsUsed+length(children)];
                        
                        for eID = 1:length(children)
                            %if it is the dominate element want to update
                            if sum(eID == dominateElem)
                                % add new element to list and mark it to be updated
                                children(eID).update = 1;
                            else %IF eID = dominateElem
                                % otherwise we have a dependent child and reuse the parent sol
                                param.elem(adaptElemId).dependentChildren = ...
                                    [param.elem(adaptElemId).dependentChildren param.numElemsUsed+eID];
                            end % IF eID = dominateElem
                        end % FOR - eID = 1:nElem
                        
                        % append children to param.elem vector
                        %param.elem = [param.elem children];
                        param.elem(param.numElemsUsed+1:param.numElemsUsed+length(children)) = children;
                        param.numElemsUsed = param.numElemsUsed+length(children);
                        
                        %mark element as inactive
                        param.elem(adaptElemId).update=0;
                        param.elem(adaptElemId).active=0;
                        
                    else % IF - find(piecewiseEN > hAlpha * param.elem(adaptElemId))
                      fprintf('------> Keeping Element %0.1i\n',adaptElemId)

                      if param.anisotropic
                      %       - find out what most important dimension(s) are
                      %      (-)increase order in those directions
                        [C,I]=setdiff(param.elem(adaptElemId).Epoly.index_set',param.elem(adaptElemId).Qpoly.index_set','rows');
                        [sortedCoefficients,ids]=sort(abs(param.elem(adaptElemId).Epoly.coefficients(I)),2,'descend');
                        inc = C(ids(1),:) > param.elem(adaptElemId).order ;
                      else
                        inc = ones(1,dim);
                      end

                      inc = inc * pIncrement;

                      param.elem(adaptElemId).order   = param.elem(adaptElemId).order + inc;
                      param.elem(adaptElemId).eorder  = eorder(param.elem(adaptElemId).order) ;
                      param.elem(adaptElemId).update  = 1;
                    end % IF - find(piecewiseEN > hAlpha * param.elem(adaptElemId))
                    
                end % IF param.elem(adaptElemId).dof(1) == 0
            end % FOR - adaptElemId = splitIndices
            
            
        else % IF - no elements proposed for splitting
            %select all elements to refine
            fprintf('----> No Proposals for splitting (refining all)\n')
            
            for eID = find([param.elem.active])
                % refine unless the element is a dependent child, then just solve
                if param.elem(eID).dof(1) ~=0

                    if param.anisotropic
                    %       - find out what most important dimension(s) are
                    %      (-)increase order in those directions
                      [C,I]=setdiff(param.elem(eID).Epoly.index_set',param.elem(eID).Qpoly.index_set','rows');
                      [sortedCoefficients,ids]=sort(abs(param.elem(eID).Epoly.coefficients(I)),2,'descend');
                      inc = C(ids(1),:) > param.elem(eID).order ;
                    else
                      inc = ones(1,dim);
                    end

                    inc = inc * pIncrement;

                    param.elem(eID).order   = param.elem(eID).order + inc;
                    param.elem(eID).eorder  = eorder(param.elem(eID).order) ;

                end % IF param.elem(adaptElemId).dof(1) ~= 0
                param.elem(eID).update  = 1;
            end % FOR - eID = find([param.elem.active])
            
        end % IF - ~isempty(splitIndices)
        
        
     end %IF - ENGlobal > EhGlobal



     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % In case we aren't refining the space that needs it
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     if (~(param.refinePhysical || param.forcePrefine) && (ENGlobal(nIter-1) <= EhGlobal(nIter-1)))
       assert(false,'ERROR: space to be refined is not marked for refinement')
     end
       
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % initialize zero error values because everything will be updated
    EtrueGlobal = [EtrueGlobal 0];
    EestGlobal = [EestGlobal 0];
    EhGlobal = [EhGlobal 0];
    ENGlobal = [ENGlobal 0];
    dofs = [dofs zeros(param.numVars,1)];
    cost = [cost zeros(param.numVars,1)];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update all flagged elements
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [param,cost(:,nIter)] = updateElements_reuse(param);
    cost(:,nIter)  = cost(:,nIter) + cost(:,nIter-1);
    
    activeIndices = find([param.elem.active]);
    % calculate global quantities
    if param.evalTrue
      EtrueGlobal(nIter)  = sum([param.elem(activeIndices).Etrue]);
    else 
      EtrueGlobal(nIter)  = 0;
    end
    EestGlobal(nIter)   = sum([param.elem(activeIndices).Eest]);
    EhGlobal(nIter)     = sum([param.elem(activeIndices).Eh]);
    ENGlobal(nIter)     = sum([param.elem(activeIndices).EN]);
    dofs(:,nIter)       = sum([param.elem.dof],2);

    nSolves      = [nSolves totalevals];
    dummyParam = {param.elem(activeIndices).order}';
    omegaDofs = [omegaDofs sum(prod([cell2mat(dummyParam)+1]'))];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PRINTING INTERMEDIATE RESULTS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update effectivity
    if param.evalTrue
      Eff(nIter)   = EestGlobal(nIter)/EtrueGlobal(nIter);
      Eff2(nIter)  = (EhGlobal(nIter)+ENGlobal(nIter))/EtrueGlobal(nIter);
    else
      Eff(nIter)   = 0;
      Eff2(nIter)  = 0;
    end
    
    fprintf('\n-------------------------------------------------------------------------------------\n')
    fprintf('nIter   EtrueGlobal   EestGlobal    EhGlobal      ENGlobal      eff.       eff2.   \n')
    fprintf(...
        ' %2i     %1.4e    %1.4e    %1.4e    %1.4e    %1.4f     %1.4f\n',...
        nIter,EtrueGlobal(nIter),EestGlobal(nIter),EhGlobal(nIter),ENGlobal(nIter),Eff(nIter),Eff2(nIter))
    fprintf('\n-------------------------------------------------------------------------------------\n\n')
    
    if EestGlobal(nIter) < Etol
        break
    end
end % ADAPTIVE LOOP

% save last iteration number in case we stopped early
maxIter =nIter;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRINTING FINAL RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n-------------------------------------------------------------------------------------\n')
fprintf('nIter   EtrueGlobal   EestGlobal   EhGlobal       ENGlobal     eff.     eff2.   \n')
for nIter=1:maxIter
    fprintf(...
        ' %2i     %1.4e    %1.4e    %1.4e    %1.4e    %1.4f     %1.4f\n',...
        nIter,EtrueGlobal(nIter),EestGlobal(nIter),EhGlobal(nIter),ENGlobal(nIter),Eff(nIter),Eff2(nIter))
end
fprintf('\n-------------------------------------------------------------------------------------\n\n')

%cumCost = [cost(1:end/2,1), cost(1:end/2,2:end)+cost(end/2+1:end,1:end-1)];
%results=[sum(dofs)',sum(cumCost,1)',EestGlobal',EhGlobal',ENGlobal'];
%results=[,EestGlobal',EhGlobal',ENGlobal'];
results = [sum(dofs)',cumsum(sum(dofs)'),EestGlobal',EhGlobal',ENGlobal'];

fid=fopen('results.dat','w');
fprintf(fid,'%% dofs  cost  Eest  Eh  EN\n');
fprintf(fid,'%0.5e %0.5e %0.5e %0.5e %0.5e\n',results');
fclose(fid);



