function [param,cost] = updateElements_reuse(param)
  activeIndices = find([param.elem.active]);
  updateIndices = activeIndices(find([param.elem(activeIndices).update]));
  num2update = length(updateIndices);

  fprintf('--> Updating flagged elements (%0.1i) \n',num2update)
  

  cost = 0;

  reusePoly = param.reusePoly;
  inputFile = param.inputFile;
  evalTrue = param.evalTrue;
  dim = param.dim;
  elem = param.elem(updateIndices);
  [ACESref] = [param.ACESref([param.elem(updateIndices).ACESref]).ACES];

  parfor eID = 1:length(updateIndices)
      tic
      currentElem = elem(eID);
      acesID = currentElem.ACESref;
      ACES = ACESref(eID);

      fprintf('----> element %0.1i (%0.1i)\n',updateIndices(eID),eID)
      fprintf('   order =  ',eID)
      fprintf('%0.1i ',currentElem.order)
      fprintf('\n',eID)
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % compute values for reuse
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      sizeOfACES = sum(ACES.settings.ndof);
      sizeOfACES = sizeOfACES +...
        size(ACES.mesh.t,1) + 2;
      sizeOfUpoly = size(currentElem.Upoly.coefficients,1);

      if reusePoly && (sizeOfACES == sizeOfUpoly)
        p=gaussian_quadrature(currentElem.s,currentElem.order+1);
        Upts =  evaluate_expansion(currentElem.Upoly,p);
        Epts = zeros(1,size(p,1));
        for k = 1:size(p,1)
            Epts(1,k) =...
              inputFile('error_estimate',p(k,:),ACES,currentElem.Upoly);
        end
    %         Epts =  evaluate_expansion(Epoly,p);
        
        ACES.newpts = p;
        ACES.Upts = Upts;
        ACES.Epts = Epts;
        ACES.checkSolve = 1;
      else
        ACES.checkSolve = 0;
      end

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Construct the PCE approximation of the response
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      u = @(t) inputFile('solve',t,ACES);
      

      [Upoly,errx] = ...
          pseudospectral(u,currentElem.s,currentElem.order,'Verbose',1);%,...
      
      NLiter = max(1,Upoly.coefficients(1,1) );
      coefficients = Upoly.coefficients(2:end,:);
      Upoly = rmfield(Upoly,'coefficients');
      
      % Extract the polynomial of the QofI
      currentElem.Qpoly = Upoly;
      currentElem.Qpoly.coefficients = coefficients(end,:);
      
      % Extract the polynomial of the spatial error
      currentElem.SEpoly = Upoly;
      currentElem.SEpoly.coefficients = coefficients(end-1,:);

      % Extract the polynomial of the error indicator
      startIdx = 1+sum(ACES.settings.ndof);
      stopIdx = startIdx-1+length(ACES.mesh.t);
      currentElem.errIndpoly = Upoly;
      currentElem.errIndpoly.coefficients = coefficients(startIdx:stopIdx,:);
      ACES.mesh.expression =...
        coefficients(startIdx:stopIdx,1);
      
      Upoly.coefficients = coefficients(1:startIdx-1,:);
      currentElem.Upoly = Upoly;
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Construct the PCE approximation of the error
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      e = @(t) inputFile('error_estimate',t,ACES,currentElem.Upoly);
      [currentElem.Epoly,errx] = ...
          pseudospectral(e,currentElem.s,currentElem.eorder,'Verbose',1);%,...
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % calculate errors
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      [pe,we] = ...
        gaussian_quadrature(currentElem.s,currentElem.qOrder,1);
        %gaussian_quadrature(currentElem.s,currentElem.eorder+2,1);
      
      QofIpts   = evaluate_expansion(currentElem.Qpoly,pe);
      SEpts     = evaluate_expansion(currentElem.SEpoly,pe);
      TEpts     = evaluate_expansion(currentElem.Epoly,pe);

      % calculate error over this element
      currentElem.Eest    = currentElem.weight * abs(TEpts)*we;
      currentElem.Eh      = currentElem.weight * abs(SEpts)*we;
      currentElem.EN      = currentElem.weight * abs(TEpts-SEpts)*we;

      % Eest2 is a higher order estimate of the error - calculated using
      % quadratrue on the true value of response not PC
      if isnan(currentElem.Eest)
        fprintf('\n\n')
        fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
        fprintf('WARNING: error estimate evaluated to NaN\n');
        fprintf('         testing with higher order quadrature\n');
        TE2pts = zeros(1,length(pe));
        for h=1:length(pe)
          TE2pts(h)    = inputFile('error_estimate',pe(h,:),ACES,currentElem.Upoly);
        end
        Eest2 = currentElem.weight * abs(TE2pts)*we;
        fprintf('Elem: %1.5f  Eest1=%1.5e Eest2=%1.5e\n',updateIndices(eID),currentElem.Eest,Eest2)
        fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
        fprintf(' \n\n')
        currentElem.Eest = Eest2;
        currentElem.EN   = currentElem.weight * abs(TE2pts-SEpts)*we;
      end

      if evalTrue && (dim == 2)
        %     Qtrue =@(s1,s2)arrayfun(@(t1,t2)inputFile('integrate',[t1,t2],currentElem.ACES),s1,s2);
        Qtrue =@(s1,s2)arrayfun(...
          @(t1,t2)inputFile('integrate',[t1,t2],ACES),s1,s2);
        QTruepts  = Qtrue(pe(:,1),pe(:,2))';
        currentElem.Etrue   = currentElem.weight * abs(QofIpts-QTruepts)*we;
      end
      
      
      % add contributions to total cost
      % if nonlinear problem 
      % cost = (numIterations * dofs_primal) + dofs_adjoint
      currentElem.dof   = ACES.settings.ndof' * prod(currentElem.order+1);
      cost = cost + [(NLiter * currentElem.dof(1:end/2)); ...
        currentElem.dof(end/2+1:end) ]  ;
      
      currentElem.update = 0;
      
      currentElem.solvetime = toc;
      ACESref(eID) = ACES;
      elem(eID) = currentElem;
      %fprintf('----> %1.5f \n',currentElem.solvetime)
  end %FOR - eID = updateIndices
  param.elem(updateIndices) = elem;
  ACESref = mat2cell(ACESref,[1],ones(1,numel(ACESref)));
  [param.ACESref([param.elem(updateIndices).ACESref]).ACES] = deal(ACESref{:}) ;

  for eID = updateIndices
    % remove current updated element from dependent list of parent (if there)
    myParent = param.elem(eID).parent;
    [param.elem(myParent).dependentChildren(find([param.elem(myParent).dependentChildren]==eID))]=[];
    if isempty([param.elem(myParent).dependentChildren]) && (myParent ~= eID)
      param.elem(myParent).dof = zeros(param.numVars,1);
    end
  end
end
