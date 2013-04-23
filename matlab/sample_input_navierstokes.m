function out = sample_input_navierstokes(dowhat,params,varargin)

switch dowhat
    
    case 'setup'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Define the geometry - see useGeometry.m and ../geometries
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        hsize = 0.03;
        if length(varargin) > 0;
          hsize=varargin{1};
        end
        ACES = useGeometry('fat_box_with_hole',hsize);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Define the physics
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Define the functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        ACES.function.diff = inline('p1');        
        ACES.function.tau_pspg = 0;
        ACES.function.tau_supg = 0;
        ACES.function.tau_lsic = 0;
        ACES.function.source1 = 0;
        ACES.function.source2 = 0;
        ACES.function.source3 = 0;
        ACES.function.bx = inline('u1');
        ACES.function.by = inline('u2');
        ACES.function.divb = inline('u1x+u2y');
        
        ACES.function.QofI1 = inline('10/pi*exp(-10*(x-1).^2 - 10*(y-0).^2)');
        ACES.function.QofI2 = 0;
        ACES.function.QofI3 = 0;
        %ACES.function.QofI3 = inline('10/pi*exp(-10*(x+1).^2 - 10*(y-0).^2)');
        ACES.function.BQofI1 = 0; % circumference of cylinder is pi/2
        ACES.function.BQofI2 = 0;
        
        ACES.function.errorrep1 = inline('source1*u4 - diff*u1x*u4x - diff*u1y*u4y + u3*u4x - bx*u1x*u4 - by*u1y*u4');
        ACES.function.errorrep2 = inline('source2*u5 - diff*u2x*u5x - diff*u2y*u5y + u3*u5y - bx*u2x*u5 - by*u2y*u5');
        ACES.function.errorrep3 = inline('source3*u6 - (u1x+u2y)*u6');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Use template to define physics
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %addpath ../work_stabilized/
        
        type = 'PSPG-LSIC-SUPG';
        
        ACES = addStabTermsStokesWithConv2(ACES,1:3,type,{'diff','bx','by','divb',{'u1x','u1y','u2x','u2y'},'1','1','source1 + bx*u1x+by*u1y','source2 + bx*u2x+by*u2y','source3'});

        ACES = addStabTermsStokesWithConv2(ACES,4:6,type,{'diff','-bx','-by','-divb',{'u1x-divb','u2x','u1y','u2y-divb'},'-1','-1','QofI1','QofI2','QofI3'});
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Define boundary conditions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        ACES.physics(1,1).boundary(4).value = inline('p2*3/32*(2-y)*(y+2)');
        ACES.physics(1,1).boundary(1).type = 'N';
        ACES.physics(2,2).boundary(1).type = 'N';
        
        for j = 1:length(ACES.mesh(1).boundary)
            ACES.physics(3,3).boundary(j).type  = 'none';
        end

        ACES.physics(4,4).boundary(1).type = 'R';
        ACES.physics(4,4).boundary(1).robincoeff = inline('(bx*nx+by*ny)');
        
        ACES.physics(5,5).boundary(1).type = 'R';
        ACES.physics(5,5).boundary(1).robincoeff = inline('(bx*nx+by*ny)');
        
        ACES.physics(5,4).boundary(5).value = inline('BQofI1');
        ACES.physics(5,5).boundary(5).value = inline('BQofI2');
        
        for j = 1:length(ACES.mesh(1).boundary)
            ACES.physics(6,6).boundary(j).type  = 'none';
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Initialize the problem - see documentation for options
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ACES = initialize(ACES,'method',11*ones(1,6),'degree',[2 2 1 3 3 2],'qrule',3*ones(1,6));
        ACES = getDefaultSolverSettings(ACES);
        [ACES,~] = integrate(ACES,inline('1'));
        
        out = ACES;
        
    case 'solve'
        
        ACES = varargin{1};
        
        for k = 1:length(params)
            ACES.function.(strcat('p',num2str(k))) = params(k);
        end
                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SOLVE THE PROBLEM
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        ACES = solve(ACES,'whichvars',1:3,'solver','NL');
        ACES = solve(ACES,'whichvars',4:6);
        
        [ACES,QofI] = integrate(ACES,inline('QofI1*u1+QofI2*u2+QofI3*u3'));
        [ACES,BQofI1] = integrate(ACES,inline('-BQofI1*(diff*(u1x*nx+u1y*ny)-u3)'),'type','boundary','bnum',5);
        [ACES,BQofI2] = integrate(ACES,inline('-BQofI2*(diff*(u2x*nx+u2y*ny)-u3)'),'type','boundary','bnum',5);
                        
        [ACES,err11] = integrate(ACES,inline('errorrep1'),'equation',1);
        errorInd=abs(ACES.mesh.expression);
        [ACES,err21] = integrate(ACES,inline('errorrep2'),'equation',2);
        errorInd=errorInd+abs(ACES.mesh.expression);
        [ACES,err31] = integrate(ACES,inline('errorrep3'),'equation',3);
        errorInd=errorInd+abs(ACES.mesh.expression);
        
        total_error_est = err11+err21+err31;

        fwdsol = [ACES.currentsol(1).sol;ACES.currentsol(2).sol;ACES.currentsol(3).sol];
        adjsol = [ACES.currentsol(4).sol;ACES.currentsol(5).sol;ACES.currentsol(6).sol];
        
%         FEMplot(ACES,'sqrt(u1^2+u2^2)')
%         FEMplot(ACES,'sqrt(u4^2+u5^2)')
        
        out = [ACES.settings.NLiter;fwdsol;adjsol;errorInd;0;total_error_est;QofI+BQofI1+BQofI2];
        
    case 'error_estimate'
        
        ACES = varargin{1};
        Upoly = varargin{2};
        
        sol = evaluate_expansion(Upoly,params);
        
        for k = 1:length(params)
            ACES.function.(strcat('p',num2str(k))) = params(k);
        end
        
        prog = 0;
        for k = 1:6
            ACES.currentsol(k).sol = sol(prog+1:prog+ACES.settings.ndof(k));
            ACES.solution(k).timeinterval(1).sol = ACES.currentsol(k).sol;
            prog = prog + ACES.settings.ndof(k);
        end
        
        [ACES,err11] = integrate(ACES,inline('errorrep1'),'equation',1);
        [ACES,err21] = integrate(ACES,inline('errorrep2'),'equation',2);
        [ACES,err31] = integrate(ACES,inline('errorrep3'),'equation',3);
        
        total_error_est = err11+err21+err31;
        
        out = total_error_est;

      case 'refine'

        ACES = varargin{1};

        for k = 1:length(params)
            ACES.function.(strcat('p',num2str(k))) = params(k);
        end


        ACES = mark(ACES,'markstrategy','uniform');
        ACES = hadapt(ACES);
        ACES = initialize(ACES,'method',11*ones(1,6),'degree',[2 2 1 3 3 2],'qrule',3*ones(1,6));
        ACES = getDefaultSolverSettings(ACES);
        [ACES,QofI] = integrate(ACES,inline('1'));
        
        out = ACES;

    case 'mean_fields'
        
        ACES = varargin{1};
        Upoly = varargin{2};
        
        sol = Upoly.coefficients(:,1);
        
        for k = 1:length(params)
            ACES.function.(strcat('p',num2str(k))) = params(k);
        end
        
        prog = 0;
        for k = 1:6
            ACES.currentsol(k).sol = sol(prog+1:prog+ACES.settings.ndof(k));
            ACES.solution(k).timeinterval(1).sol = ACES.currentsol(k).sol;
            prog = prog + ACES.settings.ndof(k);
        end
        
        
        out = ACES;
        
end
