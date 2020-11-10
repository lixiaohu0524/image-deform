clear all
V =  [  0.41089       3.2532;
       2.8251       2.7639;
      0.14673      0.96887;
       2.0997       1.4961];


F=delaunay(V)
C=[1,2;1,1.5];
bc=[0,1;1,0]
n = size(V,1);
m = size(bc,2);
L = cotmatrix(V,F)
M = massmatrix(V,F,'voronoi')
low = 0;
up = 1;
type = 'conic';
b=[3,4]


param = [];
Qi = L*(M\L)

Q = sparse(m*n,m*n)

for ii = 1:m
    d = (ii - 1)*n + 1;
    Q(d:(d + n-1), d:(d + n-1)) = Qi;
end
Q
PA = repmat(speye(n,n),1,m)
Pb = ones(n,1)

BCAi = speye(n,n)
BCAi = BCAi(b,:)
BCA = sparse(m*size(BCAi,1),m*size(BCAi,2))

for ii = 1:m
    di = (ii - 1)*size(BCAi,1) + 1;
    dj = (ii - 1)*size(BCAi,2) + 1;
    BCA(di:(di + size(BCAi,1)-1), dj:(dj + size(BCAi,2)-1)) = BCAi;
end

BCb=bc(:)

ux = up.*ones(m*n,1)
lx = low.*ones(m*n,1)


pabca=[PA;BCA]
pabca=[Pb;BCb]
W = quadprog(Q,zeros(n*m,1),[],[],[PA;BCA],[Pb;BCb],lx,ux,[],param)



W = reshape(W,n,m);

else
    error( [ ...
        'Enforcing partition of unity only support in conjunction with ' ...
        'type=''quad''']);
    end
    
    else
        
        if(strcmp(type,'quad'))
            
            Q = L*(M\L);
            
            ux = up.*ones(n,1);
            lx = low.*ones(n,1);
        elseif(strcmp(type,'least-squares'))
            
            I = speye(n);
            Z = sparse(n,n);
            Q = [Z,Z;Z,I];
            F = sqrt(M)\L;
            c = zeros(n,1);
            B = [F,-I];
            ux = [up.*ones(n,1) ;  Inf*ones(n,1)];
            lx = [low.*ones(n,1); -Inf*ones(n,1)];
        elseif(strcmp(type,'conic'))
            
            F = sqrt(M)\L;
            prob.c = [zeros(2*n,1); 1];
            I = speye(n);
            prob.a = [F,-I,zeros(n,1)];
            prob.blc = zeros(n,1);
            prob.buc = zeros(n,1);
            prob.bux = [ up.*ones(n,1);  Inf*ones(n,1);  Inf];
            prob.blx = [ low.*ones(n,1); -Inf*ones(n,1); -Inf];
            prob.cones = cell(1,1);
            prob.cones{1}.type = 'MSK_CT_QUAD';
            t_index = 2*n +1;
            z_indices = (n+1):(2*n);
            prob.cones{1}.sub = [t_index z_indices];
        else
            error('Bad type');
        end
        
        
        m = size(bc,2);
        
        W = zeros(n,m);
        tic;
        
        for i = 1:m
            if(strcmp(type,'quad'))
                
                Aeq = speye(n,n);
                Aeq = Aeq(b,:);
                if(mosek_exists)
                    fprintf('Quadratic optimization using mosek...\n');
                else
                    fprintf('Quadratic optimization using matlab...\n');
                end
                fprintf( [ ...
                    '  minimize:     x''LM\\Lx\n' ...
                    'subject to: %g <= x <= %g\n' ], ...
                    low,up);
                
                [x,fval,err] = quadprog(Q,zeros(n,1),[],[],Aeq,bc(:,i),lx,ux,[],param);
                if(err ~= 1)
                    fprintf([...
                        '----------------------------------------------------------\n' ...
                        'ERROR ('  num2str(err) ',' num2str(fval) '):' ...
                        ' solution may be inaccurate...\n' ...
                        '----------------------------------------------------------\n' ...
                        ]);
                end
            elseif(strcmp(type,'least-squares'))
                
                lx(b) = bc(:,i);
                ux(b) = bc(:,i);
                fprintf('Quadratic optimization using mosek...\n');
                fprintf([ ...
                    '  minimize:       z''z\n' ...
                    '  subject to: sqrt(M)\\Lx - z = 0\n' ...
                    '  and          %g <= x <= %g\n'] , ...
                    low,up);
                x = quadprog(Q,zeros(2*n,1),[],[],B,c,lx,ux,[],param);
            elseif(strcmp(type,'conic'))
                prob.bux(b) = bc(:,i);
                prob.blx(b) = bc(:,i);
                fprintf('Conic optimization using mosek...\n');
                fprintf([ ...
                    '  minimize:         t\n' ...
                    '  subject to: sqrt(M)\\Lx - z = 0,\n' ...
                    '             t >= sqrt(z''z),\n' ...
                    '               %f <= x <= %f\n'], ...
                    low,up);
                [r,res]=mosekopt('minimize echo(0)',prob,param);
                % check for mosek error
                if(r == 4006)
                    warning(['MOSEK ERROR. rcode: ' ...
                        num2str(res.rcode) ' ' ...
                        res.rcodestr ' ' ...
                        res.rmsg ...
                        'The solution is probably OK, but ' ...
                        'to make this error go away, increase: ' ...
                        'MSK_DPAR_INTPNT_CO_TOL_REL_GAP' ...
                        n]);
                elseif(r ~= 0)
                    error(['FATAL MOSEK ERROR. rcode: ' ...
                        num2str(res.rcode) ' ' ...
                        res.rcodestr ' ' ...
                        res.rmsg]);
                end
                
                x = res.sol.itr.xx;
            end
            
            W(:,i) = x(1:n);
            fprintf('Lap time: %gs\n',toc);
        end
        t = toc;
        fprintf('Total elapsed time: %gs\n',t);
        fprintf('Average time per handle: %gs\n',t/m);
        
