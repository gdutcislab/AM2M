function mop = testmop( testname, dimension )
%Get test multi-objective problems from a given name.
%   The method get testing or benchmark problems for Multi-Objective
%   Optimization. The implemented problems included ZDT, OKA, KNO.
%   User get the corresponding test problem, which is an instance of class
%   mop, by passing the problem name and optional dimension parameters.

mop=struct('name',[],'od',[],'pd',[],'domain',[],'func',[]);
path('gtest',path);
% switch lower(testname)
switch testname
    case 'MaOP1'
        mop=MaOP1(mop);
    case 'MaOP2'
        mop=MaOP2(mop);
    case 'MaOP3'
        mop=MaOP3(mop);
    case 'MaOP4'
        mop=MaOP4(mop);
    case 'MaOP5'
        mop=MaOP5(mop);
    case 'MaOP6'
        mop=MaOP6(mop); 
    case 'MaOP7'
        mop=MaOP7(mop);
    case 'MaOP8'
        mop=MaOP8(mop);
    case 'MaOP9'
        mop=MaOP9(mop);
    case 'MaOP10'
        mop=MaOP10(mop);
    case 'MaOP11'
        mop=MaOP11(mop);
    case 'MaOP12'
        mop=MaOP12(mop);
    case 'MaOP13'
        mop=MaOP13(mop);
    case 'MaOP14'
        mop=MaOP14(mop);
    case 'MaOP15'
        mop=MaOP15(mop);
    case 'P16'
        mop=P16(mop);
    case 'WFG1'
        mop=wfg1(mop);
    case 'WFG2'
        mop=wfg2(mop);
    case 'WFG3'
        mop=wfg3(mop);
    case 'WFG4'
        mop=wfg4(mop);
    case 'WFG5'
        mop=wfg5(mop);
    case 'WFG6'
        mop=wfg6(mop);
    case 'WFG7'
        mop=wfg7(mop);
    case 'WFG8'
        mop=wfg8(mop);
    case 'WFG9'
        mop=wfg9(mop);
    case 'DTLZ7'
        mop=DTLZ7(mop);
    otherwise 
        error('Undefined test problem name');                
end 
end

%KNO1 function generator
function p=kno1(p)
 p.name='KNO1';
 p.od = 2;
 p.pd = 2;
 p.domain= [0 3;0 3];
 p.func = @evaluate;
 
    %KNO1 evaluation function.
    function y = evaluate(x)
      y=zeros(2,1);
	  c = x(1)+x(2);
	  f = 9-(3*sin(2.5*c^0.5) + 3*sin(4*c) + 5 *sin(2*c+2));
	  g = (pi/2.0)*(x(1)-x(2)+3.0)/6.0;
	  y(1)= 20-(f*cos(g));
	  y(2)= 20-(f*sin(g)); 
    end
end

%ZDT1 function generator
function p=zdt1(p,dim)
 p.name='ZDT1';
 p.pd=dim;
 p.od=2;
 p.domain=[zeros(dim,1) ones(dim,1)];
 p.func=@evaluate;
 
    %KNO1 evaluation function.
    function y=evaluate(x)
        y=zeros(2,1);
        y(1) = x(1);
    	su = sum(x)-x(1);    
		g = 1 + 9 * su / (dim - 1);
		y(2) =g*(1 - sqrt(x(1) / g));
    end
end



 