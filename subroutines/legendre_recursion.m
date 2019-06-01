%% recursion of Legendre functions
%% By ZHU Bo at Nanjing University ( bzhu@nju.edu.cn )

function p3 = legendre_recursion(fun,varying,x,p1,p2,l,m)
% fun : 'legendre' or 'sch'
% varying : 'order' or 'degree'
% x : argument
% p1,p2 : using p1 and p2 to compute p3
% l : degree of p2, so l+1 is the degree of p3 in the case of varying degree.  
% m : order of p2, so m+1 is the order of p3 in the case of varying order.

switch fun
case {'legendre'}
	switch varying
	case 'order'
		p3 = -2*m*x/sqrt(1-x^2)*p2-(l+m)*(l-m+1)*p1;
	case 'degree'
		p3 = ( (2*l+1)*x*p2-(l+m)*p1 )/(l-m+1);
	endswitch
case {'sch'}
	switch varying
	case 'order'
		p3 = 2*m*x/sqrt( (1-x^2)*(l+m+1)*(l-m) )*p2-sqrt( (l+m)*(l-m+1)/(l+m+1)/(l-m) )*p1;
	case 'degree'
		p3 = ( (2*l+1)*x*p2-sqrt(l^2-m^2)*p1 )/sqrt( (l+1)^2-m^2 );
	endswitch
endswitch

endfunction
