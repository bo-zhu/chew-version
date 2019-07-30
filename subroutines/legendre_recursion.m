%% Legendre function recursion on degree l .
%% By ZHU Bo at Nanjing University ( bzhu@nju.edu.cn )

function p3 = legendre_recursion(fun,x,p1,p2,l,m)
% fun : 'legendre' or 'sch'
% x : argument
% p1,p2 : using p1 and p2 to compute p3
% p1 = p^m_l-1, p2 = p^m_l, and p3 = p^m_l+1 .

switch fun
case 'legendre'
	p3 = ( (2*l+1)*x*p2 - (l+m)*p1 ) / (l-m+1);
case 'sch'
	p3 = ( (2*l+1)*x*p2 - sqrt(l^2-m^2)*p1 ) / sqrt( (l+1)^2-m^2 );
end

end
