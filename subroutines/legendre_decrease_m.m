%% Legendre function recursion on decreasing m .
%% By ZHU Bo at Nanjing University ( bzhu@nju.edu.cn )

function p1 = legendre_decrease_m(fun,x,p2,p3,l,m)
% fun : 'legendre' or 'sch'
% x : argument
% p2,p3 : using p2 and p3 to compute p1
% p1 = p^m-1_l, p2 = p^m_l, and p3 = p^m+1_l .

if abs(1-abs(x))<1e-10

  printf('x = (+/-)1 is the singularity of the recurrence formula !');
  quit;

else

  switch fun
  case 'legendre'
    p1 = -1/(l+m)/(l-m+1) * (2*m*x/sqrt(1-x^2)*p2 + p3);
  case 'sch'
    p1 = 1/sqrt( (l+m)*(l-m+1) ) * ( 2*m*x/sqrt(1-x^2)*p2 - sqrt( (l+m+1)*(l-m) )*p3 );
  end

end

end
