
Comment construct one element model of stochastic modified
cubic reaction-diffusion equation with multiplicative noise,
to effects quadratic in the noise amplitude, and seeks the
normal form where the evolution involves no convolutions.
Also, transform the quadratic noise in the evolution.

Tony Roberts, 27 June 2005, Dec 2011;

% improve printing
linelength 60$
on div; off allfac; on revpri; factor a,s;
procedure beta(m); (m^2-1); % decay rate of linear modes
operator linv; linear linv;
let linv(sin(~n*x),x)=>sin(n*x)/beta(n);

% parametrise by evolving s 
depend s,t;
let df(s,t)=>g;

% linear approximation
u:=s*sin(x);
g:=0;

% iterate towards a solution
let { a^3=>0, gam=>0 };
repeat begin
    deq:=trigsimp(-df(u,t)+df(u,x,2)+u-a*u^3,combine);
    g:=g+(gd:=(deq where {sin(x)=>1,sin(~n*x)=>0}));
    u:=u+linv(deq-gd*sin(x),x);
    showtime;
end until deq=0;



% noise: tt labels the fast time of stochastic fluctuations
% let phi(n,{m1,...}) denote convolutions with exp(-beta(m1)t)...
% so df(phi(n,m.p),t) = -beta(m)phi(n,m.p) + phi(n,p)
factor sig;
depend tt,t;
depend x,xx; 
operator phi; depend phi,tt,xx;
operator secular; linear secular;
procedure gunga(m,p);
if p={} then gunga_error(m,p) else -beta(first(p))*phi(m,p)+phi(m,rest(p));
procedure gungb(n,p);
if p={} then 0 else (gungb(n,rest(p))-phi(n,p))/beta(first(p));
let { df(phi(~m,~p),t)=>df(phi(~m,~p),tt)
    , df(phi(~m,~p),tt)=>gunga(m,p)
    , linv(sin(~m*x),xx)=>sin(m*x)/beta(m)
    , linv(sin(~m*x)*phi(~n,~p),xx)=>phi(n,m.p)*sin(m*x)
    , linv(sin(x)*phi(~n,~p),xx)=>gungb(n,p)*sin(x)
    , secular(sin(~m*x),xx)=>0
    , secular(sin(~m*x)*~aa,xx)=>0
    , secular(sin(x),xx)=>1
    , secular(sin(x)*phi(~n,~p),xx)=>
      phi(n,{})/(for each m in p product beta(m))
    };

% truncate the spatial structure of noise, may go up to 12 in 2mins
noise:=phi(0,{});%for n:=1:3 sum phi(n,{})*sin(n*x);

% now derive effects linear in noise
let sig^2=>0;
it:=0$
repeat begin
    deq:=trigsimp(-df(u,t)+df(u,x,2)+u-a*u^3+sig*noise*u,combine);
    g:=g+(gd:=secular(deq,xx));
    u:=u+linv(deq-gd*sin(x),xx);
    showtime;
end until ((it:=it+1)>9)or(deq=0);

write "Slow manifold to error O(a^3,sig^2):";
noise:=noise;
dsdt:=g;
uslow:=u;

end;%%%%%%%%%%%%%%%%%%%%%%%%%%%

% zz(a,p) denotes multiple convolution of (nonlinear) a 
% zz(a,m.p) = exp[-beta(m)t]*zz(a,p)  and  zz(a,{})=a
operator zz; depend zz,tt,xx; 
% gives time derivatives of convolution zz
procedure gungc(a,p); 
if p={} then gungc_error(a,p) else -beta(first(p))*zz(a,p)+zz(a,rest(p));
% int by parts so all non-integrable convolutions are on one noise
procedure gungd(n,p,m,q); 
if (p={})or(q={}) then phi(n,p)*phi(m,q) else
(gungd(n,rest(p),m,q)+gungd(n,p,m,rest(q)))/(beta(first(p))+beta(first(q)));
% int by parts so all integrable convolutions update the approx
procedure gunge(n,p,m,q); 
if (p={})or(q={}) then 0 else
(-phi(n,p)*phi(m,q)+gunge(n,rest(p),m,q)+gunge(n,p,m,rest(q)))
/(beta(first(p))+beta(first(q)));
% similarly for zz ---integrable convolutions update the approx
procedure gungf(a,p);
if p={} then 0 else (gungf(a,rest(p))-zz(a,p))/beta(first(p));
let { zz(~a,{})=>a
    , df(zz(~a,~p),tt)=>gungc(a,p)
    , df(zz(~a,~p),t)=>df(zz(a,p),tt)
    , secular(sin(x)*zz(~a,~p),xx)
      =>secular(sin(x)*a,xx)/(for each m in p product beta(m))
    , secular(sin(~m*x)*zz(~a,~p),xx)=>0
    , secular(sin(x)*phi(~n,~p)*phi(~m,~q),xx)=>gungd(n,p,m,q)
    , secular(sin(x)*phi(~n,~p)^2,xx)=>gungd(n,p,n,p)
    , linv(sin(x)*phi(~n,~p)*phi(~m,~q),xx)=>gunge(n,p,m,q)*sin(x)
    , linv(sin(x)*phi(~n,~p)^2,xx)=>gunge(n,p,n,p)*sin(x)
    , linv(sin(~l*x)*phi(~n,~p)*phi(~m,~q),xx)
      =>sin(l*x)*zz(phi(n,p)*phi(m,q),{l})
    , linv(sin(~l*x)*phi(~n,~p)^2,xx)=>sin(l*x)*zz(phi(n,p)^2,{l})
    , linv(sin(~l*x)*zz(~a,~p),xx)=>sin(l*x)*zz(a,l.p)
    , linv(sin(x)*zz(~a,~p),xx)=>sin(x)*gungf(a,p)
    };

% now derive effects quadratic in noise
let sig^3=>0;
it:=0$
repeat begin
    deq:=trigsimp(-df(u,t)+df(u,x,2)+u-a*u^3+sig*noise*u,combine);
    g:=g+(gd:=secular(deq,xx));
    u:=u+linv(deq-gd*sin(x),xx);
    showtime;
end until ((it:=it+1)>9)or(deq=0);

end;%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Now transform the quadratic noise into equivalent new noises psi,
% but for now only up to two convolutions.
% New noises then have subscripts to uniquely identify them.
write "Now transforming the quadratic noises";
operator yy; linear yy;
operator psi; depend psi,tt,xx;
let { yy(1,tt)=>1,
      yy(phi(~i,{}),tt)=>phi(i,{}),
      yy(phi(~i,{})*phi(~j,{~k}),tt)
      => 1/2*(if i=j then 1 else 0)+psi(i,j,{k})/sqrt(2*beta(k)),
      yy(phi(~i,{})*phi(~j,{~k2,~k1}),tt)
      => (psi(i,j,{k1})/sqrt(2*beta(k1))+psi(i,j,{k2,k1})/sqrt(2*beta(k2)))
          /(beta(k1)+beta(k2))
    };
%eps:=1;
gg:=yy(g,tt);
% Root sum squares of the determined noise coefficients;
% it implicitly assumes that there is no correlation between terms.
operator yyy; linear yyy;
let { yyy(1,tt)=>0,
      yyy(psi(~i,~j,~p),tt)=>0,
      yyy(psi(~i,~j,~p)^2,tt)=>1,
      yyy(psi(~i,~j,~p)*psi(~ii,~jj,~pp),tt)=>0
    };
on rounded;
gg:=gg;
c20:=sqrt(yyy(coeffn(coeffn(gg,sig,2),a,0)^2,tt));
c21mean:=(coeffn(coeffn(gg,sig,2),a,1) where psi(~i,~j,~p)=>0);
c21:=sqrt(yyy(coeffn(coeffn(gg,sig,2),a,1)^2,tt));
off rounded;
showtime;

end;
