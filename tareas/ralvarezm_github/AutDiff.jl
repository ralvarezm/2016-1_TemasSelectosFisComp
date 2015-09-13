###  Módulo que implementa números Duales      ###
###  para obtener la primera derivada de una   ###
###  función usando diferenciación automática. ###

module AutDiff
export Dual,dual_cte,dual_var
immutable Dual{T<:Real} <:Number
    f::T
    fp::T
end

Dual(f::Real,fp::Real)=Dual(promote(f,fp)...)
dual_var(x0::Real)=Dual(x0,1.0)
dual_cte(x0::Real)=Dual(x0,0.0)

import Base.+
+(x::Dual,y::Dual)=Dual(x.f+y.f,x.fp+y.fp)
+(x::Dual,y::Real)=Dual(x.f+y,x.fp)
+(x::Real,y::Dual)=Dual(x+y.f,y.fp)

import Base.-
-(x::Dual,y::Dual)=Dual(x.f-y.f,x.fp-y.fp)
-(x::Dual,y::Real)=Dual(x.f-y,x.fp)
-(x::Real,y::Dual)=Dual(y.f-x,y.fp)

import Base.*
*(x::Dual, y::Dual) = Dual(x.f*y.f,x.fp*y.f+y.fp*x.f)
*(x::Dual, y::Real) = Dual(x.f*y,x.fp*y)
*(y::Real, x::Dual) = Dual(x.f*y,x.fp*y)

import Base./
/(x::Dual, y::Dual) = Dual(x.f/y.f,(x.fp-(x.f/y.f)y.fp)/(y.f))
/(x::Dual, y::Real) = Dual(x.f/y,x.fp/y)
/(y::Real, x::Dual) = Dual(y/x.f,(y/x.f)x.fp/(x.f))

import Base.^
^(x::Dual,y::Dual) = Dual(x.f^y.f,y.f*x.f^(y.f-1)*x.fp)
^(x::Dual,y::Integer) = Dual(x.f^y,y*x.f^(y-1)*x.fp)
^(x::Dual,y::Real) = Dual(x.f^y,y*x.f^(y-1)*x.fp)

import Base.abs

abs(x::Dual)=Dual(abs(x.f),x.f*(sign(x.f)/sign(x.f)))  #Gracias al pull request de Ceboc pude corregir este error.

import Base.sin
sin(x::Dual)=Dual(sin(x.f),cos(x.f)*x.fp)

import Base.cos
cos(x::Dual)=Dual(cos(x.f),-sin(x.f)*x.fp)

import Base.tan
tan(x::Dual)=Dual(tan(x.f),sec(x.f)^2 * x.fp)

import Base.cot
cot(x::Dual)=Dual(cot(x.f),-csc(x.f)^2 * x.fp)

import Base.sec
sec(x::Dual)=Dual(sec(x.f),sec(x.f)*tan(x.f)*x.fp)

import Base.csc
csc(x::Dual)=Dual(csc(x.f),-csc(x.f)*cot(x.f)*x.fp)

import Base.asin
asin(x::Dual)=Dual(asin(x.f),x.fp/sqrt(1-x.f^2))

import Base.acos
acos(x::Dual)=Dual(acos(x.f),-x.fp/sqrt(1-x.f^2))

import Base.atan
atan(x::Dual)=Dual(atan(x.f),x.fp/(x.f^2+1))

import Base.acot
acot(x::Dual)=Dual(acot(x.f),-x.fp/(x.f^2+1))

import Base.asec
asec(x::Dual)=Dual(asec(x.f),x.fp/(abs(x.f)*sqrt(x.f^2-1)))

import Base.acsc
acsc(x::Dual)=Dual(acsc(x.f),-x.fp/(abs(x.f)*sqrt(x.f^2-1)))

import Base.exp
exp(x::Dual)=Dual(exp(x.f),exp(x.f)*x.fp)

import Base.sqrt
sqrt(x::Dual)=Dual(sqrt(x.f),x.fp/(2*sqrt(x.f)))

import Base.log
log(x::Dual)=Dual(log(x.f),x.fp/x.f)

import Base.sinh
sinh(x::Dual)=Dual(sinh(x.f),cosh(x.f)*x.fp)

import Base.cosh
cosh(x::Dual)=Dual(cosh(x.f),sinh(x.f)*x.fp)

import Base.tanh
tanh(x::Dual)=sinh(x)/cosh(x)

import Base.coth
coth(x::Dual)=cosh(x)/sinh(x)

import Base.sech
sech(x::Dual)=1/cosh(x)

import Base.csch
csch(x::Dual)=1/sinh(x)

import Base.asinh
asinh(x::Dual)=Dual(asinh(x.f),x.fp/sqrt(x.f^2+1))

import Base.acosh
acosh(x::Dual)=Dual(acosh(x.f),x.fp/sqrt(x.f^2-1))

import Base.atanh
atanh(x::Dual)=Dual(atanh(x.f),x.fp/(1-x.f^2))

import Base.acoth
acoth(x::Dual)=Dual(acoth(x.f),x.fp/(1-x.f^2))

import Base.asech
asech(x::Dual)=Dual(asech(x.f),-x.fp/(x.f*sqrt(1-x.f^2)))

import Base.acsch
acsch(x::Dual)=Dual(acsch(x.f),-x.fp/(abs(x.f)*sqrt(1+x.f^2)))

end
