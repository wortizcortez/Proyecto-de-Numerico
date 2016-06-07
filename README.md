# Proyecto-de-Numerico
este reporte se a creado para comprender y entender lo smetodos de cuadratura, brent y de monte carlos....temas de analisi numerico 
Método de brent
function [p itmax] = brent(funct,w1,w2,TOL)

%valores de entrada:
%	funct: Es una función a la cual se desea calcular su una aproximación a su raíz.
%	w1 y w2: Son los extremos del intervalo que contiene a la raíz de la función. 
%	TOL: Es la tolerancia que deseamos tener en la aproximación de la raíz.

%Valores de salida:
%	la función devuelve un vector [p itmax]; donde p es la aproximación a la raiz y itmax es la cantidad de veces que se evaluo la función.

ITMAX=10000;% valor máximo de iteraciones.
EPS=3.0*(10^-8); % presición para punto flotante
a=w1;
b=w2;
fa=feval(funct,a);
fb=feval(funct,b);
if((fa>0) & (fb>0)) || ((fa<0) & (fb<0));
    error('Procedimiento terminado sin exito');
else
    c=b;
    fc=fb;
    iter=0; % para guardar el numero de iteraciones
    while iter<ITMAX;
        if((fb>0) & (fc>0)) || ((fb<0) & (fc<0));
            c=a; 
            fc=fa;
            d=(b-a);
            e=d;
        end
        if abs(fc)<abs(fb);
            a=b;
            b=c;
            c=a;
            fa=fb;
            fb=fc;
            fc=fa;
        end
        TOL1=2*EPS*abs(b)+ 0.5*TOL; %Observa la convergencia 
        xm=0.5*(c-b);
        if(abs(xm)<=TOL(1) || (fb==0)
            p=b;
            itmax=iter;
            break
        end
        if(abs(e)>=TOL1 && abs(fa)>abs(fb))
            s=fb/fa; %interpolacion inversa de la cuadratica 
            if(a==c);
                p=2*xm*s;
                q=1-s;
            else
                q=fa/fc;
                r=fb/fc;
                p=s*(2*xm*q*(q-r)-(b-a)*(r-1));
                q=(q-1)*(r-1)*(s-1);
            end
            if(p>0); %en los limites se comprueba
                q=-q;
                p=abs(p);
            end
            if(2*p < min(3*xm*q-abs(TOL1*q),abs(e*q)));
                e=d; %interpolacion aprobada
                d=p/q;
            else
                d=xm; %se usa la biseccion 
                e=d;
            end
        else 
            d=xm;
            e=d;
        end
        a=b; 
        fa=fb;
        if(abs(d)>TOL1); %se evalua la nueva raiz
            b=b+d;
        else
            if sign(xm)==0;
                b=b+abs(TOL1)*(-1);
            else
                b=b+abs(TOL1);
            end
        end
        fb=feval(funct,b);
        iter=iter+1;
    end
    if iter==ITMAX
    	error ('La funcion Brent excedio maximo iteraciones');
    end
    p=b;
    itmax=iter;
end
