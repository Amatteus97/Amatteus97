SUBROUTINE tqli(d,e,n,np,z)
INTEGER n,np
REAL d(np),e(np),z(np,np)

INTEGER i,iter,k,l,m
REAL b,c,dd,f,g,p,r,s,pythag

do  i=2,n 
    e(i-1)=e(i)
end do
e(n)=0.
do  l=1,n
    iter=0
1   do m=l,n-1 
        dd=abs(d(m))+abs(d(m+1)) 
        if (abs(e(m))+dd.eq.dd) goto 2
    end do 
    m=n
2   if(m.ne.l)then
        iter=iter+1
        g=(d(l+1)-d(l))/(2.*e(l)) 
        r=pythag(g,1.)
        g=d(m)-d(l)+e(l)/(g+sign(r,g)) 
        s=1.
        c=1.
        p=0.
        do  i=m-1,l,-1 
            f=s*e(i)
            b=c*e(i)
            r=pythag(f,g)
            e(i+1)=r
            if(r.eq.0.)then 
                d(i+1)=d(i+1)-p
                e(m)=0.
                goto 1
            end if
            s=f/r
            c=g/r
            g=d(i+1)-p
            r=(d(i)-g)*s+2.*c*b
            p=s*r
            d(i+1)=g+p
            g=c*r-b
            do  k=1,n 
                f=z(k,i+1)
                z(k,i+1)=s*z(k,i)+c*f
                z(k,i)=c*z(k,i)-s*f
            end do 
        end do 
        d(l)=d(l)-p
        e(l)=g
        e(m)=0.
        goto 1
    end if
end do 
return
END SUBROUTINE tqli