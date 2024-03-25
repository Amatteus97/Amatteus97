module hilbert_space
    implicit none
    integer, public, parameter :: dp = selected_real_kind(8), di = selected_int_kind(16)
PUBLIC :: input, binom, hilbert, base, index, integer2binary, binary2integer
integer(kind = di),dimension(: ), allocatable:: mu, md
integer, public :: lx, ly, nu, nd
integer(kind=di), public :: maxvalu, maxvald, d1, d2, minvalu, minvald, lmax
real (kind = dp), public :: u 
complex(10), public :: t_hop
contains

subroutine hilbert()
implicit none
integer :: l
!integer(kind=di) ::  i, v

call input()
l = lx*ly
allocate(mu(d1))
allocate(md(d2))
call base(minvalu, maxvalu, nu, d1,  mu)
call base(minvald, maxvald, nd, d2,  md)
 
   !call index(mu, md, arru, arrd)
! do v=1,d1
! print*, "mui", mu(v), "arrui", (arru(i,v), i=1,l) 
! print*, ""
! end do

! do v=1,d2
!     print*, "mdi", md(v), "arrdi",(arrd(i,v), i=1,l) 
!     print*, ""
!     end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !ricordarsi che qui dealloco solo perché è tutto il programma
   !quando verrà implementat tutto il resto deallocare
  !  !alla fine
  !   deallocate(mu)
  !   deallocate(md)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine hilbert

subroutine input()
implicit none
integer :: l, n 
integer (kind=di)::dimh
integer(kind=di) :: i
real(dp) :: tr, ti 
do
    
 print *, "Lx ="
 read *, lx
 print *, "Ly ="
 read *, ly
 print *, "Numero di spin up ="
 read *, nu
 print *, "Numero di spin down ="
 read *, nd
 l = lx*ly
 do
 if (l < 17) then
    exit
 else 
    print *, "Dimensione del reticolo troppo grande. Diminuire i valori di input"
    print *, "Lx ="
 read *, lx
 print *, "Ly ="
 read *, ly
 l = lx*ly
 end if
end do
 n=nu+nd
 do 
    if (n < l+1) then
        exit
     else 
        print *, "Numero di spin troppo grande. Diminuire i valori di input"
        print *, "Numero di spin up ="
    read *, nu
    print *, "Numero di spin down ="
    read *, nd
    n=nu+nd
     end if
 end do
 call binom(l,nu,d1) 
 call binom(l,nd, d2)
 dimh = d1*d2
 if (dimh < 166000000) then
    exit
 else 
    print *, "Dimensione dell spazio di Hilbert troppo grande. Diminuire i valori di input"
 end if
end do
maxvalu = 0
do i=l-nu,l-1
    maxvalu = maxvalu + 2**(i)
end do
maxvald = 0
do i=l-nd,l-1
    maxvald = maxvald + 2**(i)
end do
minvalu = 0
do i=0,nu-1
    minvalu = minvalu + 2**(i)
end do
minvald = 0
do i=0,nd-1
    minvald = minvald + 2**(i)
end do

print *, "Valore interazione ="
    read *, u

print *, "Numero massimo di passi Lanczos (minore o uguale a d1*d2, consigliato nell'ordine di 100)="
read *, lmax
lmax = min(lmax, d1*d2)

print *, "Parte reale hopping="
read *, tr
print *, "Parte immaginaria hopping="
read *, ti

t_hop = cmplx(tr, ti, 10)

print*, "lx ", lx, "ly ", ly, "l ", l,"nu ", nu, "nd ", nd, "u", u, "lmax", lmax, "dimh", dimh 

end subroutine input


subroutine index(mu, md, arru, arrd)
implicit none
integer(kind = di),dimension(d1),intent(out) :: mu
integer(kind = di),dimension(d2),intent(out) :: md
integer, dimension (lx*ly,d1), intent(in) :: arru
integer, dimension (lx*ly,d2), intent(in) :: arrd
integer (kind = di):: k
do k=1,d1
mu(k) = binary2integer(arru(:,k))
end do
do k=1,d2
md(k) = binary2integer(arrd(:,k))
end do
end subroutine index

function integer2binary(i, l) result(b)
    integer(kind = di),intent(in) :: i
    integer, intent(in) :: l
    integer, dimension(l) :: b
    integer(kind = di) :: k,j
    b=0
    j=i
    do k=1,l
      b(k)=int(mod(j,2_di))
      j=j/2
    enddo
  end function

  function binary2integer(b) result(i)
    integer (kind = di):: i
    integer,intent(in) :: b(:)
    integer k,j
    i=0
    j=1
    do k=1,size(b)
      i=i+b(k)*j
      j=j+j
    enddo
  end function


subroutine base(minval, maxval, n, d, ind)
implicit none

integer(kind = di),intent(in) :: minval, maxval, d
integer, dimension ( lx*ly ) :: temp 
integer(kind = di), dimension ( d ), intent(inout) :: ind
integer :: l, v, x, i
integer, intent(in) :: n 
integer(kind = di) :: m

l=lx*ly
v=1
do m=minval, maxval

temp = integer2binary(m, l)
x = sum(temp)
if ( x == n) then
    do i=1,l
      !  arr(i, v) = temp(i)
       ! print *, "n", n, "arr(i,v)",arr(i,v)," temp(i)",  temp(i)
        ind(v) =m
        
    end do
    v=v+1
    !print *, "v", v-1, "arr(,v)",(arr(i,v-1), i=1,l)," temp",  temp

end if
end do
end subroutine base


subroutine binom(n, k, bin)

    implicit none
    integer(kind = di) :: temp1, temp2, i
    integer, intent(in) :: n, k
    integer(kind=di), intent(inout) :: bin
    
    temp1=1
     temp2=1
    IF (k == 0) THEN
        BIN= 1
    end if
        IF (k ==1) then
            BIN = n
          END IF
            do i = k+1, n
                temp1= temp1*i
            end do
            do i = 2, n-k
                temp2= temp2*i
            end do
          bin = temp1/temp2
         
              
  end subroutine binom



end module hilbert_space

module hamilton_matrix
  use hilbert_space
  implicit none
  PUBLIC :: potential, getNeighbours, hoplist, quick_search, eimac

contains

subroutine getNeighbours(site, neighbours) !!!pbc
  INTEGER, intent(in):: site
  INTEGER,DIMENSION(4), intent(out) :: neighbours
  integer :: l
  l=lx*ly

  ! neighbours(1) = MODULO(site + lx,l +1) + FLOOR(REAL(site+lx-1)/l)
  ! neighbours(2) = MODULO(site - lx + l,l +1)+ CEILING(REAL(site-lx)/l)
  ! neighbours(3) = (site + 1) - lx*(CEILING(REAL(site+1)/lx)-CEILING(REAL(site)/lx))
  ! neighbours(4) = (site - 1) - lx*(CEILING(REAL(site-1)/lx)-CEILING(REAL(site)/lx))


  neighbours(1) = site + 1
if ( mod(site, lx) == 0 ) neighbours(1) = site + 1 - lx
neighbours(3) = site - 1
if ( mod(site, lx) == 1 ) neighbours(3) = site - 1 + lx
neighbours(2) = site + lx
if ( site > lx*(ly-1) ) neighbours(2) = site - lx*(ly-1)
neighbours(4) = site - lx
if ( site < lx+1) neighbours(4) = lx*(ly-1) + site 
  return
END subroutine

subroutine potential(mu,  md, pot)
  implicit none
  integer(kind = di), dimension ( d1 ), intent(in) :: mu
  integer(kind = di), dimension ( d2 ), intent(in) :: md
  integer, dimension(lx*ly) :: confu, confd
  integer, dimension(d1*d2), intent(out) :: pot 
  integer :: j, l
  integer(kind = di) :: ind, i,k
  l=lx*ly


  do i=1,d1
      do k=1,d2
          confu=integer2binary(mu(i), l)
          confd=integer2binary(md(k), l)
          ind = (k-1)*d1 + i 
          do j=1,l
            if(confu(j) == 1 .and. confd(j) == 1 )  then
              pot(ind) = pot(ind) + 1
              !print *, "pot-i-k", i, k 
            end if
          end do
      end do
  end do
    !print *, "potential", (pot(i), i=1,d1*d2)
  
  end subroutine potential
  

  subroutine hoplist(m, d,n, conn, cin)
  implicit none
  integer(kind = di), intent(in) :: d
  integer :: l, k, j, jnew, icount, n, r
  integer(kind = di), dimension ( d ), intent(in) :: m
  integer(kind = di), dimension(d,4*n), intent(out) :: conn
  complex(10), dimension(d,4*n), intent(out) :: cin
  integer, dimension(lx*ly) :: conf
  INTEGER,DIMENSION(4) :: neigh
  integer(kind = di) :: i , mnew, f
  complex(10) :: t
    t = t_hop
  l=lx*ly
  
  
    do i=1,d
      icount = 1
      conf=integer2binary(m(i), l)
      do j=1,l
        if(conf(j) == 1) then
          call getNeighbours(j, neigh)
          do k = 1, 4
            if ( conf(neigh(k)) == 0 ) then
              jnew = neigh(k)
              mnew = m(i) + (2**(jnew-1)) - (2**(j-1))
              !print *, "hoplist-mnew-mi", mnew, m(i)
              call quick_search(mnew, m, d, f)
              conn(i, icount) = f
              call eimac(j, jnew, conf, r)
              if ( k > 2 ) then
                cin(i,icount) = ((-1)**r)*conjg(t)
              else
                cin(i,icount) = ((-1)**r)*t 
              end if
              icount = icount + 1
              !print *, "hoplist-i-j-icount", i, j, icount
            end if
          end do
        end if
      end do
    end do


  end subroutine hoplist


  subroutine quick_search(target, list, d, pos)
  implicit none
  integer(kind = di), intent(in) :: target, d
  integer(kind = di), dimension ( d ), intent(in) :: list 
  integer(kind = di), intent(out) :: pos 
  integer(kind = di) ::  k, imax, imin
  
    k = nint(real(d)/2)
    imin=1
    imax=d

 do while (list(k) .ne. target)
   if ( list(k) < target ) then
    imin=k
   
    k = nint(real(imax+imin)/2)
  end if
  if ( list(k) > target ) then
  
    imax=k
    k = nint(real(imax+imin)/2)
  end if
  if ( (imax-imin) < 3 ) then
    do k = imin, imax
      if (list(k) == target) exit
    end do
  end if
  !print *, "q-s-imin-imax", imin, imax, list(k), target 
end do
!print *, "quick search", list(k), target 
pos = k 

  end subroutine quick_search

  subroutine eimac(start, finish, conf, part)
  implicit none
  integer, intent(in) :: start, finish
  integer, dimension(lx*ly), intent(in) :: conf
  integer, intent(out) :: part
  integer :: s, f, i

  part = 0
  s = min(start, finish) +1
  f = max(start, finish) -1

  do i=s,f
    if ( conf(i) == 1 ) then
      part = part +1
      
    end if
  end do

  end subroutine eimac



  subroutine hamilton(connu, cinu, connd, cind, pot)
  implicit none
  integer(kind = di), dimension(d1,4*nu), intent(out) :: connu
  complex(10), dimension(d1,4*nu), intent(out):: cinu
  integer(kind = di), dimension(d2,4*nd), intent(out) :: connd
  complex(10), dimension(d2,4*nd), intent(out) :: cind
  integer, dimension(d1*d2), intent(out):: pot

    call potential(mu,  md, pot)
    call hoplist(mu, d1,nu, connu, cinu)
    call hoplist(md, d2,nd, connd, cind)

  
  end subroutine hamilton
  
end module hamilton_matrix

module lanczos_diagonal
  use hilbert_space
  use hamilton_matrix
  implicit none
  public :: init_random_seed, lanczos_alg, hpsi, PRINT_MATRIX, permu, eigenvalue
  
contains

subroutine init_random_seed()
  use iso_fortran_env, only: int64
  implicit none
  integer, allocatable :: seed(:)
  integer :: i, n, un, istat, dt(8), pid
  integer(int64) :: t

  call random_seed(size = n)
  allocate(seed(n))
  ! First try if the OS provides a random number generator
  open(newunit=un, file="/dev/urandom", access="stream", &
       form="unformatted", action="read", status="old", iostat=istat)
  if (istat == 0) then
     read(un) seed
     close(un)
  else
     ! Fallback to XOR:ing the current time and pid. The PID is
     ! useful in case one launches multiple instances of the same
     ! program in parallel.
     call system_clock(t)
     if (t == 0) then
        call date_and_time(values=dt)
        t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
             + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
             + dt(3) * 24_int64 * 60 * 60 * 1000 &
             + dt(5) * 60 * 60 * 1000 &
             + dt(6) * 60 * 1000 + dt(7) * 1000 &
             + dt(8)
     end if
     pid = getpid()
     t = ieor(t, int(pid, kind(t)))
     do i = 1, n
        seed(i) = lcg(t)
     end do
  end if
  call random_seed(put=seed)
contains
  ! This simple PRNG might not be good enough for real work, but is
  ! sufficient for seeding a better PRNG.
  function lcg(s)
    integer :: lcg
    integer(int64) :: s
    if (s == 0) then
       s = 104729
    else
       s = mod(s, 4294967296_int64)
    end if
    s = mod(s * 279470273_int64, 4294967291_int64)
    lcg = int(mod(s, int(huge(0), int64)), kind(0))
  end function lcg
end subroutine init_random_seed


subroutine lanczos_alg(connu, cinu, connd, cind, pot)
implicit none
  integer(kind = di), dimension(d1,4*nu), intent(in) :: connu
  complex(kind = 10), dimension(d1,4*nu), intent(in):: cinu
  integer(kind = di), dimension(d2,4*nd), intent(in) :: connd
  complex(kind = 10), dimension(d2,4*nd), intent(in) :: cind
  integer, dimension(d1*d2), intent(in):: pot
  complex(kind = 10), dimension(d1*d2) :: v
  !integer :: i

call rnd_cplx16_norm(d1*d2,  v)

!print*, "v", (v(i), i=1,10)

call lanczos_eigenvalue(v, connu, cinu, connd, cind, pot)

end subroutine lanczos_alg

subroutine rnd_cplx16_norm(n,  v)
implicit none

       INTEGER(di), intent(in) :: n
       COMPLEX(10), dimension(n), intent(out) ::   v
       REAL(dp) :: norm
       INTEGER(di) :: i
       real(dp), dimension(:), allocatable :: u
    allocate(u(2*n))
       call init_random_seed()
       call random_number(u)
!print *, "n", n
     norm = 0
             DO i = 1, n
                v( i ) = cmplx( u(i), u(i+n) , 10) 
             end do
deallocate(u)
             DO i = 1, n
              norm = norm + ((real(v(i), dp))**2) + ((aimag(v(i)))**2)
             end do
           
             do i = 1, n
              v(i) = v(i)/sqrt(norm)
            end do

end subroutine rnd_cplx16_norm


subroutine hpsi(v1, v2, connu, cinu, connd, cind, pot)
implicit none
  integer(kind = di), dimension(d1,4*nu), intent(in) :: connu
  complex(kind = 10), dimension(d1,4*nu), intent(in):: cinu
  integer(kind = di), dimension(d2,4*nd), intent(in) :: connd
  complex(kind = 10), dimension(d2,4*nd), intent(in) :: cind
  integer, dimension(d1*d2), intent(in):: pot
  complex(kind = 10), dimension(d1*d2) :: v1, v2
  integer(kind = di) :: i, m, k, j
  integer :: l

v2 = 0
do i=1,d1*d2
  v2(i) = v2(i) + u*pot(i)*v1(i)
end do
do j = 1, d2
  do i = 1, d1
    k = (j-1)*d1 + i
    do l = 1, 4*nu
      if ( connu(i,l) .ne. 0 ) then
        m = (j-1)*d1 + connu(i,l)
        v2(m) = v2(m) + cinu(i,l)*v1(k)
      end if     
    end do
    do l = 1, 4*nd
      if ( connd(j,l) .ne. 0 ) then
        m = (connd(j,l) - 1)*d1 + i
        v2(m) = v2(m) + cind(j,l)*v1(k)
      end if     
    end do
  end do
end do


end subroutine hpsi

subroutine lanczos_eigenvalue(v, connu, cinu, connd, cind, pot)
implicit none

!!!!!!!!!!!!!!!!!!!!!!!allocare tutto insieme!!!!!!!!!!!!!!
  integer(kind = di), dimension(d1,4*nu), intent(in) :: connu
  complex(kind = 10), dimension(d1,4*nu), intent(in):: cinu
  integer(kind = di), dimension(d2,4*nd), intent(in) :: connd
  complex(kind = 10), dimension(d2,4*nd), intent(in) :: cind
  integer, dimension(d1*d2), intent(in):: pot
  complex(kind = 10), dimension(d1*d2), intent(inout) :: v
  complex(kind = 10), dimension(d1*d2) :: t, w
  real (kind = dp), dimension(lmax) :: alfa
  real (kind = dp), dimension(lmax) :: beta
  real (kind = dp), dimension(:), allocatable ::alfa_tmp, beta_tmp
  integer(kind = di) :: k, i
  real (kind = dp), dimension(:, :), allocatable :: eig
  real (kind = dp) :: etemp
  real(dp) ::  epsilon 
    beta = 0
    alfa = 0
    w = 0
    etemp = 0
    epsilon = 0.0000000001
  
    DO k=1, lmax
       call hpsi(v, t, connu, cinu, connd, cind, pot)
       alfa(k) = real(dot_product(v,t), dp)
       if ( k > 1 ) then
        t = t - alfa(k)*v - beta(k-1)*w
       else
        t = t - alfa(k)*v
       end if
    
        IF (k < lmax) then
          beta(k) = real(dot_product(t, t), dp)
          beta(k) = sqrt(beta(k))
          w = v
          write(12, *) "beta", beta
          write(13, *) "alfa", alfa
            v = t/beta(k)
         END IF
 
        !call permu(v, t, w)

        if ( k>1 ) then
          allocate(alfa_tmp(k), beta_tmp(k-1), eig(k,k))
          alfa_tmp = alfa
          do i = 1,k-1
            beta_tmp(i) = beta(i)
          end do
          CALL eigenvalue(alfa_tmp, beta_tmp, eig, k)
          write(11, *) "eig", alfa_tmp(1), abs(alfa_tmp(1) - etemp)
           if ( abs(alfa_tmp(1) - etemp) < epsilon) then
            twopi = 4*asin(1.0_dp)/lx
            do i = -lx, lx
                test(i) = cos(twopi*i)
                write(11, *) i, -2*real(t_hop)*test(i)
            end do
            do i = 1, k
                write(12, *) i, alfa_tmp(i)
            end do
             exit 
           end if
           !print *, "eig", alfa_tmp(1)
           etemp = alfa_tmp(1)
         deallocate(alfa_tmp, beta_tmp, eig)
        end if
       
       
       end do
      


end subroutine lanczos_eigenvalue

subroutine eigenvalue(alfa, beta, H, n)
implicit none
      INTEGER(di), intent(in) :: n
      real(dp), dimension(n, n), intent(inout) :: H
      real(dp), dimension(n), intent(inout):: alfa
      real(dp), dimension(n-1), intent(inout) :: beta
      INTEGER(di) ::  INFO
      real(dp), dimension(:), allocatable:: WORK 
      Allocate (work(2*n-2))
print *, "eccomi1", n 
!print*, alfa
!print *, beta
      call dstev('v', n, alfa, beta, H, n, WORK, INFO)
      print *, "info", info 
      print *, "eccomi2"
      if ( INFO == 0 ) then
        print *, "va bene" 
      else
        print *, "non va bene"
      end if
deallocate(work)


     
end subroutine eigenvalue


! subroutine eigenvalue(alfa, beta, z, k)
!   implicit none
!         INTEGER(di), intent(in) :: k
!         real(dp), dimension(lmax, lmax), intent(inout) :: z
!         real(dp), dimension(lmax), intent(inout):: alfa
!         real(dp), dimension(lmax), intent(inout) :: beta
!         !integer(di) :: i
!         ! print *, "eccomi1", k 
!         ! print*, alfa(k)
!         ! print *, beta(k-1)
!         call tqli(alfa,beta,k,lmax,z)
!         !print *, "eccomi2"
  
!   end subroutine eigenvalue
  

! SUBROUTINE tqli(d,e,n,np,z)


!   ! QL algorithm with implicit shifts, to determine the eigenvalues and eigenvectors of a real,
!   ! symmetric, tridiagonal matrix, or of a real, symmetric matrix previously reduced by tred2
!   ! §11.2. d is a vector of length np. On input, its first n elements are the diagonal elements of
!   ! the tridiagonal matrix. On output, it returns the eigenvalues. The vector e inputs the subdiagonal
!   ! elements of the tridiagonal matrix, with e(1) arbitrary. On output e is destroyed.
!   ! When finding only the eigenvalues, several lines may be omitted, as noted in the comments.
!   ! If the eigenvectors of a tridiagonal matrix are desired, the matrix z (n by n matrix stored
!   ! in np by np array) is input as the identity matrix. If the eigenvectors of a matrix that has
!   ! been reduced by tred2 are required, then z is input as the matrix output by tred2. In
!   ! either case, the kth column of z returns the normalized eigenvector corresponding to d(k).


!   implicit none
!   INTEGER(di), intent(in) :: n,np
!   REAL(dp), dimension(np), intent(inout) :: d
!   REAL(dp), dimension(np), intent(inout) :: e
!   REAL(dp), dimension(np,np), intent(inout) :: z
!   INTEGER(di) :: i,iter,k,l,m
!   REAL(dp) :: b,c,dd,f,g,p,r,s
 
!   do  l=1,n
!       iter=0
!   1   do m=l,n-1 
!           dd=abs(d(m))+abs(d(m+1)) 
!           if (abs(e(m))+dd.eq.dd) exit 
!       end do 
!       if (abs(e(m))+dd.ne.dd) m=n
!       if(m.ne.l)then
!           iter=iter+1
!           g=(d(l+1)-d(l))/(2.*e(l)) 
!           r=sqrt((g**2) + 1)
!           g=d(m)-d(l)+e(l)/(g+sign(r,g)) 
!           s=1.
!           c=1.
!           p=0.
!           do  i=m-1,l,-1 
!             print*, "ieig", i 
!               f=s*e(i)
!               b=c*e(i)
!               print*, "geig", g
!               print*, "feig", f
!               r=sqrt((g**2) + (f**2))
!               e(i+1)=r
!               if(r.eq.0.)then 
!                   d(i+1)=d(i+1)-p
!                   e(m)=0.
!                   goto 1
!               end if
!               s=f/r
!               c=g/r
!               g=d(i+1)-p
!               r=(d(i)-g)*s+2.*c*b
!               p=s*r
!               d(i+1)=g+p
!               g=c*r-b
! !!!!!!!!!!!!!!!!!!!!!autovettori!!!!!!
!               do  k=1,n 
!                   f=z(k,i+1)
!                   z(k,i+1)=s*z(k,i)+c*f
!                   z(k,i)=c*z(k,i)-s*f
!               end do 
!               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           end do 
!           d(l)=d(l)-p
!           e(l)=g
!           e(m)=0.
!           goto 1
!       end if
!   end do 
!   END SUBROUTINE tqli


SUBROUTINE PRINT_MATRIX( DESC, M, N, A, LDA )
  CHARACTER*(*)    DESC
  INTEGER(di)          M, N, LDA
  COMPLEX(10)       A( LDA, * )

  INTEGER(di)          I, J

  WRITE(*,*)
  WRITE(*,*) DESC
  DO I = 1, M
     WRITE(*,9998) ( A( I, J ), J = 1, N )
  END DO

9998 FORMAT( 11(:,1X,'(',F6.2,',',F6.2,')') )
  
  END subroutine PRINT_MATRIX

  subroutine permu(v1, v2, v3)
  implicit none
 
  complex(kind = 10), dimension(d1*d2), intent(inout) :: v1, v2, v3
  complex(kind = 10), dimension(:), allocatable :: t1, t2, t3


    allocate(t1(d1*d2), t2(d1*d2), t3(d1*d2))
    t1 = v1
    t2 = v2
    t3 = v3

    v1 = t3
    v2 = t1
    v3 = t2

    deallocate(t1, t2, t3)

  end subroutine permu


end module lanczos_diagonal

program hilb
    
    use hilbert_space
    use hamilton_matrix
    use lanczos_diagonal
    integer(kind = di), dimension(:,:), allocatable :: connu
    complex(10), dimension(:,:), allocatable:: cinu
    integer(kind = di), dimension(:,:), allocatable :: connd
    complex(10), dimension(:,:), allocatable :: cind
    integer, dimension(:), allocatable:: pot
    call hilbert()
    allocate(connu(d1,4*nu), connd(d2,4*nd), cinu(d1,4*nu), cind(d2,4*nd), pot(d1*d2))
    !!!!!!!!!!allocare vettori lanczos!!!!!!!!!!!!!!!!!!!
    call hamilton(connu, cinu, connd, cind, pot)
    call lanczos_alg(connu, cinu, connd, cind, pot)


  end program hilb