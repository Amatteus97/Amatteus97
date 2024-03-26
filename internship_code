module hilbert_space
    implicit none
integer, public, parameter :: dp = selected_real_kind(8), di = selected_int_kind(16), bi = selected_int_kind(33)
PUBLIC :: input, binom, interchange_sort
public :: lanczos_alg, hpsi, eigenvalue, nn_int
PUBLIC :: getNeighbours, hoplist
integer(kind = di),dimension(: ), allocatable, public:: mu
integer, public :: lx, nu, derind
integer(kind=di), public :: d1, lmax
real (kind = dp), public :: twopi, vnn
complex(8), public :: mag_phase
real(dp), public :: cappa, deltac
real(dp), dimension(:), allocatable, public :: deri
contains


subroutine input()
implicit none
integer :: l, n 
integer (kind=di)::dimh
do
    
 print *, "Lx ="
 read *, lx
!lx = 4
  l = lx
  do
 if (l < 101) then
    exit
 else 
    print *, "Dimensione del reticolo troppo grande. Diminuire i valori di input"
    print *, "Lx ="
 read *, lx
 l = lx
 end if
end do

print *, "Numero di spin up ="
read *, nu

 n=nu
 do 
    if (n < l+1) then
        exit
     else 
        print *, "Numero di spin troppo grande. Diminuire i valori di input"
        print *, "Numero di spin up ="
    read *, nu
    n=nu
     end if
 end do
 call binom(l,nu,d1) 
 dimh = d1
 if (dimh < 3921225) then
    exit
 else 
    print *, "Dimensione dell spazio di Hilbert troppo grande. Diminuire i valori di input"
 end if
end do
    print *, "Valore interazione primi vicini="
    read *, vnn
!vnn = 1

    print*, "Valore deltac"
    read *, deltac
    !deltac = 0.000001

     twopi = 4*asin(1.0_dp)
     
!print *, "Numero massimo di passi Lanczos (minore o uguale a d1*d2, consigliato nell'ordine di 100)="
!read *, lmax
     lmax = 500
lmax = min(lmax, d1)

print*, "lx ", lx,"nu ", nu, "vnn", vnn, "lmax", lmax, "dimh", dimh

end subroutine input


  subroutine dritto(nsite,nel,ivp,m,pos)
    implicit none
    INTEGER(4), intent(in) :: nsite,nel
    INTEGER(di), intent(in) :: ivp(nsite,nsite,nel), m
    integer(di) :: vp, vp1, a
    INTEGER(4) :: ll,mm,j,p
    INTEGER(4), intent(out) :: pos(nel)

    a=m
    ll=nsite
    mm=nel
    j=1

    do while(mm.gt.0)
     vp=0
     do p=1,ll-mm+1
      vp1=vp
      vp=ivp(p,ll,mm)
      if(a.le.vp) then
       if(mm.eq.1) vp1=vp
       pos(j)=p
       if(j.ge.2) pos(j)=pos(j)+pos(j-1)
       go to 4
      endif
     enddo
4      continue
     ll=nsite-pos(j)
     mm=mm-1
     j=j+1
     a=a-vp1
    enddo

    return
    end subroutine dritto 


    subroutine rovescio(nsite,nel,ivp,pos,a, order, permu)
        implicit none
        INTEGER(4), intent(in) :: nsite,nel
        INTEGER(4) :: ll,mm,j,p
        integer, dimension(nel) :: npos 
        integer(di) :: vp 
        INTEGER(di), intent(in) :: ivp(nsite,nsite,nel)
        INTEGER(4), intent(inout) :: pos(nel), order(nel)
        INTEGER(di), intent(out) :: a
        integer, intent(out) :: permu

        npos = pos 

        call interchange_sort(1, nel, npos, pos, order, permu)
  
        p=pos(1)
        ll=nsite
        mm=nel
        j=1
        a=0
  
        do while(mm.gt.0)
         vp=0
         if(p.gt.1) vp=ivp(p-1,ll,mm)
         if(mm.eq.1) vp=vp+1
         a=a+vp
         if(j.lt.nel) p=pos(j+1)-pos(j)
         ll=nsite-pos(j)
         mm=mm-1
         j=j+1
        enddo
  
        return
    end subroutine rovescio 

    SUBROUTINE interchange_sort(left_end, right_end, list_in, list_out, order, permu)
    
    INTEGER, INTENT(IN) :: left_end, right_end
    integer, dimension(right_end), intent(in) :: list_in
    integer, DIMENSION (right_end), INTENT(OUT)  :: list_out
    INTEGER, DIMENSION (right_end), intent(inout)  :: order
    integer, intent(out) :: permu
    
    !     Local variables
    INTEGER             :: i, j, itemp
    integer                :: temp

  
    list_out = list_in
    DO i = 1, right_end
      order(i) = i
    END DO

    permu = 0
    
    DO i = left_end, right_end - 1
      DO j = i+1, right_end
        IF (list_out(i) > list_out(j)) THEN
          temp = list_out(i)
          list_out(i) = list_out(j) 
          list_out(j) = temp
          itemp = order(i) 
          order(i) = order(j)
          order(j) = itemp
          permu = permu + 1
        END IF
      END DO
    END DO
    
    END SUBROUTINE interchange_sort

    subroutine binom(n, k, bin)

    implicit none
    integer(kind = bi) :: temp1, temp2, i
    integer, intent(in) :: n, k
    integer(kind=di), intent(inout) :: bin
    
    temp1=1
     temp2=1
    
        do i = n-k + 1, n
            temp1= temp1*i
        end do
        do i = 2, k
            temp2= temp2*i
        end do
   
    bin = int(temp1/temp2, di)
         
              
  end subroutine binom


subroutine getNeighbours(site, neighbours) !!!pbc
    implicit none
  INTEGER, intent(in):: site
  INTEGER,DIMENSION(2), intent(out) :: neighbours

neighbours(1) = site + 1
if ( mod(site, lx) == 0 ) neighbours(1) = site + 1 - lx
neighbours(2) = site - 1
if ( mod(site, lx) == 1 ) neighbours(2) = site - 1 + lx
if (neighbours(1) == neighbours(2)) neighbours(2) = 0
return
 
END subroutine


  subroutine hoplist(d,n, conn, cin, ivp)
  implicit none
  integer(kind = di), intent(in) :: d
  integer, intent(in) :: n
  integer(kind = di), dimension(lx,lx,n), intent(in) :: ivp 
  integer :: l, j, icount, r, b
  integer(kind = di), dimension(d,n), intent(out) :: conn
  real(dp), dimension(d,n), intent(out) :: cin
  integer, dimension(:), allocatable :: conf, nconf
  INTEGER,DIMENSION(2) :: neigh
  integer(kind = di) :: i , f
  logical :: uscita 
  integer, dimension(n) :: order 
  !print *, "hoplist start"
  l=lx
  allocate(conf(n), nconf(n))
  
    do i=1,d
      icount = 1
      call dritto(l,n,ivp,i,conf)
      do j=1,n
          call getNeighbours(conf(j), neigh)
            nconf = conf 
            nconf(j) = neigh(1)
            uscita = .true.
            do b = 1, n
                if ( b == j ) cycle
                if ( conf(b) == nconf(j) ) then
                    uscita = .false.
                end if 
            end do
            if ( uscita ) then
                call rovescio(l,n,ivp,nconf,f, order, r)
                conn(i, icount) = f
                cin(i,icount) = ((-1)**r)
                icount = icount + 1
            end if
      end do
    end do
    !print *, "hoplist end"

  deallocate(conf, nconf)
  end subroutine hoplist


  subroutine hamilton(connu, cinu, vint, dvint, ivp)
  implicit none
  integer(kind = di), dimension(lx,lx,nu), intent(inout) :: ivp 
  integer(kind = di), dimension(d1,nu), intent(out) :: connu
  real(dp), dimension(d1,nu), intent(out):: cinu
  integer(di), intent(in) :: dvint
  integer, dimension(dvint), intent(out):: vint

    call nn_int(vint, dvint, ivp)
    call hoplist(d1,nu, connu, cinu, ivp)
  
  end subroutine hamilton
  

subroutine lanczos_alg(connu, cinu, v, t, w, alfa, beta, vint, dvint, ivp)
implicit none
  integer(di), intent(in) :: dvint
  integer(kind = di), dimension(d1,nu), intent(in) :: connu
  real(dp), dimension(d1,nu), intent(in):: cinu
  complex(kind = 8), dimension(d1), intent(out) :: v, t, w
  real (kind = dp), dimension(lmax), intent(out) :: alfa, beta
  integer, dimension(dvint), intent(out) :: vint
  integer(kind = di), dimension(lx,lx,nu), intent(inout) :: ivp

  

call rnd_cplx16_norm(d1,  v)

call lanczos_eigenvalue(v, connu, cinu, t, w, alfa, beta, vint, dvint, ivp)


end subroutine lanczos_alg

subroutine rnd_cplx16_norm(n,  v)
implicit none

       INTEGER(di), intent(in) :: n
       COMPLEX(8), dimension(n), intent(out) ::   v
       REAL(dp) :: norm
       INTEGER(di) :: i
       integer :: z
       integer, dimension(:), allocatable :: seed 
       real(dp), dimension(:), allocatable :: u
       allocate(u(2*n))
       call random_seed(size = z)
       allocate(seed(z))
       call random_seed(get = seed)
       call random_number(u)
     norm = 0
             DO i = 1, n
                v( i ) = cmplx( u(i), u(i+n) , 8) 
             end do
deallocate(u)
             DO i = 1, n
              norm = norm + ((real(v(i), dp))**2) + ((aimag(v(i)))**2)
             end do
           
             do i = 1, n
              v(i) = v(i)/sqrt(norm)
            end do

end subroutine rnd_cplx16_norm


subroutine hpsi(v1, v2, connu, cinu, vint, dvint, ivp)
implicit none
  integer(kind = di), dimension(d1,nu), intent(in) :: connu
  real(dp), dimension(d1,nu), intent(in):: cinu
  integer(di), intent(in) :: dvint
  complex(kind = 8), dimension(d1), intent(in) :: v1
  complex(kind = 8), dimension(d1), intent(out) :: v2
  integer(kind = di) :: i, m
  integer :: l
  integer, dimension(dvint), intent(in) :: vint
  complex(8) :: tx
  integer, dimension(:), allocatable :: conf
  integer(kind = di), dimension(lx,lx,nu), intent(inout) :: ivp

  allocate(conf(nu))

  tx = (1.0_dp, 0.0_dp)
v2 = 0
do i=1,d1
  v2(i) = v2(i) + vnn*vint(i)*v1(i)
end do
  do i = 1, d1
    do l = 1, nu
      if ( connu(i,l) .ne. 0 ) then
        m = connu(i,l)
        call dritto(lx,nu,ivp,i,conf)
        if ( conf(nu) == lx) then
            call dritto(lx,nu,ivp,m,conf)
            if ( conf(1) == 1) tx = mag_phase 
        end if
        v2(m) = v2(m) + cinu(i,l)*v1(i)*tx
        v2(i) = v2(i) + cinu(i,l)*v1(m)*conjg(tx)
      end if     
    end do
  end do

  deallocate(conf)



end subroutine hpsi

subroutine lanczos_eigenvalue(v, connu, cinu, t, w, alfa, beta, vint, dvint, ivp)
implicit none

  integer(kind = di), dimension(d1,nu), intent(in) :: connu
  real(dp), dimension(d1,nu), intent(in):: cinu
  integer(di), intent(in) :: dvint
  complex(kind = 8), dimension(d1), intent(inout) :: v
  complex(kind = 8), dimension(d1), intent(out) :: t, w
  real (kind = dp), dimension(lmax), intent(out) :: alfa, beta
  real (kind = dp), dimension(:), allocatable ::alfa_tmp, beta_tmp
  integer(di) :: k, i
  real (kind = dp), dimension(:, :), allocatable :: eig
  real (kind = dp) :: etemp
  real(dp) ::  epsilon 
  complex(kind = 8), dimension(d1) :: v_rnd
  integer, dimension(dvint), intent(inout) :: vint
  integer(kind = di), dimension(lx,lx,nu), intent(inout) :: ivp


    beta = 0
    alfa = 0
    w = 0
    etemp = 0
    epsilon = 0.000000000001
    !epsilon = 0.000001
    v_rnd = v
  
    DO k=1, lmax
        call hpsi(v, t, connu, cinu, vint, dvint, ivp)
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
           v = t/beta(k)
          END IF
  
        if ( k>1 ) then
          allocate(alfa_tmp(k), beta_tmp(k-1), eig(k,k))
          alfa_tmp = alfa
          do i = 1,k-1
            beta_tmp(i) = beta(i)
          end do
          CALL eigenvalue(alfa_tmp, beta_tmp, eig, k)
          !print*, k, alfa_tmp(1)
           if ( abs(alfa_tmp(1) - etemp) < epsilon) then
            print *, "gs", cappa, alfa_tmp(1)
         !   write(100+lx, *) twopi*cappa, alfa_tmp(1), -2*(cos(twopi*cappa) + cos(twopi/lx + twopi*cappa) &
           ! &+ cos(twopi/lx - twopi*cappa))
            deri(derind) = alfa_tmp(1)
             exit 
           end if
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
      INTEGER ::  INFO, ndstev
      real(dp), dimension(:), allocatable:: WORK 
      ndstev = int(n, 4)
      Allocate (work(2*n-2))
      call dstev('v', ndstev, alfa, beta, H, n, WORK, INFO)
      if ( INFO .ne. 0 ) then
        print *, "non va bene"
        print *, "info", info
      end if
deallocate(work)
    
end subroutine eigenvalue

subroutine nn_int(vint, dvint, ivp)
    implicit none
    integer(di) :: i
    integer :: l, j, dens, k, b
    integer, dimension(2) :: neigh
    integer, dimension(nu):: confu
    integer(di), intent(in) :: dvint
    integer(kind = di), dimension(lx,lx,nu), intent(in) :: ivp
    integer, dimension(dvint), intent(out) :: vint
    
    l=lx
    vint = 0
    
    do i=1,d1
        call dritto(l,nu,ivp,i,confu)
        dens = 0
        do j = 1, nu 
          call getNeighbours(confu(j), neigh) 
          do k = 1, 2
            do b = 1, nu
                if ( b == j ) cycle
                if ( lx == 2 .and. modulo(confu(j), 2) == 1 .and. k == 1) cycle 
                if ( confu(b) == neigh(k) ) vint(i) = vint(i) + 1
            end do
          end do 
        end do
    end do
    
    
    end subroutine nn_int


end module hilbert_space

program hilb
    
    use hilbert_space
    implicit none
    integer(kind = di), dimension(:,:), allocatable :: connu
    real(dp), dimension(:,:), allocatable:: cinu
    integer, dimension(:), allocatable:: vint
    complex(kind = 8), dimension(:), allocatable :: v, t, w
    real (kind = dp), dimension(:), allocatable :: alfa, beta
    integer :: i, ipmax, ill, imm, ip, ib, temp, b
    integer(di) :: dvint, bin
    INTEGER(di), dimension(:,:,:), allocatable :: ivpup
    real(dp) :: deri2, deri2temp2, deri2temp1
    
    
    call input()
    allocate(mu(d1)) 
    ALLOCATE(ivpup(lx,lx,nu))
    allocate(connu(d1,nu), cinu(d1,nu))
    allocate(alfa(lmax), beta(lmax), v(d1), t(d1), w(d1))
    dvint = d1
    allocate(vint(dvint))
    do ill=1,lx
        do imm=1,nu
         ipmax=ill-imm+1
         if(ipmax.ge.1)then
          do ip=1,ipmax
            ivpup(ip,ill,imm)=0
           do ib=1,ip
           if(ill-ib.gt.lx.or.imm-1.gt.nu)then
            write(6,*)'what',ill-ib,imm-1
            stop
           endif
           call binom(ill-ib, imm-1, bin)
            ivpup(ip,ill,imm)=ivpup(ip,ill,imm)+bin
           enddo
          enddo
         endif
        enddo
       enddo
      connu = 0
     ! open(unit = 20, file = 'ene.dat', status = 'replace', action = 'write' )
      call hamilton(connu, cinu, vint, dvint, ivpup)
      !temp = nint(1.0_dp/deltac)
      temp = 3
      allocate( deri(temp))
      do b = 0, 0
        !vnn = 0
          !if ( b == 0 ) vnn = 0
          !if ( b == 1 ) vnn = 1
          !if ( b == 2 ) vnn = 10
         ! write(100+lx, *) "#vnn=", vnn
          do i = 0,temp-1
            derind = i+1
            cappa = deltac*i 
            mag_phase = exp(cmplx(0.0_dp, - twopi * cappa, 8))
            call lanczos_alg(connu, cinu, v, t, w, alfa, beta, vint, dvint, ivpup)
          end do
         ! close(20)
          open(unit = 21, file = 'deri2.dat', status = 'replace', action = 'write' )
          open(unit =24, file = 'druweivcv0.dat', status = 'old', access = 'sequential', position='append')
          open(unit =22, file = 'druweivc.dat', status = 'old', access = 'sequential', position='append')
          open(unit =23, file = 'druweiline.dat', status = 'old', access = 'sequential', position='append')
          write(21, *) "#deltac =", deltac 
          do i = 1, temp-2
            deri2temp1 = deri(i+2)-deri(i+1)
            !if ( abs(deri2temp1) < 0.000000001) deri2temp1 = 0
            deri2temp2 = deri(i+1)-deri(i)
            !if ( abs(deri2temp2) < 0.000000001) deri2temp2 = 0
            deri2 = deri2temp1-deri2temp2
            !if ( i == 1 ) write(22, *) "#lx =", lx, "nu =", nu, "deltac=", deltac,"vnn", vnn, "drude weight=", &
            !&twopi*deri2/(2*twopi*twopi*deltac*deltac*lx)
            !if ( i == 1) write(24, *) "drude weight=", twopi*deri2/(2*twopi*twopi*deltac*deltac*lx), "expe", &
            !&twopi*(1+2*cos(twopi/lx))/lx
            !if ( i == 1) write(23, *) deltac, abs(twopi*deri2/(2*twopi*twopi*deltac*deltac*lx) - &
            !&twopi*(1+2*cos(twopi/lx))/lx)&
            !&/(twopi*deri2/(2*twopi*twopi*deltac*deltac*lx))
            !write(21, *) deltac*(i-1), deri2
          end do
          !write(23, *) "# lx =", lx, "nu=", nu 
        !  write(23,*)  real(nu)/lx, twopi*deri2/(2*twopi*twopi*deltac*deltac*lx) 
          close(21)
          close(23)
          close(24)
          close(22)
      end do
      

    deallocate(alfa, ivpup) 
    deallocate(beta, v, t, w) 
    deallocate(mu, vint)
    deallocate(connu)
    deallocate(cinu, deri) 

  end program hilb
