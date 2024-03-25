module hilbert_space
    implicit none
    integer, public, parameter :: dp = selected_real_kind(8), di = selected_int_kind(16), bi = selected_int_kind(33)
PUBLIC :: input, binom, interchange_sort
public :: lanczos_alg, hpsi, eigenvalue, nn_int
PUBLIC :: hoplist, smeno
integer(kind = di),dimension(: ), allocatable, public:: mu
integer, public :: lx, nu, Jsign
integer(kind=di), public :: d1, lmax
logical, public :: pbc 
contains


subroutine input()
implicit none
integer ::  n

Jsign = -1
do

print *, "Lx ="
 read *, lx
  !lx = 18
  do
 if (lx < 31 .and. lx > 0) then
    exit
 else 
    print *, "Dimensione del reticolo troppo grande. Diminuire i valori di input"
    print *, "Lx ="
 read *, lx
 end if
end do

print *, "Numero di spin up ="
read *, nu
!nu = lx/2
 n=nu
 do 
    if (n < lx+1) then
        exit
     else 
        print *, "Numero di spin troppo grande. Diminuire i valori di input"
        print *, "Numero di spin up ="
    read *, nu
    n=nu
     end if
 end do

 print *, "Pbc?"
 read *, pbc 
! pbc = .false.
!d1 = binom(lx,nu)
 call binom(lx,nu, d1) 
 if (d1 < 3921225) then
    exit
 else 
    print *, "Dimensione dell spazio di Hilbert troppo grande. Diminuire i valori di input", d1
 end if
end do
     
!print *, "Numero massimo di passi Lanczos (minore o uguale a d1*d2, consigliato nell'ordine di 100)="
!read *, lmax
lmax = 500
lmax = min(lmax, d1)

print*, "lx ", lx, "nu ", nu, "lmax", lmax, "dimh", d1

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

! recursive function binom(n,k) result (bun)
!     ! binomial coefficient n!/k!/(n-k)! using recursion
!     implicit none
!     integer :: bun, n, k
!     !
!     if ( n < 1 .or. k < 0 .or. k > n ) then
!        bun = 0
!     else if ( n == 1 ) then
!        bun = 1
!     else
!        bun = binom(n-1, k-1) + binom(n-1,k)
!     end if
!     !
! end function binom

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


  subroutine hoplist(d,n, conn, cin, ivp)
  implicit none
  integer(kind = di), intent(in) :: d
  integer, intent(in) :: n
  integer(kind = di), dimension(lx,lx,n), intent(in) :: ivp 
  integer :: l, k, j, icount, r, b
  integer(kind = di), dimension(d,2*n), intent(out) :: conn
  real(kind = dp), dimension(d,2*n), intent(out) :: cin
  integer, dimension(:), allocatable :: conf, nconf
  INTEGER :: neigh
  integer(kind = di) :: i , f
  logical :: uscita 
  integer, dimension(n) :: order 
  !print *, "hoplist start"
  l=lx
  allocate(conf(n), nconf(n))
  
    do i=1,d
      icount = 1
      call dritto(l,n,ivp,i,conf)
      !print*, "conf", i, conf
      do j=1,n
        do k = 0, 1
          call getNeighbours(conf(j), neigh, k)
          !print *, "wewe", conf(j), neigh, k
          if ( neigh == 0 ) cycle
          !print *, "neigh", conf(j), neigh 
            nconf = conf 
            nconf(j) = neigh
            uscita = .false.
            do b = 1, n
                if ( b == j ) cycle
                if ( conf(b) == nconf(j) ) then
                    uscita = .true.
                end if 
            end do
            if ( uscita ) cycle
            call rovescio(l,n,ivp,nconf,f, order, r)
            conn(i, icount) = f

            !print *, "nuun", i, f
            cin(i,icount) = - 0.5_dp*Jsign
            icount = icount + 1
        end do
          end do
    end do
    !print *, "hoplist end"

  deallocate(conf, nconf)
  end subroutine hoplist


  subroutine hamilton(connu, cinu, vint, ivp)
  implicit none
  integer(kind = di), dimension(lx,lx,nu), intent(inout) :: ivp 
  integer(kind = di), dimension(d1,2*nu), intent(out) :: connu
  real(kind = dp), dimension(d1,2*nu), intent(out):: cinu
  integer, dimension(d1), intent(out):: vint

    call nn_int(vint, ivp)
    call hoplist(d1,nu, connu, cinu, ivp)
  
  end subroutine hamilton
  

subroutine lanczos_alg(connu, cinu, v, t, w, alfa, beta, vint)
implicit none
  integer(kind = di), dimension(d1,2*nu), intent(in) :: connu
  real(kind = dp), dimension(d1,2*nu), intent(in):: cinu
  real(kind = dp), dimension(d1), intent(out) :: v, t, w
  real (kind = dp), dimension(lmax), intent(out) :: alfa, beta
  integer, dimension(d1), intent(out) :: vint

  

call rnd_cplx16_norm(d1,  v)

call lanczos_eigenvalue(v, connu, cinu, t, w, alfa, beta, vint)


end subroutine lanczos_alg

subroutine rnd_cplx16_norm(n,  v)
implicit none

       INTEGER(di), intent(in) :: n
       real(kind = dp), dimension(n), intent(out) ::   v
       REAL(dp) :: norm
       INTEGER(di) :: i
       integer :: z
       integer, dimension(:), allocatable :: seed 
       call random_seed(size = z)
       allocate(seed(z))
       call random_seed(get = seed)
       call random_number(v)
      norm = 0

             DO i = 1, n
              norm = norm + v(i)**2
             end do
           
             do i = 1, n
              v(i) = v(i)/sqrt(norm)
            end do

end subroutine rnd_cplx16_norm


subroutine hpsi(v1, v2, connu, cinu, vint)
implicit none
  integer(kind = di), dimension(d1,2*nu), intent(in) :: connu
  real(kind = dp), dimension(d1,2*nu), intent(in):: cinu
  real(kind = dp), dimension(d1) :: v1, v2
  integer(kind = di) :: i, m
  integer :: l
  integer, dimension(d1), intent(in) :: vint


v2 = 0
do i=1,d1
  v2(i) = v2(i) - Jsign*0.25_dp*vint(i)*v1(i)
end do
  do i = 1, d1
    do l = 1, 2*nu
      if ( connu(i,l) .ne. 0 ) then
        m = connu(i,l)
        v2(m) = v2(m) + cinu(i,l)*v1(i)
      end if     
    end do
  end do



end subroutine hpsi

subroutine lanczos_eigenvalue(v, connu, cinu, t, w, alfa, beta, vint)
implicit none

  integer(kind = di), dimension(d1,4*nu), intent(in) :: connu
  real(kind = dp), dimension(d1,4*nu), intent(in):: cinu
  real(kind = dp), dimension(d1), intent(inout) :: v
  real(kind = dp), dimension(d1), intent(out) :: t, w
  real (kind = dp), dimension(lmax), intent(out) :: alfa, beta
  real (kind = dp), dimension(:), allocatable ::alfa_tmp, beta_tmp
  integer(di) :: k, i
  real (kind = dp), dimension(:, :), allocatable :: eig
  real (kind = dp) :: etemp
  real(dp) ::  epsilon 
  real(kind = dp), dimension(d1) :: v_rnd
  integer, dimension(d1), intent(inout) :: vint


    beta = 0
    alfa = 0
    w = 0
    etemp = 0
    !epsilon = 0.000000000001
    epsilon = 0.000001
    v_rnd = v
  
    DO k=1, lmax
        call hpsi(v, t, connu, cinu, vint)
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
          print*, k, alfa_tmp(1)/lx
           if ( abs(alfa_tmp(1) - etemp) < epsilon) then
            print *, "gs totale, gs per sito", k, alfa_tmp(1)/(lx)
            !open(unit =40, file ="ene.dat", position = "append", status = "old", action ="write")
            !write(40, *) lx, alfa_tmp(1)/lx
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

subroutine nn_int(vint, ivp)
    implicit none
    integer(di) :: i
    integer :: j
    integer :: neigh
    integer, dimension(nu):: confu
    integer(kind = di), dimension(lx,lx,nu), intent(in) :: ivp
    integer, dimension(d1), intent(out) :: vint
    integer, dimension(lx) :: lat
    
    vint = 0
    
    do i=1,d1
        call dritto(lx,nu,ivp,i,confu)
        call smeno(confu, lat)
        do j = 1, lx 
          call getNeighbours(j, neigh, 1) 
          !print *, "op", j, neigh, pbc
          if ( neigh == 0 ) cycle
            if ( lat(j) == lat(neigh) ) vint(i) = vint(i) + 1 
            if ( lat(j) /= lat(neigh) ) vint(i) = vint(i) - 1 
        end do
    end do
    
    end subroutine nn_int

    subroutine getNeighbours(confuj, neigh, dir)
    implicit none
    integer,intent(in) :: confuj, dir
    integer,intent(out) ::  neigh
    integer :: i 

    i = 0
    !print *, "dir", dir 

if ( dir == 1 .and. pbc .eqv. .true.) then
  if ( confuj == lx  ) i=lx
  ! neigh(j) is the index of the neighbor of the j-th spin
  !dir = 1 avanti, dir = 0 dietro
  neigh=confuj+1-i
  !print *, "n1", confuj, i, dir, pbc, neigh
end if
if ( dir == 0 .and. pbc .eqv. .true.) then
  if ( confuj == 1  ) i=lx
  ! neigh(j) is the index of the neighbor of the j-th spin
  !dir = 1 avanti, dir = 0 dietro
  neigh=confuj-1+i
  if (neigh == (confuj+1) ) neigh = 0
  !print *, "n2", confuj, i, dir, pbc, neigh
end if

if ( pbc .neqv. .true. ) then
  if ( dir == 1  ) then
    if ( confuj == lx  ) i=lx+1
    ! neigh(j) is the index of the neighbor of the j-th spin
    !dir = 1 avanti, dir = 0 dietro
    neigh=confuj+1-i
    !print *, "n3", confuj, i, dir, pbc, neigh
  end if
end if 
if ( pbc .neqv. .true. ) then
  if ( dir == 0  ) then
    if ( confuj == 1  ) i=0
    ! neigh(j) is the index of the neighbor of the j-th spin
    !dir = 1 avanti, dir = 0 dietro
    neigh=confuj-1+i
    !print *, "n4", confuj, i, dir, pbc, neigh
  end if
end if 


    
    end subroutine getNeighbours

    subroutine smeno(conf,  lat)
    implicit none
    integer, dimension(nu),intent(in) :: conf
    integer, dimension(lx),intent(out) ::  lat
    integer :: i

      lat = -1
      do i = 1, nu
        lat(conf(i)) = 1
      end do
    
    end subroutine smeno


end module hilbert_space

program hilb
    
    use hilbert_space
    implicit none
    integer(kind = di), dimension(:,:), allocatable :: connu
    real(kind = dp), dimension(:,:), allocatable:: cinu
    integer, dimension(:), allocatable:: vint
    real(kind = dp), dimension(:), allocatable :: v, t, w
    real (kind = dp), dimension(:), allocatable :: alfa, beta
    integer :: ipmax, ill, imm, ip, ib
    integer(di) ::  bin
    INTEGER(di), dimension(:,:,:), allocatable :: ivpup
    
    call input()
    allocate(mu(d1)) 
    ALLOCATE(ivpup(lx,lx,nu))
    allocate(connu(d1,2*nu), cinu(d1,2*nu))
    allocate(alfa(lmax), beta(lmax), v(d1), t(d1), w(d1))
    allocate(vint(d1))
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
           !bin = binom(ill-ib, imm-1)
           CALL binom(ill-ib, imm-1, bin)
            ivpup(ip,ill,imm)=ivpup(ip,ill,imm)+bin
           enddo
          enddo
         endif
        enddo
       enddo
      connu = 0
      call hamilton(connu, cinu, vint, ivpup)
      call lanczos_alg(connu, cinu, v, t, w, alfa, beta, vint)
    
    deallocate(alfa, ivpup) 
    deallocate(beta, v, t, w) 
    deallocate(mu, vint)
    deallocate(connu)
    deallocate(cinu) 

  end program hilb