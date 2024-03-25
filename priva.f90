program cannon
    implicit none
    include "mpif.h"
    integer :: cannon_comm, status(mpi_status_size), error_number, rank, number_of_processes, nl
    real(8), allocatable, dimension (:,:) :: a, b
    real(8), allocatable, dimension (:,:) :: a1, a2, a3, a4, b1, b2, b3, b4
    real(8), allocatable, dimension (:,:) :: a11, a21, a31, a41, b11, b21, b31, b41
    real(8), allocatable, dimension (:,:) :: c1, c2, c3, c4, c, tmp, tmp1, tmp2, tmp3
    integer :: n = 8
    integer :: i, j, k, ordina
    integer :: start
    integer :: end, iseed
    integer :: recv_start
    real(8) :: iniz, fin, rnd, alfa =1, beta = 0
    
    
    call mpi_init(error_number)
    call mpi_comm_size(mpi_comm_world, number_of_processes, error_number)
    call mpi_comm_rank(mpi_comm_world, rank, error_number)

    nl = n/2
    allocate(a(n, n), b(n, n), c(n, n), a1(nl, nl), b1(nl, nl), c1(nl, nl))
    allocate(a2(nl, nl), b2(nl, nl), c2(nl, nl), a3(nl, nl), b3(nl, nl), c3(nl, nl))
    allocate(a4(nl, nl), b4(nl, nl), c4(nl, nl), tmp(nl,nl))

    if ( number_of_processes == 4 ) then
        iniz=MPI_Wtime()
    if ( rank == 0 ) then

    do i = 1, n
        do j = 1, n
            a(i,j) = real(10*j,8)
            b(i,j) = 0.01/real(i,8)
            c(i,j) = 0
        end do
    end do

    do ordina = 1, 3
     CALL MPI_send(a, n*n, mpi_real8, ordina, 0, mpi_comm_world, error_number)
     CALL MPI_send(b, n*n, mpi_real8, ordina, 0, mpi_comm_world, error_number)
    end do


        CALL MPI_send(c2, nl*nl, mpi_real8, 1, 0, mpi_comm_world, error_number)
        CALL MPI_send(c3, nl*nl, mpi_real8, 2, 0, mpi_comm_world, error_number)
        CALL MPI_send(c4, nl*nl, mpi_real8, 3, 0, mpi_comm_world, error_number)
        do i = 1, nl 
            do j = 1, nl 
                a1(i,j) = a(i,j)
                a2(i,j) = a(i,j+nl)
                a3(i,j) = a(i+nl,j)
                a4(i,j) = a(i+nl,j+nl)
                b1(i,j) = b(i,j)
                b2(i,j) = b(i,j+nl)
                b3(i,j) = b(i+nl,j)
                b4(i,j) = b(i+nl,j+nl)
                
            end do
            
        end do

        c1=0
            call dgemm('N', 'N', nl, nl, nl, alfa, a1, nl, b1, nl, beta, tmp, nl )
            c1 = c1 + tmp
            call dgemm('N', 'N', nl, nl, nl, alfa, a2, nl, b3, nl, beta, tmp, nl )
            c1 = c1 + tmp
        deallocate(tmp)
        CALL MPI_RECV(c2, nl*nl, mpi_real8, 1, 0, mpi_comm_world, status, error_number)
        CALL MPI_RECV(c3, nl*nl, mpi_real8, 2, 0, mpi_comm_world, status, error_number)
        CALL MPI_RECV(c4, nl*nl, mpi_real8, 3, 0, mpi_comm_world, status, error_number)


        do i = 1, nl 
            do j = 1, nl 
                c(i,j) = c1(i,j)
                c(i+nl,j) = c3(i,j)
                c(i,j+nl) = c2(i,j)
                c(i+nl,j+nl) = c4(i,j)
                
            end do
            
        end do
        open(unit = 12, file = "Matrice_C_parallelo.txt")
        write(12,*) "Matrix C"
        do i = 1, n
       
            write(12, *) (c(i, :))
           
        end do
        deallocate(a, b, c, a1, a2, a3, a4, b1, b2, b3, b4, c1,c2, c3, c4)





    else

        CALL MPI_RECV(a, n*n, mpi_real8, 0, 0, mpi_comm_world, status, error_number)
        CALL MPI_RECV(b, n*n, mpi_real8, 0, 0, mpi_comm_world, status, error_number)
        allocate(  a11(nl, nl), b11(nl, nl))
        allocate(a21(nl, nl), b21(nl, nl), a31(nl, nl), b31(nl, nl))
        allocate(a41(nl, nl), b41(nl, nl))
        do i = 1, nl 
            do j = 1, nl 

                a11(i,j) = a(i,j)
                a31(i,j) = a(i+nl,j)
                a21(i,j) = a(i,j+nl)
                a41(i,j) = a(i+nl,j+nl)
                b11(i,j) = b(i,j)
                b31(i,j) = b(i+nl,j)
                b21(i,j) = b(i,j+nl)
                b41(i,j) = b(i+nl,j+nl)
                
            end do
        enddo
        if ( rank == 1 ) then
            allocate(tmp1(nl,nl))
            CALL MPI_RECV(c2, nl*nl, mpi_real8, 0, 0, mpi_comm_world, status, error_number)
            c2=0
            call dgemm('N', 'N', nl, nl, nl, alfa, a11, nl, b21, nl, beta, tmp1, nl )
            c2 = c2 + tmp1
            call dgemm('N', 'N', nl, nl, nl, alfa, a21, nl, b41, nl, beta, tmp1, nl )
            c2 = c2 + tmp1

            CALL MPI_send(c2, nl*nl, mpi_real8, 0, 0, mpi_comm_world, error_number)

            deallocate(tmp1)
        end if
        
        if ( rank == 2 ) then
            allocate(tmp2(nl,nl))
        CALL MPI_RECV(c3, nl*nl, mpi_real8, 0, 0, mpi_comm_world, status, error_number)
        c3=0
        call dgemm('N', 'N', nl, nl, nl, alfa, a31, nl, b11, nl, beta, tmp2, nl )
        c3 = c3 + tmp2
        call dgemm('N', 'N', nl, nl, nl, alfa, a41, nl, b31, nl, beta, tmp2, nl )
        c3 = c3 + tmp2

        CALL MPI_send(c3, nl*nl, mpi_real8, 0, 0, mpi_comm_world, error_number)
        deallocate(tmp2)
        end if

        if ( rank == 3 ) then
            allocate(tmp3(nl,nl))
        CALL MPI_RECV(c4, nl*nl, mpi_real8, 0, 0, mpi_comm_world, status, error_number)
        c4=0
        call dgemm('N', 'N', nl, nl, nl, alfa, a31, nl, b21, nl, beta, tmp3, nl )
        c4 = c4 + tmp3
        call dgemm('N', 'N', nl, nl, nl, alfa, a41, nl, b41, nl, beta, tmp3, nl )
        c4 = c4 + tmp3
        CALL MPI_send(c4, nl*nl, mpi_real8, 0, 0, mpi_comm_world, error_number)
        deallocate(tmp3)

        end if

            
        deallocate(a1, a2, a3, a4, b1, b2, b3, b4, c2, c3, c4)
    end if

    
    fin = MPI_Wtime()
    if ( rank == 0 ) then
        open(unit=21, file = "Tempo_parallelo.txt")
        write(21,*)  fin - iniz
        close(21)
    end if


    call  MPI_Finalize(error_number)
    end if
    if ( number_of_processes == 3 ) then
        iniz=MPI_Wtime()
    if ( rank == 0 ) then

    do i = 1, n
        do j = 1, n
            a(i,j) = real(10*j,8)
            b(i,j) = 0.01/real(i,8)
            c(i,j) = 0
        end do
    end do

    do ordina = 1, 2
     CALL MPI_send(a, n*n, mpi_real8, ordina, 0, mpi_comm_world, error_number)
     CALL MPI_send(b, n*n, mpi_real8, ordina, 0, mpi_comm_world, error_number)
    end do


        CALL MPI_send(c2, nl*nl, mpi_real8, 1, 0, mpi_comm_world, error_number)
        CALL MPI_send(c3, nl*nl, mpi_real8, 2, 0, mpi_comm_world, error_number)
        CALL MPI_send(c4, nl*nl, mpi_real8, 2, 0, mpi_comm_world, error_number)
        do i = 1, nl 
            do j = 1, nl 
                a1(i,j) = a(i,j)
                a2(i,j) = a(i,j+nl)
                a3(i,j) = a(i+nl,j)
                a4(i,j) = a(i+nl,j+nl)
                b1(i,j) = b(i,j)
                b2(i,j) = b(i,j+nl)
                b3(i,j) = b(i+nl,j)
                b4(i,j) = b(i+nl,j+nl)
                
            end do
            
        end do


        c1=0
            call dgemm('N', 'N', nl, nl, nl, alfa, a1, nl, b1, nl, beta, tmp, nl )
            c1 = c1 + tmp
            call dgemm('N', 'N', nl, nl, nl, alfa, a2, nl, b3, nl, beta, tmp, nl )
            c1 = c1 + tmp
        deallocate(tmp)
        CALL MPI_RECV(c2, nl*nl, mpi_real8, 1, 0, mpi_comm_world, status, error_number)
        CALL MPI_RECV(c3, nl*nl, mpi_real8, 2, 0, mpi_comm_world, status, error_number)
        CALL MPI_RECV(c4, nl*nl, mpi_real8, 2, 0, mpi_comm_world, status, error_number)


        do i = 1, nl 
            do j = 1, nl 
                c(i,j) = c1(i,j)
                c(i+nl,j) = c3(i,j)
                c(i,j+nl) = c2(i,j)
                c(i+nl,j+nl) = c4(i,j)
                
            end do
            
        end do
        open(unit = 12, file = "Matrice_C_parallelo.txt")
        write(12,*) "Matrix C"
        do i = 1, n
       
            write(12, *) (c(i, :))
           
        end do
        deallocate(a, b, c, a1, a2, a3, a4, b1, b2, b3, b4, c1,c2, c3, c4)





    else

        CALL MPI_RECV(a, n*n, mpi_real8, 0, 0, mpi_comm_world, status, error_number)
        CALL MPI_RECV(b, n*n, mpi_real8, 0, 0, mpi_comm_world, status, error_number)
        allocate(  a11(nl, nl), b11(nl, nl))
        allocate(a21(nl, nl), b21(nl, nl), a31(nl, nl), b31(nl, nl))
        allocate(a41(nl, nl), b41(nl, nl))
        do i = 1, nl 
            do j = 1, nl 

                a11(i,j) = a(i,j)
                a31(i,j) = a(i+nl,j)
                a21(i,j) = a(i,j+nl)
                a41(i,j) = a(i+nl,j+nl)
                b11(i,j) = b(i,j)
                b31(i,j) = b(i+nl,j)
                b21(i,j) = b(i,j+nl)
                b41(i,j) = b(i+nl,j+nl)
                
            end do
        enddo
        if ( rank == 1 ) then
            allocate(tmp1(nl,nl))
            CALL MPI_RECV(c2, nl*nl, mpi_real8, 0, 0, mpi_comm_world, status, error_number)
            c2=0
            call dgemm('N', 'N', nl, nl, nl, alfa, a11, nl, b21, nl, beta, tmp1, nl )
            c2 = c2 + tmp1
            call dgemm('N', 'N', nl, nl, nl, alfa, a21, nl, b41, nl, beta, tmp1, nl )
            c2 = c2 + tmp1

            CALL MPI_send(c2, nl*nl, mpi_real8, 0, 0, mpi_comm_world, error_number)

            deallocate(tmp1)
        end if
        
        if ( rank == 2 ) then
            allocate(tmp2(nl,nl))
        CALL MPI_RECV(c3, nl*nl, mpi_real8, 0, 0, mpi_comm_world, status, error_number)

        CALL MPI_RECV(c4, nl*nl, mpi_real8, 0, 0, mpi_comm_world, status, error_number)
        c3=0
        call dgemm('N', 'N', nl, nl, nl, alfa, a31, nl, b11, nl, beta, tmp2, nl )
        c3 = c3 + tmp2
        call dgemm('N', 'N', nl, nl, nl, alfa, a41, nl, b31, nl, beta, tmp2, nl )
        c3 = c3 + tmp2


        c4=0
        call dgemm('N', 'N', nl, nl, nl, alfa, a31, nl, b21, nl, beta, tmp2, nl )
        c4 = c4 + tmp2
        call dgemm('N', 'N', nl, nl, nl, alfa, a41, nl, b41, nl, beta, tmp2, nl )
        c4 = c4 + tmp2
        DEALLOCATE(tmp2)
        CALL MPI_send(c4, nl*nl, mpi_real8, 0, 0, mpi_comm_world, error_number)

        CALL MPI_send(c3, nl*nl, mpi_real8, 0, 0, mpi_comm_world, error_number)
        end if

            
        deallocate(a1, a2, a3, a4, b1, b2, b3, b4, c2, c3, c4)
    end if

    
    fin = MPI_Wtime()
    if ( rank == 0 ) then
        open(unit=21, file = "Tempo_parallelo.txt")
        write(21,*)  fin - iniz
        close(21)
    end if


    call  MPI_Finalize(error_number)

    end if

    if ( number_of_processes == 2 ) then
        iniz=MPI_Wtime()
    if ( rank == 0 ) then

    do i = 1, n
        do j = 1, n
            a(i,j) = real(10*j,8)
            b(i,j) = 0.01/real(i,8)
            c(i,j) = 0
        end do
    end do

     CALL MPI_send(a, n*n, mpi_real8, 1, 0, mpi_comm_world, error_number)
     CALL MPI_send(b, n*n, mpi_real8, 1, 0, mpi_comm_world, error_number)

        CALL MPI_send(c3, nl*nl, mpi_real8, 1, 0, mpi_comm_world, error_number)
        CALL MPI_send(c4, nl*nl, mpi_real8, 1, 0, mpi_comm_world, error_number)
        do i = 1, nl 
            do j = 1, nl 
                a1(i,j) = a(i,j)
                a2(i,j) = a(i,j+nl)
                a3(i,j) = a(i+nl,j)
                a4(i,j) = a(i+nl,j+nl)
                b1(i,j) = b(i,j)
                b2(i,j) = b(i,j+nl)
                b3(i,j) = b(i+nl,j)
                b4(i,j) = b(i+nl,j+nl)
                
            end do
            
        end do


        c1=0
            call dgemm('N', 'N', nl, nl, nl, alfa, a1, nl, b1, nl, beta, tmp, nl )
            c1 = c1 + tmp
            call dgemm('N', 'N', nl, nl, nl, alfa, a2, nl, b3, nl, beta, tmp, nl )
            c1 = c1 + tmp

            c2=0
            call dgemm('N', 'N', nl, nl, nl, alfa, a1, nl, b2, nl, beta, tmp, nl )
            c2 = c2 + tmp
            call dgemm('N', 'N', nl, nl, nl, alfa, a2, nl, b4, nl, beta, tmp, nl )
            c2 = c2 + tmp


        deallocate(tmp)
        CALL MPI_RECV(c3, nl*nl, mpi_real8, 1, 0, mpi_comm_world, status, error_number)
        CALL MPI_RECV(c4, nl*nl, mpi_real8, 1, 0, mpi_comm_world, status, error_number)


        do i = 1, nl 
            do j = 1, nl 
                c(i,j) = c1(i,j)
                c(i+nl,j) = c3(i,j)
                c(i,j+nl) = c2(i,j)
                c(i+nl,j+nl) = c4(i,j)
                
            end do
            
        end do
        open(unit = 12, file = "Matrice_C_parallelo.txt")
        write(12,*) "Matrix C"
        do i = 1, n
       
            write(12, *) (c(i, :))
           
        end do

        deallocate(a, b, c, a1, a2, a3, a4, b1, b2, b3, b4, c1,c2, c3, c4)





    else

        CALL MPI_RECV(a, n*n, mpi_real8, 0, 0, mpi_comm_world, status, error_number)
        CALL MPI_RECV(b, n*n, mpi_real8, 0, 0, mpi_comm_world, status, error_number)
        allocate( b11(nl, nl))
        allocate(b21(nl, nl), a31(nl, nl), b31(nl, nl))
        allocate(a41(nl, nl), b41(nl, nl))
        do i = 1, nl 
            do j = 1, nl 

                a31(i,j) = a(i+nl,j)
                a41(i,j) = a(i+nl,j+nl)
                b11(i,j) = b(i,j)
                b31(i,j) = b(i+nl,j)
                b21(i,j) = b(i,j+nl)
                b41(i,j) = b(i+nl,j+nl)
                
            end do
        enddo

 
            allocate(tmp2(nl,nl))
        CALL MPI_RECV(c3, nl*nl, mpi_real8, 0, 0, mpi_comm_world, status, error_number)

        CALL MPI_RECV(c4, nl*nl, mpi_real8, 0, 0, mpi_comm_world, status, error_number)
        c3=0
        call dgemm('N', 'N', nl, nl, nl, alfa, a31, nl, b11, nl, beta, tmp2, nl )
        c3 = c3 + tmp2
        call dgemm('N', 'N', nl, nl, nl, alfa, a41, nl, b31, nl, beta, tmp2, nl )
        c3 = c3 + tmp2


        c4=0
        call dgemm('N', 'N', nl, nl, nl, alfa, a31, nl, b21, nl, beta, tmp2, nl )
        c4 = c4 + tmp2
        call dgemm('N', 'N', nl, nl, nl, alfa, a41, nl, b41, nl, beta, tmp2, nl )
        c4 = c4 + tmp2
        DEALLOCATE(tmp2)
        CALL MPI_send(c4, nl*nl, mpi_real8, 0, 0, mpi_comm_world, error_number)

        CALL MPI_send(c3, nl*nl, mpi_real8, 0, 0, mpi_comm_world, error_number)



            
        deallocate( a31, a41, b11, b21, b31, b41, c3, c4)
    end if

    
    fin = MPI_Wtime()
    if ( rank == 0 ) then
        open(unit=21, file = "Tempo_parallelo.txt")
        write(21,*)  fin - iniz
        close(21)
    end if


    call  MPI_Finalize(error_number)
    end if


end program cannon