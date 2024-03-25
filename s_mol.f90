program serial
    implicit none

    integer, parameter :: dp = selected_real_kind(8)
    real(dp), dimension(:,:), allocatable :: m_a, m_b, m_c
    real(dp) :: alfa, beta, start, finish
    integer(8) :: r_a, c_a, c_b, i, j, n

    n=8
    r_a=n
    c_a=n
    c_b=n
    OPEN(unit=10, file="Matrice_A.txt")
    OPEN(unit=11, file="Matrice_B.txt")
    OPEN(unit=12, file="Matrice_C_seriale.txt")
    OPEN(unit=13, file="Tempo_seriale.txt")

    call cpu_time(start)

    allocate(m_a(r_a, c_a), m_b(c_a, c_b), m_c(r_a, c_b))
    
    do i = 1, r_a
        do j = 1, c_a
            m_a(i, j) = real(10*j,dp)
        end do
    end do
    do i = 1, c_a
        do j = 1, c_b
            m_b(i, j) = 0.01/real(i,dp)
        end do
    end do

    write(10, *) "Matrix A"

do i = 1, r_a 

    write(10, *) m_a(i, :)
    
end do

write(11, *) "Matrix B"

do i = 1, c_a 

    write(11, *) m_b(i, :)
    
end do

alfa = 1
beta = 0

    call dgemm('N', 'N', r_a, c_b, c_a, alfa, m_a, r_a, m_b, c_b, beta, m_c, r_a )

    call cpu_time(finish)
    write(13,*) finish - start

 write(12,*) "Matrix C"
 do i = 1, r_a

     write(12, *) (m_c(i, :))
    
 end do

 close(10)
 close(11)
 close(12)
 close(13)

 

end program serial