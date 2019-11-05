module mathf
    implicit none
contains
    ! Integration of function f using Trapezoidal rule
    ! f   : Integrand array
    ! N   : Dimension of space excluding the first element
    ! dh  : Step distance of space
    ! sum : The result of the integration
    subroutine integrate(f, N, dh, sum)
        integer,intent(in)           :: N
        double precision,intent(in)  :: f(0:N), dh
        double precision,intent(out) :: sum
        integer i
        sum = 0d0
        do i = 0, N
            if (i == 0 .or. i == N) then
                sum = sum + 0.5*f(i)*dh
            else
                sum = sum + f(i)*dh
            end if
        end do
    end subroutine

    ! Normalize the given function f
    ! f   : Function to be normalized
    ! N   : Dimension of f
    ! dh  : step distance of space
    subroutine normalize(f, N, dh)
        integer,intent(in)                :: N
        double precision,intent(in)       :: dh
        complex(kind(0d0)),intent(inout)  :: f(0:N)
        double precision sum
        call Integrate(abs(f(:))**2d0, N, dh, sum)
        f(:) = f(:) / sqrt(sum)
    end subroutine normalize

    ! Solve Ax = lambda*x for lambda and x
    subroutine solve_eigen(A, z, lambda, dim)
        integer,intent(in)               :: dim                     ! ubound of matrix A
        double precision,intent(out)     :: lambda(0:dim)           ! Eigenvalue
        complex(kind(0d0)),intent(out)   :: z(0:dim, 0:dim)         ! Eigenvectors
        double precision,intent(in)      :: A(0:dim, 0:dim)         ! Matrix
        double precision                 :: d(0:dim), e(0:dim-1)    ! d:Diagonal/e:Non-diagonal elements
        double precision                 :: abstol                  ! Obsolete feature of LAPACK
        integer                          :: m                       ! Total number of eigenvalues found
        integer                          :: i                       ! Loop variable
        integer                          :: isuppz(0:2*dim+1)       ! The indices indicating the nonzero elements in Z
        double precision,allocatable     :: work(:)                 ! Workspace
        integer                          :: lwork                   ! The dimension of the array work
        integer,allocatable              :: iwork(:)                ! Workspace
        integer                          :: liwork                  ! The dimension of the array iwork
        integer                          :: info                    ! Success/Error indicator

        ! Substitute diagonal/non-diagonal elements of input matrix A
        do i = 0, dim
            d(i) = A(i, i)
            if (i < dim) then
                e(i) = A(i, i+1)
            end if
        end do
        
        allocate (work(0:dim), iwork(0:dim))
        call zstegr('V', 'I', dim+1, d, e, 1, 1, abstol, m, lambda, z, dim+1, isuppz, work, -1, iwork, -1, info)
        stop
        lwork = int(work(1))
        liwork = int(work(1))
        deallocate(work, iwork)
        allocate(work(lwork), iwork(liwork))
        call zstegr('V', 'I', dim+1, d, e, 1, 1, abstol, m, lambda, z, dim+1, isuppz, work, lwork, iwork, liwork, info)
        deallocate(work, iwork)
        if (info < 0) then
            write (*, *) "ERROR : Argument ", -info, " has illegal value"
        else if (info == 1) then
            write (*, *) "ERROR : The dqds algorithm failed to converge"
            stop
        else if (info == 2) then
            write (*, *) "ERROR : Inverse iteration failed to converge"
            stop
        end if
    end subroutine solve_eigen
end module mathf