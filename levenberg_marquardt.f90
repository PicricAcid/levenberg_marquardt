program levenberg_marquardt
  implicit none
!Parameters
  integer, parameter:: n = 2
  integer, parameter:: m = 7
  integer, parameter:: itaration = 10
!variables
  double precision:: rcond
  integer:: i, j, ifail, info, lda, ldb, nrhs, ita
  double precision:: beta1
  double precision:: beta2
  double precision:: snew, sold
  double precision:: lambda = 0.0001
!arrays
  double precision, allocatable :: a(:, :), b(:)
  integer, allocatable :: ipiv(:)
  double precision, allocatable:: jacobian(:, :), e(:)
  real:: x0(7) = (/0.038, 0.194, 0.425, 0.626, 1.253, 2.500, 3.740/)
  real:: y0(7) = (/0.050, 0.127, 0.094, 0.2122, 0.2729, 0.2665, 0.3317/)
!init variables
  write(*, *) 'init beta value'
  read(*, *) beta1, beta2
  nrhs = 1
  lda = n
  ldb = m
  allocate (a(lda, n), b(ldb), ipiv(n))
  allocate(jacobian(n, m), e(m))
  do i=1, m
    snew = snew + (y0(i) - ((beta1*x0(i))/(beta2 + x0(i))))**2
  end do

!main loop
  do ita=1, itaration
    write(*, *) 'itaration =', ita, '/', itaration

!update matrix
    jacobian = 0
    sold = snew
    snew = 0
    a = 0
    b = 0

    do i=1, m
      jacobian(1, i) = -x0(i)/(beta2 + x0(i))
      jacobian(2, i) = (beta1*x0(i))/((beta2 + x0(i))**2)
      e(i) = y0(i) - ((beta1*x0(i))/(beta2 + x0(i)))
    end do

    a = matmul(jacobian, transpose(jacobian))

    do i=1, lda
      do j=1, n
        if(i==j) then
          a(i, j) = a(i, j) + lambda
        else
          a(i, j) = a(i, j)
        end if
      end do
    end do

    b = -matmul(e, transpose(jacobian))

    do i=1, m
      snew = snew + (y0(i) - ((beta1*x0(i))/(beta2 + x0(i))))**2
    end do

!solve the equations Ax = b for x
    call dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)

    if (info == 0) then
      write(*,*) 'DGESV works fine'
    else
        write(*,*) 'INFO = ', info
        write(*,*) 'Something is wrong'
    end if
!update beta
    beta1 = beta1 + b(1)
    beta2 = beta2 + b(2)

    if(snew >= sold) then
      lambda = lambda*10
    else
      lambda = lambda/10
    end if

    write(*, *) 'S = ', snew
    write(*, *) '----------------------------------------'

  end do

  write(*, *) beta1
  write(*, *) beta2

end program
