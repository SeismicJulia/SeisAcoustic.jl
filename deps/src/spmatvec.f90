! !--------------------------------------------------------------------
! !     Float32: y = A * x
! !--------------------------------------------------------------------
subroutine a_mul_b_ss(m, n, A, jA, iA, x, y)

!DIR$ ATTRIBUTES DLLEXPORT :: a_mul_b_rr
!DIR$ ATTRIBUTES ALIAS: 'a_mul_b_rr_':: a_mul_b_rr

implicit none

integer(kind=8),intent(in)   :: m  ! # of rows in A
integer(kind=8),intent(in)   :: n  ! # of columns in A

real(kind=4)   ,intent(in)   :: A(*)
integer(kind=8),intent(in)   :: jA(*), iA(n+1)
real(kind=4)   ,intent(inout):: y(m)
real(kind=4)   ,intent(in)   :: x(n)

integer i, j, j1, j2, jaj

y = 0.0
include "a_mul_b.fi"

return
end subroutine a_mul_b_ss

! !--------------------------------------------------------------------
! !     Float64: y = A * x
! !--------------------------------------------------------------------
subroutine a_mul_b_rr(m, n, A, jA, iA, x, y)

!DIR$ ATTRIBUTES DLLEXPORT :: a_mul_b_rr
!DIR$ ATTRIBUTES ALIAS: 'a_mul_b_rr_':: a_mul_b_rr

implicit none

integer(kind=8),intent(in)   :: m  ! # of rows in A
integer(kind=8),intent(in)   :: n  ! # of columns in A
real(kind=8)   ,intent(in)   :: A(*)
integer(kind=8),intent(in)   :: jA(*), iA(n+1)
real(kind=8)   ,intent(inout):: y(m)
real(kind=8)   ,intent(in)   :: x(n)

integer i, j, j1, j2, jaj

y = 0.d0
include "a_mul_b.fi"

return
end subroutine a_mul_b_rr


! !--------------------------------------------------------------------
! !     Float32: y = A' * x
! !--------------------------------------------------------------------
subroutine ac_mul_b_ss(m, n, A, jA, iA, x, y)
!DIR$ ATTRIBUTES DLLEXPORT :: ac_mul_b_rr
!DIR$ ATTRIBUTES ALIAS: 'ac_mul_b_rr_':: ac_mul_b_rr

implicit none

integer(kind=8),intent(in)   :: m  ! # of rows in A
integer(kind=8),intent(in)   :: n  ! # of columns in A

real(kind=4)   ,intent(in)   :: A(*)
integer(kind=8),intent(in)   :: jA(*), iA(n+1)
real(kind=4)   ,intent(in)   :: x(m)
real(kind=4)   ,intent(inout):: y(n)

integer(kind=4) i, j, j1, j2


y = 0.0
#include "ac_mul_b.fi"

return

end subroutine ac_mul_b_ss

! !--------------------------------------------------------------------
! !     Float64: y = A' * x
! !--------------------------------------------------------------------
subroutine ac_mul_b_rr(m, n, A, jA, iA, x, y)
!DIR$ ATTRIBUTES DLLEXPORT :: ac_mul_b_rr
!DIR$ ATTRIBUTES ALIAS: 'ac_mul_b_rr_':: ac_mul_b_rr

implicit none

integer(kind=8),intent(in)   :: m  ! # of rows in A
integer(kind=8),intent(in)   :: n  ! # of columns in A

real(kind=8)   ,intent(in)   :: A(*)
integer(kind=8),intent(in)   :: jA(*), iA(n+1)
real(kind=8)   ,intent(in)   :: x(m)
real(kind=8)   ,intent(inout):: y(n)

integer i, j, j1, j2

! make a vector to be zero
y = 0.d0

#include "ac_mul_b.fi"

return
end subroutine ac_mul_b_rr
