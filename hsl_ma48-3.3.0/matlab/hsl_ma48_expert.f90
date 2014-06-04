!
! This routine offers direct access to hsl_ma48 functionality
!     [handle, info] = ma48_expert('factor', A[, control, P, Q])
!     [x, info] = ma48_expert('solve', handle, b[, control])
!     [x, info, handle] = ma48_expert('backslash', A, b[, control, P, Q])
!     ma48_expert('destroy', handle)
! P and Q are optional row and column orderings. Either both or neither should
! be present.
!

! Structure of code:
! The module ma48_matlab_main provides a series of functions, one for each
!    possible 'action'.
! The module ma48_handles provides a series of functions for dealing with
!    integer handles associated on the Fortran side with factors and order
! The mexFunction unwraps handles and calls the corrected subroutine
!    for each 'action'.

module ma48_matlab_main
!$ use omp_lib
   use hsl_matlab
   use hsl_ma48_double, only: ma48_factors,                    &
                              ma48_control,                    &
                              ma48_ainfo,                      &
                              ma48_finfo,                      &
                              ma48_sinfo,                      &
                              ma48_initialize,                 &
                              ma48_analyse,                    &
                              ma48_factorize,                  &
                              ma48_solve,                      &
                              ma48_finalize
   use hsl_zd11_double, only: zd11_type
   implicit none

   integer, parameter :: wp = kind(0d0)
   integer, parameter :: long = selected_int_kind(18)

contains
   subroutine copy_control_in(pm, control)
      integer(mwp_) :: pm
      type(ma48_control), intent(inout) :: control

      integer(mwp_) :: pc
      integer(int4_) :: fnum
      character(80) :: fname
      character(200) :: warnmsg

      do fnum = 1, MATLAB_get_no_fields(pm)
         fname = MATLAB_get_field_name_by_no(pm, fnum)
         select case(trim(fname))
         case("btf")
            call MATLAB_get_value(pm, fname, pc, control%btf)
         case("drop")
            call MATLAB_get_value(pm, fname, pc, control%drop)
         case("factor_blocking")
            call MATLAB_get_value(pm, fname, pc, control%factor_blocking)
         case("maxit")
            call MATLAB_get_value(pm, fname, pc, control%maxit)
         case("switch")
            call MATLAB_get_value(pm, fname, pc, control%switch)
         case("u")
            call MATLAB_get_value(pm, fname, pc, control%u)
         case default
            write(warnmsg, "(3a)") "Ignored unrecognised control parameter '", &
               trim(fname), "'"
            call MATLAB_warning(warnmsg)
         end select
      end do
   end subroutine copy_control_in

   subroutine unknownError(routine, flag)
      character(len=*) :: routine
      integer :: flag
      character(len=200) :: errmsg

      write(errmsg, "(3a,i3)") "Error return from ", routine, ". flag = ", flag
      call MATLAB_error(errmsg)
   end subroutine unknownError

   ! [handle, info] = ma48_expert('factor', A[, P])
   ! handle and 'factor' already dealt with
   subroutine ma48_matlab_analyse_factor(nlhs_in, plhs, nrhs_in, prhs, &
         matrix, factors)
      integer*4 :: nlhs_in, nrhs_in
      integer(mwp_) :: plhs(*), prhs(*)
      type(zd11_type), intent(inout) :: matrix
      type(ma48_factors), intent(inout) :: factors

      integer, parameter :: A_in = 1, &
                            control_in = 2, &
                            P_in = 3, &
                            Q_in = 4
      integer, parameter :: info_out = 1

      integer(mws_) :: mwm, mwn, temp
      integer :: i, flag, st

      type(ma48_control) :: control
      type(ma48_ainfo) :: ainfo
      type(ma48_finfo) :: finfo

      character(len=2000) :: errmsg

      real(wp) :: atime, ftime
      integer :: t_start, t_stop, t_rate
      integer, dimension(:), allocatable :: perm, row_perm, col_perm

      ! Check number of arguments
      if(nrhs_in.lt.1 .or. nrhs_in .gt.4) &
         call MATLAB_error("Wrong number of input arguments")
      if(nlhs_in.ne.0 .and. nlhs_in.ne.1) &
         call MATLAB_error("Wrong number of output arguments")

      ! Check matrix is sparse
      call matlab_to_fortran(prhs(A_in), mwm, mwn, matrix%ptr, matrix%row, &
         matrix%val, 'A')
      if(mwm.ne.mwn) call MATLAB_error("The matrix must be square")
      matrix%m = mwn
      matrix%n = mwn
      matrix%ne = matrix%ptr(matrix%n+1)-1

      ! Convert from CSC to COORD
      allocate(matrix%col(matrix%ne), stat=st)
      if(st.ne.0) goto 100
      do i = 1, matrix%n
         matrix%col(matrix%ptr(i):matrix%ptr(i+1)-1) = i
      end do
      deallocate(matrix%ptr, stat=st)

      ! Initalise
      call ma48_initialize(factors=factors, control=control)

      ! Read control in
      if(nrhs_in.ge.control_in) then
         ! User supplied control structure
         if(MATLAB_is_structure(prhs(control_in))) &
            call copy_control_in(prhs(control_in), control)
      endif

      if(nrhs_in.ge.P_in) then
         if(nrhs_in.eq.P_in) then
            call MATLAB_error("If P is present, Q must also be present.")
         else
            call matlab_to_fortran(prhs(P_in), row_perm, temp, 'P')
            if(mwm .ne. temp) then
               write(errmsg, "(4(a,i5))") &
                  "Dimensions of A and P inconsistent: A=", &
                  matrix%m, "x", matrix%n, ", size(P)=", temp
               call MATLAB_error(errmsg)
            endif
            call matlab_to_fortran(prhs(Q_in), col_perm, temp, 'Q')
            if(mwm .ne. temp) then
               write(errmsg, "(4(a,i5))") &
                  "Dimensions of A and Q inconsistent: A=", &
                  matrix%m, "x", matrix%n, ", size(Q)=", temp
               call MATLAB_error(errmsg)
            endif
         endif
         allocate(perm(matrix%m+matrix%n), stat=st)
         if(st.ne.0) goto 100
         perm(1:matrix%m) = row_perm(:)
         perm(matrix%m+1:matrix%m+matrix%n) = col_perm(:)
         deallocate(row_perm, stat=st)
         deallocate(col_perm, stat=st)
      endif

      ! Call analyse
      call system_clock(t_start)
      if(allocated(row_perm)) then
         call ma48_analyse(matrix, factors, control, ainfo)
      else
         call ma48_analyse(matrix, factors, control, ainfo, perm=perm)
      endif
      call system_clock(t_stop, t_rate)
      select case (ainfo%flag)
      case(0)
         ! Success. Do nothing
      case(-6)
         ! Problem with order
         if(allocated(perm)) then
            write(errmsg, *) &
               "Problem with P or Q."
         else
            write(errmsg, "(a, 10i4)") "Unkown bug with P and Q. Please report."
         endif
         call MATLAB_error(errmsg)
      case(4)
         call MATLAB_warning("Matrix is rank deficient")
      case default
         call unknownError("ma48_analyse", ainfo%flag)
      end select
      atime = real(t_stop-t_start)/t_rate

      ! Call factor
      call system_clock(t_start)
      call ma48_factorize(matrix, factors, control, finfo)
      call system_clock(t_stop, t_rate)
      select case (finfo%flag)
      case(0:1)
         ! Success. Do nothing.
      case(4)
         call MATLAB_warning("Matrix is rank deficient")
      case default
         call unknownError("ma48_factorize", finfo%flag)
      end select
      ftime = real(t_stop-t_start)/t_rate

      ! Setup info structure if required
      if(nlhs_in.ge.info_out) then
         call copy_info_out(plhs(info_out), ainfo, finfo, atime, ftime)
      endif

      return
      100 continue
      call MATLAB_error("Insufficient memory")

   contains
      subroutine copy_info_out(ml_info, ainfo, finfo, atime, ftime)
         integer(mwp_), intent(inout) :: ml_info
         type(ma48_ainfo), intent(in) :: ainfo
         type(ma48_finfo), intent(in) :: finfo
         real(wp), intent(in) :: atime, ftime

         character(len=50), dimension(7) :: fields = (/ &
            "drop        ", &
            "ops         ", &
            "rank        ", &
            "size_factor ", &
            "struc_rank  ", &
            "analyse_time", &
            "factor_time " /)

         ml_info = MATLAB_create_structure(fields)

         call MATLAB_set_field(ml_info, 'drop', finfo%drop)
         call MATLAB_set_field(ml_info, 'ops', finfo%ops)
         call MATLAB_set_field(ml_info, 'rank', finfo%rank)
         call MATLAB_set_field(ml_info, 'size_factor', finfo%size_factor)
         call MATLAB_set_field(ml_info, 'struc_rank', ainfo%struc_rank)
         call MATLAB_set_field(ml_info, 'analyse_time', atime)
         call MATLAB_set_field(ml_info, 'factor_time', ftime)
      end subroutine copy_info_out

   end subroutine ma48_matlab_analyse_factor

   ! [x, info] = ma48_expert('solve', handle, b)
   ! 'solve' and handle already dealt with
   subroutine ma48_matlab_solve(nlhs_in, plhs, nrhs_in, prhs, matrix, factors)
      integer*4 :: nlhs_in, nrhs_in
      integer(mwp_) :: plhs(*), prhs(*)
      type(zd11_type), intent(inout) :: matrix
      type(ma48_factors), intent(inout) :: factors

      integer, parameter :: b_in = 1, &
                            control_in = 2
      integer, parameter :: x_out = 1, &
                            info_out = 2

      integer(mws_) :: mwn, mwnrhs
      real(wp), dimension(:,:), allocatable :: rhs, x
      integer :: st

      type(ma48_control) :: control
      type(ma48_sinfo) :: info

      character(len=200) :: errmsg
      integer :: nrhs, i

      real(wp) :: stime
      integer :: t_start, t_stop, t_rate

      ! Check number of arguments
      if(nrhs_in.lt.1) call MATLAB_error("Insufficient input arguments")
      if(nrhs_in.gt.2) call MATLAB_error("Too many input arguments")
      if(nlhs_in.lt.1) call MATLAB_error("Insufficient output arguments")
      if(nlhs_in.gt.2) call MATLAB_error("Too many output arguments")

      ! Get rhs
      call matlab_to_fortran(prhs(b_in), rhs, mwn, mwnrhs, 'b')
      if(mwn .ne. matrix%m) then
         write(errmsg, "(4(a,i5))") "Dimensions of A and b inconsistent: A=", &
            matrix%m, "x", matrix%n, ", b=", mwn, "x", mwnrhs
         call MATLAB_error(errmsg)
      endif
      nrhs = mwnrhs

      ! Copy control in
      call ma48_initialize(control=control)
      if(nrhs_in.ge.control_in) then
         ! User supplied control structure
         if(MATLAB_is_structure(prhs(control_in))) &
            call copy_control_in(prhs(control_in), control)
      endif
      
      ! Call solve
      allocate(x(matrix%n,nrhs), stat=st)
      if(st.ne.0) goto 100
      call system_clock(t_start)
      do i = 1, nrhs
         call ma48_solve(matrix, factors, rhs(:, i), x(:, i), control, info)
      end do
      call system_clock(t_stop, t_rate)
      select case (info%flag)
      case(0)
         ! Success. Do nothing.
      case default
         call unknownError("ma48_solve", info%flag)
      end select
      stime = real(t_stop-t_start)/t_rate

      ! Setup info structure if required
      if(nlhs_in.ge.info_out) then
         call copy_info_out(plhs(info_out), info, stime)
      endif

      ! Copy solution out
      plhs(x_out) = fortran_to_matlab(x)

      return
      100 continue
      call MATLAB_error("Insufficient memory")

   contains
      subroutine copy_info_out(ml_info, sinfo, stime)
         integer(mwp_), intent(inout) :: ml_info
         type(ma48_sinfo), intent(in) :: sinfo
         real(wp), intent(in) :: stime

         character(len=50), dimension(1) :: fields = (/ &
            "solve_time    " /)

         ml_info = MATLAB_create_structure(fields)

         call MATLAB_set_field(ml_info, 'solve_time', stime)
      end subroutine copy_info_out

   end subroutine ma48_matlab_solve

   ! [x, info, handle] = ma48_expert('backslash', A, b[, control, P])
   ! handle and 'backslash' already dealt with

   ! ma48_expert('destroy', handle)
   ! arguments 'destroy' and handle already removed
   subroutine ma48_matlab_finalise(factors)
      type(ma48_factors), intent(inout) :: factors

      type(ma48_control) :: control
      integer :: info

      call ma48_initialize(control=control)
      call ma48_finalize(factors, control, info)
   end subroutine ma48_matlab_finalise

end module ma48_matlab_main

! This module looks after a SAVEd set of variables mapping integer handles
! to Fortran factors and order variables
module ma48_handles
   use hsl_ma48_double, only: ma48_factors
   use hsl_zd11_double, only: zd11_type
   use ma48_matlab_main
   implicit none

   ! Data associated with the handle
   ! Considered to be empty if associated(factors) is .false.
   type ma48_hdl
      type(zd11_type), pointer :: matrix => null()
      type(ma48_factors), pointer :: factors => null()
   end type ma48_hdl

   ! How many handles initally and how much increase once exhausted
   integer, parameter :: initial_handles = 5
   double precision, parameter :: multiplier = 2.0

   ! SAVEd data
   integer, save :: next_handle = 1
   integer, save :: total_handles = 0
   type(ma48_hdl), dimension(:), allocatable, save :: handles

contains

   integer function ma48_new_handle()

      type(ma48_hdl), dimension(:), allocatable :: temp
      integer :: i

      ! Do we need to expand the number of available handles?
      if (next_handle .gt. total_handles) then
         if(total_handles.ne.0) then
            ! Need to expand existing handle selection
            allocate(temp(total_handles))
            do i = 1, total_handles
               temp(i)%matrix => handles(i)%matrix
               temp(i)%factors => handles(i)%factors
            end do
            deallocate(handles)
            total_handles = max(int(multiplier*total_handles), total_handles)
            allocate(handles(total_handles))
            do i = 1, size(temp)
               handles(i)%matrix => temp(i)%matrix
               handles(i)%factors => temp(i)%factors
            end do
            deallocate(temp)
         else
            ! First call since module loaded
            total_handles = initial_handles
            allocate(handles(total_handles))

            ! Register clean function
            call mexAtExit(cleanup_all_handles)
         endif
      endif

      ma48_new_handle = next_handle
      allocate(handles(next_handle)%matrix)
      allocate(handles(next_handle)%factors)
      next_handle = next_handle + 1
   end function ma48_new_handle

   ! This routine is called at unload of this module from MATLAB.
   ! It shuold cleanup all SAVEd data
   subroutine cleanup_all_handles()
      integer :: i

      do i = 1, next_handle-1
         call cleanup_handle(i)
      end do
   end subroutine cleanup_all_handles

   ! Destroy the data associated with a handle.
   ! Recover all free pointers at end of handle list.
   subroutine cleanup_handle(handle)
      integer, intent(in) :: handle

      integer :: current

      if(associated(handles(handle)%factors)) then
         deallocate(handles(handle)%matrix)
         nullify(handles(handle)%matrix)
         call ma48_matlab_finalise(handles(handle)%factors)
         deallocate(handles(handle)%factors)
         nullify(handles(handle)%factors)
      endif

      do current = handle, 1, -1
         if(current.ne.next_handle-1) exit
         if(.not.associated(handles(current)%factors)) then
            ! Current "last" element is unallocated, make it next available
            ! element.
            next_handle = next_handle - 1
         endif
      end do
   end subroutine cleanup_handle

end module ma48_handles


! Gateway routine
! Strip first argument and converts any handles, then calls relevant routine
subroutine mexFunction(nlhs_in, plhs, nrhs_in, prhs)
   use hsl_matlab
   use hsl_ma48_double
   use ma48_handles
   use ma48_matlab_main
   implicit none

   integer*4 :: nlhs_in, nrhs_in
   integer(mwp_) :: plhs(*), prhs(*)
   integer(mws_) :: mwstemp

   character(len=200) :: act
   integer :: handle

   character(len=200) :: errmsg

   if(nrhs_in.lt.1) call MATLAB_error("Insufficient input arguments")

   if(.not. MATLAB_is_character(prhs(1))) &
      call MATLAB_error("First argument must be string")
   mwstemp = len(act)
   call MATLAB_get_string(prhs(1), act, mwstemp)

   select case(trim(act))
   case('factor')
      ! [handle, info] = ma48_expert('factor', A[, P])
      ! At least handle is required for output

      ! Setup handle and store in first output argument
      if(nlhs_in.lt.1) call MATLAB_error("Insufficient output arguments")
      handle = ma48_new_handle()
      plhs(1) = fortran_to_matlab(handle)

      ! Call work routine
      call ma48_matlab_analyse_factor(nlhs_in-1_int4_, plhs(2), &
         nrhs_in-1_int4_, prhs(2), handles(handle)%matrix, &
         handles(handle)%factors)

   case('solve')
      ! [x, info] = ma48_expert('solve', handle, b)

      ! Obtain and check handle
      if(nrhs_in.lt.2) call MATLAB_error("Insufficient input arguments")
      call matlab_to_fortran(prhs(2), handle, 'handle')
      if(handle.lt.1 .or. handle.gt.total_handles) &
         call MATLAB_error("Invalid handle")

      ! Call worker routine based on allocated parts of handle
      if(associated(handles(handle)%factors)) then
         call ma48_matlab_solve(nlhs_in, plhs(1), nrhs_in-2_int4_, &
            prhs(3), handles(handle)%matrix, handles(handle)%factors)
      else
         ! factors not allocated, probably a destroyed or unallocated handle
         call MATLAB_error("Invalid handle")
      endif
      
   case('destroy')
      ! ma48_expert('destroy', handle)

      ! Obtain and check handle
      if(nrhs_in.lt.2) call MATLAB_error("Insufficient input arguments")
      call matlab_to_fortran(prhs(2), handle, 'handle')
      if(handle.lt.1 .or. handle.gt.total_handles) &
         call MATLAB_error("Invalid handle")

      ! Destroy anything that is there
      call cleanup_handle(handle)

   case default
      write(errmsg, "(3a)") "Unrecognised action: '", trim(act), "'"
      call MATLAB_error(errmsg)
   end select

end subroutine mexFunction
