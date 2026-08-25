module mod_dream_input
    use constants, only: MU_ZERO
    implicit none

    private
    public :: set_dream_output, interp_dream_jre, n_dream_psin, load_dream_output_from_file
    public :: dream_psin, dream_jre_par, dream_bmin

    integer, save :: n_dream_psin = 0
    real*8, allocatable, save :: dream_psin(:)
    real*8, allocatable, save :: dream_jre_par(:)
    real*8, allocatable, save :: dream_bmin(:)

contains

    subroutine load_dream_output_from_file(filename, ierr)
        use hdf5_io_module
        character(len=*), intent(in) :: filename
        integer, intent(out) :: ierr
        integer(HID_T) :: file_id
        real*8, allocatable :: psin_file(:), jre_file(:), bmin_file(:)

        call HDF5_open(filename, file_id, ierr)
        if (ierr /= 0) return

        call HDF5_allocatable_array1D_reading(file_id, psin_file, 'jorek/psi_n')
        call HDF5_allocatable_array1D_reading(file_id, jre_file,  'jorek/jre_on_jorek_grid')
        call HDF5_allocatable_array1D_reading(file_id, bmin_file, 'jorek/bmin_on_jorek_grid')

        call HDF5_close(file_id)

        if (.not. allocated(psin_file) .or. .not. allocated(jre_file) &
            .or. .not. allocated(bmin_file)) then
            write(*,*) 'ERROR: failed to read jorek/* datasets from ', trim(filename)
            ierr = 1
            return
        end if

        call set_dream_output( &
            psin_arr = psin_file, &
            jre_arr  = jre_file, &
            bmin_arr = bmin_file, &
            n_pts    = size(psin_file) )

        ierr = 0
    end subroutine load_dream_output_from_file

    subroutine set_dream_output(psin_arr, jre_arr, bmin_arr, n_pts)
        integer, intent(in) :: n_pts
        real*8,  intent(in) :: psin_arr(n_pts)
        real*8,  intent(in) :: jre_arr(n_pts)
        real*8,  intent(in) :: bmin_arr(n_pts)

        n_dream_psin = n_pts

        if (allocated(dream_psin))    deallocate(dream_psin)
        if (allocated(dream_jre_par)) deallocate(dream_jre_par)
        if (allocated(dream_bmin))    deallocate(dream_bmin)

        allocate(dream_psin(n_pts))
        allocate(dream_jre_par(n_pts))
        allocate(dream_bmin(n_pts))

        dream_psin    = psin_arr
        dream_jre_par = jre_arr
        dream_bmin    = bmin_arr

    end subroutine set_dream_output


    function interp_dream_jre(psi_norm_jorek, B_local, R_local, F0_local) result(jre)
        real*8, intent(in) :: psi_norm_jorek
        real*8, intent(in) :: B_local
        real*8, intent(in) :: R_local   ! local major radius, matches zj0's R-scaling convention
        real*8, intent(in) :: F0_local  ! R*B_phi (toroidal field function) at the local point
        real*8 :: jre

        integer :: i
        real*8  :: frac, jre_dream, bmin_interp

        if (n_dream_psin < 1) then
            jre = 0.d0
            return
        end if

        if (psi_norm_jorek <= dream_psin(1)) then
            jre_dream   = dream_jre_par(1)
            bmin_interp = dream_bmin(1)
        else if (psi_norm_jorek >= dream_psin(n_dream_psin)) then
            jre = 0.d0   ! outside LCFS — no runaways
            return
        else
            do i = 1, n_dream_psin - 1
                if (psi_norm_jorek <= dream_psin(i+1)) then
                    frac        = (psi_norm_jorek  - dream_psin(i)) &
                                / (dream_psin(i+1) - dream_psin(i))
                    jre_dream   = dream_jre_par(i) * (1.d0 - frac) &
                                + dream_jre_par(i+1) * frac
                    bmin_interp = dream_bmin(i) * (1.d0 - frac) &
                                + dream_bmin(i+1) * frac
                    exit
                end if
            end do
        end if

        if (bmin_interp > 0.d0) then
            !jre = (B_local / bmin_interp) * jre_dream * R_local * MU_ZERO
            jre = jre_dream * MU_ZERO * abs(F0_local) / (bmin_interp * R_local)
        else
            jre = 0.d0
        end if

    end function interp_dream_jre

end module mod_dream_input
