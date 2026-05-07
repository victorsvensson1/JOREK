module mod_dream_input
    implicit none

    private
    public :: set_dream_output, interp_dream_jre, n_dream_psin
    public :: dream_psin, dream_jre_par, dream_bmin

    integer, save :: n_dream_psin = 0
    real*8, allocatable, save :: dream_psin(:)
    real*8, allocatable, save :: dream_jre_par(:)
    real*8, allocatable, save :: dream_bmin(:)

contains

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


    function interp_dream_jre(psi_norm_jorek, B_local) result(jre)
        real*8, intent(in) :: psi_norm_jorek
        real*8, intent(in) :: B_local
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
            jre = (B_local / bmin_interp) * jre_dream
        else
            jre = 0.d0
        end if

    end function interp_dream_jre

end module mod_dream_input