module mod_dream_input
    implicit none

    private
    public :: read_dream_output, interp_dream_jre, n_dream_psin
    public :: dream_psin, dream_jre_par

    integer, save :: n_dream_psin = 0
    real*8, allocatable, save :: dream_psin(:)
    real*8, allocatable, save :: dream_jre_par(:)

contains

    subroutine read_dream_output(fname, ierr)
        character(len=*), intent(in)  :: fname
        integer,          intent(out) :: ierr

        ! DUMMY: just populate with a flat zero profile on 3 points
        ! so interp_dream_jre has something to work with
        ierr = 0
        n_dream_psin = 3

        if (allocated(dream_psin))    deallocate(dream_psin)
        if (allocated(dream_jre_par)) deallocate(dream_jre_par)
        allocate(dream_psin(n_dream_psin))
        allocate(dream_jre_par(n_dream_psin))

        dream_psin    = [0.0d0, 0.5d0, 1.0d0]
        dream_jre_par = [0.0d0, 0.0d0, 0.0d0]

        write(*,*) 'DUMMY read_dream_output: returning zero j_re profile (fname ignored: '//trim(fname)//')'

    end subroutine read_dream_output


    function interp_dream_jre(psi_norm_jorek) result(jre)
        real*8, intent(in) :: psi_norm_jorek
        real*8 :: jre
        integer :: i
        real*8  :: frac

        if (n_dream_psin < 1) then
            jre = 0.d0
            return
        end if

        if (psi_norm_jorek <= dream_psin(1)) then
            jre = dream_jre_par(1)
            return
        end if
        if (psi_norm_jorek >= dream_psin(n_dream_psin)) then
            jre = dream_jre_par(n_dream_psin)
            return
        end if

        do i = 1, n_dream_psin - 1
            if (psi_norm_jorek <= dream_psin(i+1)) then
                frac = (psi_norm_jorek  - dream_psin(i)) &
                     / (dream_psin(i+1) - dream_psin(i))
                jre  = dream_jre_par(i) * (1.d0 - frac) &
                     + dream_jre_par(i+1) * frac
                return
            end if
        end do

        jre = dream_jre_par(n_dream_psin)

    end function interp_dream_jre

end module mod_dream_input