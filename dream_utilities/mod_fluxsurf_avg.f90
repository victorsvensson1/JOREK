!> Module for flux surface averaging, inspired by the 'average' command in exec_commands.f90.
!> Provides a subroutine to average a given expression over a flux surface at a specified psi.

module mod_fluxsurf_avg
    use mod_fields_linear  ! For fields_linear type
    use mod_expression     ! For expression evaluation
    use mod_position
    use mod_four_filter
    use mod_diag_output
    use equil_info
    implicit none

    real*8, allocatable, private, save :: result(:,:,:,:), res1d(:,:)

    private
    public :: avg_fluxsurf_list

contains

    !> Subroutine to compute averages of an expression over multiple flux surfaces.
    !! @param fields The fields object containing the equilibrium data.
    !! @param expr The expression string to evaluate (e.g., 'temperature').
    !! @param psi_list Array of poloidal flux values for the flux surfaces.
    !! @param avg_vals Output array of average values (length N_psi).
    !! @param ierr Error flag (0 on success).
    subroutine avg_fluxsurf_list(node_list, element_list, expr_list, psi_list, avg_vals, flux_av, ierr)

    ! --- Routine parameters
    type(type_node_list), intent(in) :: node_list
    type(type_element_list), intent(in) :: element_list
    type(t_expr_list), intent(in) :: expr_list
    real*8, intent(in) :: psi_list(:)
    integer, intent(out) :: ierr
    real*8, allocatable, intent(out) :: avg_vals(:)
    logical, optional, intent(in) :: flux_av     !< Perform proper flux average

    ! --- Local variables
    integer :: units, npts, nsmall, i_exp, nmaxstep, NTht, n_plane
    real*8  :: deltaphi, PsiNmin, PsiNmax
    character(len=1024) :: filename, comment
    type(t_pol_pos_list), save :: pol_pos_list
    type(t_tor_pos_list), save :: tor_pos_list
    
    npts = 100
    nsmall = 25
    nmaxstep = 100
    deltaphi = 1.d-3
    PsiNmin = 0.d0
    PsiNmax = 1.d0
    nTht = 150
    

    pol_pos_list = pol_pos(node_list, element_list, ES, nPsiN=npts, nTht=nTht, nsmallsteps=nsmall, nmaxsteps=nmaxstep, deltaphi=deltaphi, PsiNmax=PsiNmax, PsiNmin=PsiNmin )
    tor_pos_list = tor_pos(nphi=max(n_plane,2))

    call eval_expr(ES, units, expr_list, pol_pos_list, tor_pos_list, result, ierr, flux_av)
    call apply_four_filter(result, simple_filter(m=0,n=0), expr_list%n_coord, ierr)
    call reduce_result_to_1d(ierr, result, res1d, i1=1, i2=1)

    do i_exp=1, expr_list%n_expr
        res1d(:,i_exp) = res1d(:,i_exp) / res1d(:,expr_list%n_expr)  ! Need to normalize for flux average
    enddo

    allocate( avg_vals(size(res1d,1)) )
    avg_vals = res1d(:,1)
    
    end subroutine avg_fluxsurf_list
    
end module mod_fluxsurf_avg