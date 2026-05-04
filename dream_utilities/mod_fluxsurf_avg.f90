!> Module for flux surface averaging, inspired by the 'average' command in exec_commands.f90.
!> Provides a subroutine to average a given expression over a flux surface at a specified psi.

module mod_fluxsurf_avg
    use mod_fields_linear  ! For fields_linear type
    use mod_expression     ! For expression evaluation
    use mod_position
    use mod_four_filter
    use mod_diag_output
    use equil_info
    use settings
    use constants
    implicit none


    private
    public :: avg_fluxsurf_list

contains

    !> Subroutine to compute averages of an expression over multiple flux surfaces.
    !! @param fields The fields object containing the equilibrium data.
    !! @param expr The expression string to evaluate (e.g., 'temperature').
    !! @param psi_list Array of poloidal flux values for the flux surfaces.
    !! @param avg_vals Output array of average values (length N_psi).
    !! @param ierr Error flag (0 on success).
    subroutine avg_fluxsurf_list(ES, node_list, element_list, psi_list, avg_vals, flux_av, ierr)

    ! --- Routine parameters
    type(t_equil_state), intent(in) :: ES
    type(type_node_list), intent(in) :: node_list
    type(type_element_list), intent(in) :: element_list
    real*8, intent(in) :: psi_list(:)
    integer, intent(out) :: ierr
    real*8, allocatable, intent(out) :: avg_vals(:,:)
    logical, optional, intent(in) :: flux_av     !< Perform proper flux average

    ! --- Local variables
    integer :: units, npts, nsmall, i_exp, nmaxstep, NTht
    real*8  :: deltaphi, PsiNmin, PsiNmax
    type(t_pol_pos_list) :: pol_pos_list
    type(t_tor_pos_list) :: tor_pos_list
    type(t_expr_list) :: expr_list
    character(len=16)  :: name_T, name_ne, name_Epar
    character(len=128) :: desc_T, desc_ne, desc_Epar
    real*8 :: conv_ne, conv_T, conv_E
    real*8, allocatable :: result(:,:,:,:), res1d(:,:)

    if (allocated(result)) deallocate(result)
    if (allocated(res1d)) deallocate(res1d)


    ierr = 0
    
    npts = size(psi_list)
    nsmall = 3
    nmaxstep = 2500
    deltaphi = 0.3
    PsiNmin = psi_list(1)
    PsiNmax = psi_list(size(psi_list))
    nTht = max(150,6*n_plane)
    
    name_T = 'T'
    name_ne = 'ne'
    name_Epar = 'E_||'
    desc_T = 'temperature'
    desc_ne = 'ne'
    desc_Epar = 'parallel E field'

    call add(expr_list, name_T, desc_T)
    call add(expr_list, name_ne, desc_ne)
    call add(expr_list, name_Epar, desc_Epar)

    if (present(flux_av)) then
      call add(expr_list, 'unity       ', 'Just unity, used to get R^2 average                   ')
    endif

    pol_pos_list = pol_pos(node_list, element_list, ES, nPsiN=npts, nTht=nTht, nsmallsteps=nsmall, nmaxsteps=nmaxstep, deltaphi=deltaphi, PsiNmax=PsiNmax, PsiNmin=PsiNmin )
    tor_pos_list = tor_pos(nphi=max(n_plane,2))

    call eval_expr(ES, units, expr_list, pol_pos_list, tor_pos_list, result, ierr, flux_av)
    call apply_four_filter(result, simple_filter(m=0,n=0), expr_list%n_coord, ierr)
    call reduce_result_to_1d(ierr, result, res1d, i1=1, i2=1)


    if (present(flux_av)) then 
      if (flux_av) then
        do i_exp=1, expr_list%n_expr
          res1d(:,i_exp) = res1d(:,i_exp) / res1d(:,expr_list%n_expr)  ! Need to normalize for flux average
        enddo
      endif      
    endif

    conv_T = 1.d0 / ( MU_zero * central_density * 1.d20 * EL_CHG )
    conv_ne = 1.d20
    conv_E = 1.d0 / sqrt(MU_zero*central_density *1.d20 * central_mass * ATOMIC_MASS_UNIT)
    res1d(:,1) = res1d(:,1) * conv_T
    res1d(:,2) = res1d(:,2) * conv_ne
    res1d(:,3) = res1d(:,3) * conv_E



    allocate( avg_vals(size(res1d,1), size(res1d,2)) )
    avg_vals = res1d
    
    end subroutine avg_fluxsurf_list
    
end module mod_fluxsurf_avg