!> Module to determine flux surfaces, inspired by the 'fluxsurface' command in exec_commands.f90.
!> Provides a subroutine to determine (R,Z) values of the flux surfaces to be given as input to DREAM

module mod_fluxsurf_compute
    use mod_fields_linear  ! For fields_linear type
    use mod_expression     ! For expression evaluation
    use mod_position
    use mod_four_filter
    use mod_diag_output
    use equil_info
    use settings
    use convert_character
    implicit none


    private
    public :: fluxsurface

contains

  subroutine fluxsurface(ES, node_list, element_list, psi_list, R_mat, Z_mat, ierr)
  
    ! --- Routine parameters
    type(t_equil_state), intent(in)      :: ES
    type(type_node_list), intent(in)     :: node_list
    type(type_element_list), intent(in)  :: element_list
    real*8, intent(in)                   :: psi_list(:)   !< Input array of normalized psi
    integer, intent(out)                 :: ierr          !< Error flag
    real*8, allocatable, intent(out)     :: R_mat(:,:)    !< Output Matrix (Points, Psi_Index)
    real*8, allocatable, intent(out)     :: Z_mat(:,:)    !< Output Matrix (Points, Psi_Index)
    
    ! --- Local variables
    integer :: i, j, k, i_elm, ip, nplot, n_psi, i_pt, max_pieces
    type (type_surface_list) :: surface_list
    real*8 :: u, si, dsi, ti, dti, target_psi
    real*8 :: R, R_s, R_t, R_st, R_ss, R_tt
    real*8 :: Z, Z_s, Z_t, Z_st, Z_ss, Z_tt
    real*8 :: ss1, dss1, ss2, dss2, tt1, dtt1, tt2, dtt2

    ierr = 0
    n_psi = size(psi_list)
    nplot = 5  ! Number of plot points per segment

    ! Validate Input
    if (any(psi_list < 0.0d0) .or. any(psi_list > 1.0d0)) then
        write(*,*) 'Error: All target psi values must be between 0 and 1.'
        ierr = -1
        return
    endif

    ! Initialize surface_list and convert normalized psi to actual psi
    surface_list%n_psi = n_psi
    allocate(surface_list%psi_values(n_psi))
    
    do k = 1, n_psi
        surface_list%psi_values(k) = ES%psi_axis + psi_list(k) * (ES%psi_bnd - ES%psi_axis)
    end do

    ! Find all flux surfaces
    call find_flux_surfaces(0, xpoint, xcase, node_list, element_list, surface_list)

    ! Determine dimensions for the matrices
    ! Different psi values might have different numbers of pieces. We find the max.
    max_pieces = 0
    do k = 1, n_psi
        max_pieces = max(max_pieces, surface_list%flux_surfaces(k)%n_pieces)
    end do

    ! Allocate matrices: Rows = Total Points, Cols = Which Psi
    allocate(R_mat(max_pieces * nplot, n_psi), stat=ierr)
    allocate(Z_mat(max_pieces * nplot, n_psi), stat=ierr)
    
    ! Initialize with 0.0 or NaN in case some surfaces have fewer pieces than max_pieces
    R_mat = 0.0d0
    Z_mat = 0.0d0

    ! Fill the matrices
    do k = 1, n_psi
        i_pt = 0 ! Reset point counter for each column
        
        do j = 1, surface_list%flux_surfaces(k)%n_pieces
            
            ! Get segment information
            i_elm = surface_list%flux_surfaces(k)%elm(j)
            ss1  = surface_list%flux_surfaces(k)%s(1,j)
            dss1 = surface_list%flux_surfaces(k)%s(2,j)
            ss2  = surface_list%flux_surfaces(k)%s(3,j)
            dss2 = surface_list%flux_surfaces(k)%s(4,j)
            
            tt1  = surface_list%flux_surfaces(k)%t(1,j)
            dtt1 = surface_list%flux_surfaces(k)%t(2,j)
            tt2  = surface_list%flux_surfaces(k)%t(3,j)
            dtt2 = surface_list%flux_surfaces(k)%t(4,j)

            ! Interpolate points along the segment
            do ip = 1, nplot
                u = -1.0d0 + 2.0d0 * real(ip-1, 8) / real(nplot-1, 8)
                
                call CUB1D(ss1, dss1, ss2, dss2, u, si, dsi)
                call CUB1D(tt1, dtt1, tt2, dtt2, u, ti, dti)
                
                call interp_RZ(node_list, element_list, i_elm, si, ti, &
                                R, R_s, R_t, R_st, R_ss, R_tt, &
                                Z, Z_s, Z_t, Z_st, Z_ss, Z_tt)
                
                i_pt = i_pt + 1
                R_mat(i_pt, k) = R
                Z_mat(i_pt, k) = Z
            end do
        end do
    end do

    ! 6. Clean up internal surface_list
    if (allocated(surface_list%psi_values)) deallocate(surface_list%psi_values)

  end subroutine fluxsurface


end module mod_fluxsurf_compute
