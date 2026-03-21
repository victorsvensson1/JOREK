



module mod_computeB
    use mod_fields_linear  ! For fields_linear type
    use mod_expression     ! For expression evaluation
    use mod_position
    use mod_four_filter
    use mod_diag_output
    use equil_info
    use settings
    implicit none

    real*8, allocatable, private, save :: result(:,:,:,:), res0d(:)

    private
    public :: comp_B_field

contains

    !> Subroutine to evaluate the magnetic field components Bphi, Br, Bz

    subroutine comp_B_field(ES, node_list, element_list, R_mat, Z_mat, Br_mat, Bz_mat, Bphi_mat, ierr)
    
    ! --- Routine parameters
    type(t_equil_state), intent(in) :: ES
    type(type_node_list), intent(in) :: node_list
    type(type_element_list), intent(in) :: element_list
    real*8, intent(in) :: R_mat(:,:) 
    real*8, intent(in) :: Z_mat(:,:) 
    integer, intent(out) :: ierr
    ! --- New Output Matrices
    real*8, allocatable, intent(out) :: Br_mat(:,:), Bz_mat(:,:), Bphi_mat(:,:)
    
    ! --- Local variables
    integer :: units, i, j, n_rows, n_cols
    type(t_expr_list) :: expr_list
    character(len=16)  :: name_Br, name_Bz, name_Bphi
    character(len=128) :: desc_Br, desc_Bz, desc_Bphi
    
    ierr = 0
    n_rows = size(R_mat, 1)
    n_cols = size(R_mat, 2)

    ! Allocate output matrices to match input size
    allocate(Br_mat(n_rows, n_cols), Bz_mat(n_rows, n_cols), Bphi_mat(n_rows, n_cols))

    ! Standard Setup (Units and Expressions)
    units = get_int_setting('units', ierr)
    name_Br = 'BR'
    name_Bz = 'BZ'
    name_Bphi = 'Btor'
    desc_Br = 'radial B field'
    desc_Bz = 'Z B field'
    desc_Bphi = 'toroidal B field'

    call add(expr_list, name_Br, desc_Br)
    call add(expr_list, name_Bz, desc_Bz)
    call add(expr_list, name_Bphi, desc_Bphi)

    ! Double Loop over the Matrices
    do j = 1, n_cols      ! Loop over flux surfaces (Psi)
       do i = 1, n_rows   ! Loop over points on the surface
          
          ! Skip points that are zero/uninitialized (if your matrix is padded)
          if (R_mat(i,j) == 0.0d0) cycle 
          call eval_expr(ES, units, expr_list, &
                         pol_pos(node_list, element_list, ES, R=R_mat(i,j), Z=Z_mat(i,j)), &
                         tor_pos(nphi=4*n_plane), result, ierr)
          
          call apply_four_filter(result, simple_filter(n=0), expr_list%n_coord, ierr)
          call reduce_result_to_0d(ierr, result, res0d, 1, 1, 1)
          
          ! Store in our output matrices
          Br_mat(i,j)   = res0d(1)
          Bz_mat(i,j)   = res0d(2)
          Bphi_mat(i,j) = res0d(3)
          
       end do
    end do

  end subroutine comp_B_field
    
end module mod_computeB
