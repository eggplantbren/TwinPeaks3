module twin_peaks_routines ! Robert Baldock (rjnbaldock@gmail.com)

  implicit none
  integer, parameter :: dp = kind(1.0d0) ! double precision

contains

! Code sections 
! 1. MC MOVES
! 2. ROUTINE FOR CREATING AN INITIAL CONFIGURATION
! It should not be neccessary to read the last two sections
! 3. ROUTINES FOR CALCULATING THE ENTHALPY OF A LENNARD JONES SYSTEM
! 4. VARIOUS NUMERICAL FUNCTIONS AND ROUTINES REQUIRED BY THE ABOVE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! MC MOVES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !   ************************************************************
  !   * Single atom MC move                                      *
  !   ************************************************************
  subroutine roam_oneatom_onestep(s_in, s_out, h0, V, delta_enth, step_1a1s, iat, cutoff, press, n_call_enth, skip_enth_calc)

    real(dp), dimension(:,:), intent(in)   :: s_in  ! fractional coords in
    real(dp), dimension(:,:), intent(out)  :: s_out ! fractional coords out
    real(dp), dimension(3,3), intent(in)   :: h0    ! normalised cell matrix
    real(dp), intent(in)                   :: V     ! volume
    real(dp), intent(out)                  :: delta_enth ! enthalpy change due to single atom step
    real(dp), intent(in)      :: step_1a1s          ! step length
                                                    ! defined as a fraction of unit_cell_volume**(1/3)
    integer, intent(out)      :: iat                ! index of random atom that was chosen moved 
                                                    ! element of (1,2,...,N)
    real(dp), intent(in)      :: cutoff             ! radial cutoff of potential, defined in units of sigma
    real(dp), intent(in)      :: press              ! pressure defined in units of eps/(sigma**3)
    real(dp), intent(inout)   :: n_call_enth        ! number of entahlpy calls (tally) - updated by +1 
    logical, optional, intent(in) :: skip_enth_calc ! OPTIONAL. if .true. do not calculate delta_enthalpy 
    ! internal arguments
    real(dp), dimension(4) :: rnd4
    real(dp), dimension(3) :: translation_roos
    integer :: i, N

    N = size(s_in)/3

    call random_number(rnd4)
    iat = int(rnd4(1)*real(N)) + 1
    ! choose random particle between 1 and N for taking the step

    rnd4 = rnd4 * 2.0_dp * step_1a1s - step_1a1s  

    s_out = s_in
    do i = 2,4
       ! do move in fractional coordinates
       s_out(i-1,iat) = s_out(i-1,iat) + rnd4(i)
    end do

    ! map back into cell if atom has stepped out of periodic boundaries
    translation_roos = nint(s_out(:,iat) - 0.5_dp, dp)
    s_out(:,iat) = s_out(:,iat) - translation_roos

    if (.not. (present(skip_enth_calc).and.skip_enth_calc)) then 
      delta_enth = enthalpy( s_out, h0, V, press, cutoff, atom=iat, n_call_enthalpy=n_call_enth )
      delta_enth = delta_enth - enthalpy( s_in, h0, V, press, cutoff, atom=iat ) 
      ! only count one enthalpy call for this move
    end if

  end subroutine roam_oneatom_onestep

  !   ************************************************************
  !   * Change the unit lattice cell by applying random shears   *
  !   * MC accept/reject ensures uniform random sampling of the  *
  !   * manifold det(h0) = 1 and cell_depth > min_height         *
  !   ************************************************************
  subroutine roam_latticeshape_shear(h0_in, h0_out, s, V, latt_step, min_height, &
                             & pressure, cutoff, change_flag, enth_out, n_call_enth, skip_enth_calc )

    real(dp), dimension(3,3), intent(in) :: h0_in   ! normalised cell matrix in
    real(dp), dimension(3,3), intent(out) :: h0_out ! normalised cell matrix out
    real(dp), dimension(:,:), intent(in) :: s       ! fractional coordinates
    real(dp), intent(in) :: V                       ! volume
    real(dp), intent(in)       :: latt_step         ! step length - element of (0.0,1.0)
    real(dp), intent(in) :: min_height              ! minimum cell height (see paper)
    real(dp), intent(in) :: pressure                ! pressure defined in units of eps/(sigma**3)
    real(dp), intent(in) :: cutoff                  ! radial cutoff of potential, defined in units of sigma
    integer , intent(inout) :: change_flag          ! flags whether the minimum cell height condition was met.
                                                    ! 0 : YES 
                                                    ! 1 : NO 
                                                    ! enth_out has ONLY been calculated if change_flag=0

    real(dp), intent(out)   :: enth_out             ! enthalpy of trial h0
    real(dp), intent(inout) :: n_call_enth          ! tally of enthalpy calls.
    logical, optional, intent(in) :: skip_enth_calc ! if .true. do not perform enthalpy calculation

    real(dp) :: volume_curt, rand_nm, rnd2(2), vol_save
    real(dp), dimension(9) :: rnd_matrix
    real(dp), dimension(3) :: delta_vec, orth1, orth2
    real(dp), dimension(3,3) :: latt_save
    integer  :: i,i3, chosen_lv, roam_count, vs(2)

    h0_out = h0_in

    call random_number(rnd_matrix)

    change_flag = 1 ! default value indicating that the move was not successful
    vol_save = lattice_Cell_Volume(h0_out)
    volume_curt = vol_save**(1.0_dp/3.0_dp)
    h0_out = h0_out / volume_curt 
    ! we are exploring the unit lattice cell h0
    latt_save = h0_out
    do roam_count = 1,3
      ! choose random lattice vector to roam
      i3=(roam_count-1)*3
      rand_nm = rnd_matrix(i3 + 1)
      chosen_lv = (rand_nm*3) + 1 ! from 1 to 3
      if (chosen_lv==1) then
        vs = (/ 2, 3 /)
      else if (chosen_lv==2) then
        vs = (/ 1, 3 /)
      else
        vs = (/ 1, 2 /)
      end if

      ! construct orthonormal basis from other two lattice vectors
      call orthonormal_basis(h0_out(:,vs(1)), &
                           & h0_out(:,vs(2)), orth1, orth2)

      rnd2(1) = rnd_matrix(i3 + 2)
      rnd2(2) = rnd_matrix(i3 + 3)
      ! Insist that these are from inside a unit circle, 
      ! because the direction of orth1, orth2 depend on the lattice vectors.
      ! Makes detailed balance good.
      if (sum(rnd2**2).gt.1.0_dp) then
        do
          call random_number(rnd2)
          rnd2 = 2.0_dp*(rnd2-0.5_dp)
          if (sum(rnd2**2).le.1.0_dp) exit
        end do
      end if
      rnd2 = 2.0_dp*(rnd2-0.5_dp)*latt_step
      ! uniform random numbers in a circle centred at the origin

      do i = 1,3
        delta_vec(i) = rnd2(1)*orth1(i) + rnd2(2)*orth2(i)
        h0_out(i,chosen_lv) = h0_out(i,chosen_lv) + &
                                        & delta_vec(i)
      end do

      ! test whether this has given us an allowed lattice cell shape:
      rand_nm = lattice_Cell_Volume(h0_out) ** (1.0_dp/3.0_dp)
      rand_nm=minval(lattice_height(h0_out))/rand_nm

      if ( rand_nm > min_height ) then
        change_flag = 0 ! the MC step results in a lattice shape change
        latt_save = h0_out
      else
        h0_out = latt_save ! revert to previous saved lattice
      end if
    end do
    if (change_flag == 0 .and. (.not. (present(skip_enth_calc).and.skip_enth_calc)) ) then 
      ! calc new enth
      enth_out = enthalpy( s, h0_out, V, pressure, cutoff, n_call_enthalpy=n_call_enth ) 
    end if

  end subroutine roam_latticeshape_shear

  !   ************************************************************
  !   * Change the unit lattice cell by applying random stretches*
  !   * MC accept/reject ensures uniform random sampling of the  *
  !   * manifold det(h0) = 1 and cell_depth > min_height         *
  !   ************************************************************
  subroutine roam_latticeshape_stretch(h0_in, h0_out, s, V, latt_step, min_height, &
                             & pressure, cutoff, change_flag, enth_out, n_call_enth, skip_enth_calc )

    real(dp), dimension(3,3), intent(in) :: h0_in   ! normalised cell matrix in
    real(dp), dimension(3,3), intent(out) :: h0_out ! normalised cell matrix out
    real(dp), dimension(:,:), intent(in) :: s       ! fractional coordinates
    real(dp), intent(in) :: V                       ! volume
    real(dp), intent(in)       :: latt_step         ! step length - element of (0.0,1.0)
    real(dp), intent(in) :: min_height              ! minimum cell height (see paper)
    real(dp), intent(in) :: pressure                ! pressure defined in units of eps/(sigma**3)
    real(dp), intent(in) :: cutoff                  ! radial cutoff of potential, defined in units of sigma
    integer , intent(inout) :: change_flag          ! flags whether the minimum cell height condition was met.
                                                    ! 0 : YES 
                                                    ! 1 : NO 
                                                    ! enth_out has ONLY been calculated if change_flag=0

    real(dp), intent(out)   :: enth_out             ! enthalpy of trial h0
    real(dp), intent(inout) :: n_call_enth          ! tally of enthalpy calls.
    logical, optional, intent(in) :: skip_enth_calc ! if .true. do not perform enthalpy calculation

    real(dp) :: rand_nm, rnd1, rnd2
    real(dp), dimension(6) :: rnd_matrix
    real(dp), dimension(3,3) :: latt_save
    integer  :: i,i2, chosen_lv, roam_count, vs(2)

    h0_out = h0_in 

    call random_number(rnd_matrix)

    change_flag = 1 ! default value indicating that the move was not successful
    latt_save = h0_out ! save
    do roam_count = 1,3
      ! choose random lattice vector to roam
      i2=(roam_count-1)*2
      rand_nm = rnd_matrix(i2 + 1)
      chosen_lv = (rand_nm*3) + 1 ! from 1 to 3
      if (chosen_lv==1) then
        vs = (/ 2, 3 /)
      else if (chosen_lv==2) then
        vs = (/ 1, 3 /)
      else
        vs = (/ 1, 2 /)
      end if

      rnd1 = rnd_matrix(i2 + 2)
      rnd1 = 2.0_dp*(rnd1-0.5_dp)*latt_step
      rnd1 = exp(rnd1)
      rnd2 = 1.0_dp/rnd1

      do i = 1,3
        h0_out(i,vs(1)) = h0_out(i,vs(1))*rnd1
        h0_out(i,vs(2)) = h0_out(i,vs(2))*rnd2
      end do

      ! test whether this has given us an allowed lattice cell shape:
      rand_nm = lattice_Cell_Volume(h0_out) ** (1.0_dp/3.0_dp)
      rand_nm=minval(lattice_height(h0_out))/rand_nm

      if ( rand_nm > min_height ) then
        change_flag = 0 ! the MC step results in a lattice shape change
        latt_save = h0_out ! save
      else
        h0_out = latt_save ! revert to previous saved lattice
      end if
    end do

    if (change_flag == 0 .and. (.not. (present(skip_enth_calc).and.skip_enth_calc)) ) then 
      enth_out = enthalpy( s, h0_out, V, pressure, cutoff, n_call_enthalpy=n_call_enth ) 
    endif

  end subroutine roam_latticeshape_stretch

  !   ************************************************************
  !   * Change the volume by applying random dilation/contraction*
  !   * MC accept/rejection guarantees V ~ V^(1/3) and V<v_max .
  !   ************************************************************
  subroutine roam_volume(V_in, V_out, s, h0, vol_step, v_max, &
                             & pressure, cutoff, change_flag, enth_out, n_call_enth, skip_enth_calc )

    real(dp), intent(in) :: V_in                ! volume in
    real(dp), intent(out) :: V_out              ! volume out
    real(dp), dimension(:,:), intent(in) :: s   ! fractional coordinates
    real(dp), dimension(3,3), intent(in) :: h0  ! normalised cell matrix
    real(dp), intent(in)       :: vol_step      ! step length
    real(dp), intent(in) :: v_max               ! maximum allowed volume. "V_0" in the paper.
    real(dp), intent(in) :: pressure            ! pressure defined in units of eps/(sigma**3)
    real(dp), intent(in) :: cutoff              ! radial cutoff of potential, defined in units of sigma
    integer , intent(inout) :: change_flag      ! flags whether the trial move was accepted, 
                                                ! according to the V**(1/3) distribution
                                                ! 0 : YES 
                                                ! 1 : NO 
                                                ! enth_out has ONLY been calculated if change_flag=0
    real(dp), intent(out)   :: enth_out         ! enthalpy of trial V
    real(dp), intent(inout) :: n_call_enth      ! tally of enthalpy calls.
    logical, optional, intent(in) :: skip_enth_calc ! if .true. do not perform enthalpy calculation

    real(dp) :: trial_volume, a_0, a_1, p
    integer  :: N
    real(dp) :: rnd2(2)
    real(dp), parameter :: v_min = 0.0d0

    call random_number(rnd2)

    V_out = V_in
    N = size(s)/3

    ! calculate the closest boundary:
    a_0 = min(v_max - V_in, V_in - v_min )
    if (a_0 .gt. vol_step) then
       a_0 = vol_step
    end if

    trial_volume = V_in
    if ( v_max - V_in .lt. V_in - v_min) then
       trial_volume = trial_volume - 0.5_dp*(vol_step - a_0)
    else if (v_max - V_in .gt. V_in - v_min) then
       trial_volume = trial_volume + 0.5_dp*(vol_step - a_0)
    end if ! This naturaly leaves trial_volume = volume if 
    ! the volume is exactly in the middle of the allowed volume range

    trial_volume = trial_volume + (rnd2(1) - 0.5_dp)*(vol_step + a_0)

    ! calculate the closest boundary to trail_volume:
    a_1 = min(v_max - trial_volume, trial_volume - v_min )
    if (a_1 .gt. vol_step) then
       a_1 = vol_step
    end if

    p =  (trial_volume/V_in)**N
    p = p * (vol_step + a_0)/(vol_step + a_1)

    change_flag = 1 ! signals that the MC resulted in no volume change - default value

    if ( rnd2(2) .lt. min( 1.0_dp, p ) ) then
       change_flag = 0 ! MC step resulted in a volume change
       V_out = trial_volume
       if (.not. (present(skip_enth_calc).and.skip_enth_calc)) &
         & enth_out = enthalpy( s, h0, V_out, pressure, cutoff, n_call_enthalpy=n_call_enth ) 
    else
      enth_out = huge(1.0d0)
      V_out = V_in

    end if

  end subroutine roam_volume

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! ROUTINE FOR CREATING AN INITIAL CONFIGURATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine initialise_config(s,h0,V,enth,min_height,Vmax,N,pressure,cutoff)

  real(dp), dimension(:,:), allocatable, intent(inout) :: s ! uniformly random fractional coordinates
                                                            ! s is an element of (0,1)^(3N)
  real(dp), dimension(3,3), intent(out) :: h0               ! uniformly random normalised cell matrix,
                                                            ! chosen from det(h0)=1, and satisfying
                                                            ! cell_depth(h0) > min_height (See paper).
  real(dp), intent(out)                 :: V                ! Volume
  real(dp), intent(out)                 :: enth             ! enthalpy of this initial configuration

  real(dp), intent(in) :: Vmax                              ! maximum allowed volume. "V_0" in the paper.
  real(dp), intent(in) :: min_height                        ! minimum cell height (see paper)
  real(dp), intent(in) :: pressure                          ! pressure defined in units of eps/(sigma**3)
  real(dp), intent(in) :: cutoff                            ! radial cutoff of potential, defined in units of sigma
  integer, intent(in) :: N                                  ! Number of atoms

  ! internal arguments
  integer :: i, spare_int
  real(dp) :: spare,spare2
  real(dp), dimension(3,3) :: h0_temp
  
  integer, parameter :: N_latt_moves_eq = 10000

  ! initialise s
  if (allocated(s)) deallocate(s)
  allocate(s(3,N))
  call random_number(s) ! univariate intial coordinate s on (0,1)^(3N)

  ! initialise h0 uniformly with min cell height > min_height
  h0 = 0.0d0
  do i = 1,3
    h0(i,i) = 1.0d0
  end do
  spare2 = 0.0d0
  do i = 1, N_latt_moves_eq
    call random_number(spare)
    if (spare > 0.5d0) then ! stretch move
      call roam_latticeshape_stretch(h0, h0_temp, s, 1.0d0, 0.35d0, min_height, &
                             & 1.0d0, 1.d0, spare_int, spare, spare2, skip_enth_calc=.true. )
      h0=h0_temp
    else ! shear move
      call roam_latticeshape_shear(h0, h0_temp, s, 1.0d0, 0.5d0, min_height, &
                             & 1.0d0, 1.d0, spare_int, spare, spare2, skip_enth_calc=.true. )
      h0=h0_temp
    end if
  end do  

  ! initialise volume with prior V^N over (0,Vmax)
  call polynomial_rand(real(N,dp),Vmax,single=V)

  enth = enthalpy( s, h0, V, pressure, cutoff )

end subroutine initialise_config

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! ROUTINES FOR CALCULATING THE ENTHALPY OF A LENNARD JONES SYSTEM !!!!!!!!!!!!!!!!!!!!!!!!!!
!! EPSILON = 1.0 , SIGMA = 1.0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! STANDARD MEAN FIELD CORRECTION FOR r > cutoff !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !***********************************************************************************
  ! function to calculate the enthalpy PE + PV + longrange (mean field) correction
  ! if atom is specified, PE only contains the contribution due to the interactions of atom "atom"
  !***********************************************************************************
  function enthalpy( s, h0, V, press, cutoff, atom, n_call_enthalpy )

    real(dp) :: enthalpy

    real(dp), dimension(:,:), intent(in) :: s ! fractional coordinates
    real(dp), dimension(3,3), intent(in) :: h0 ! normalised cell matrix
    real(dp), intent(in) :: V ! volume
    real(dp), intent(in) :: press ! units of eps/(sigma**3)
    real(dp), intent(in) :: cutoff ! units of sigma
    integer, optional, intent(in)    :: atom
    real(dp), optional, intent(inout)    :: n_call_enthalpy

    real(dp), parameter :: pi = 3.14159265358979323846264338327950288419716939937510_dp
    real(dp) :: mfc

    if ( present(atom) )  then
       enthalpy = lj_energy_correct_longrange( &
         & s,h0,V,cutoff,iat=atom)
    else ! atom not specified
       enthalpy = lj_energy_correct_longrange( &
         & s,h0,V,cutoff)
    end if

    ! add mean field correction beyond cutoff
    mfc = (2.0_dp/3.0_dp)*pi*4.0d0*dble((size(s)/3)**2)
    mfc = mfc*((((1.0d0/cutoff)**9)/3.0_dp)-((1.0d0/cutoff)**3))/V
    enthalpy = enthalpy + mfc

    enthalpy = enthalpy + ( press * V )

    if(present(n_call_enthalpy)) n_call_enthalpy = n_call_enthalpy + 1.0

  end function enthalpy

  function lj_energy_correct_longrange(s,h0,V,cutoff,iat)
    ! Calculate explicit interaction part of LJ potential energy 
    ! (not including the mean field contribution)
    ! If <iat> is peresent, only the contribution of that atom is calculated***

    real(dp) :: lj_energy_correct_longrange

    real(dp), dimension(:,:), intent(in) :: s
    real(dp), dimension(3,3), intent(in) :: h0
    real(dp), intent(in) :: V
    real(dp), intent(in) :: cutoff
    integer, optional, intent(in)    :: iat

    real(dp), allocatable, dimension(:,:) :: spatial_pos 
    real(dp), allocatable, dimension(:,:,:,:) :: lattice_vector 
    real(dp), allocatable, dimension(:) :: pair_energy 

    real(dp) :: ener_ij, dist_ij2
    real(dp), dimension(3) :: diff_ij, diff_ij_shift, spatial_diff_ij, latt_vec
    integer, dimension(3) :: n_shift
    integer :: i, j, k, ix, iy, iz
    integer :: nx_mod_low, nx_mod_high, ny_mod_low, ny_mod_high, nz_mod_low, nz_mod_high
    ! particle i (or j, respectively) is of type A, or B (also respectively).
    real(dp) :: coeff_1,coeff_2, d_2, d_6, cutoff2
    integer :: pair_energy_pos_ij, N
    real(dp), dimension(3) :: threevec

    N = size(s)/3

    if (present(iat) .and. iat > N ) then
       write(*,*) "lj_energy: ERROR, optional argument", &
            &" iat is larger than the number of atoms in the system. STOP"
       STOP
    endif

    n_shift = CEILING( cutoff / (lattice_height(h0)*(V**(1.0d0/3.0d0))) ) 
  
    lj_energy_correct_longrange = 0.0_dp    

    if(allocated(lattice_vector)) deallocate( lattice_vector)
    allocate( lattice_vector(3, -n_shift(1):n_shift(1), -n_shift(2):n_shift(2), &
         & -n_shift(3):n_shift(3)) )

    do iz = -n_shift(3) , +n_shift(3)
       do iy = -n_shift(2) , +n_shift(2)
          do ix = -n_shift(1) , +n_shift(1)
             threevec = mat3x3_vec3(h0,real((/ix, iy, iz/),dp)) 
                      ! careful - h0 is the normalised unit cell matrix
             do i = 1, 3
                lattice_vector(i, ix, iy, iz) = threevec(i) 
             end do
          end do
       end do
    end do
    lattice_vector = lattice_vector*(V**(1.0d0/3.0d0)) ! map to spatial dimensions (from h0)

    if(allocated(pair_energy)) deallocate( pair_energy)
    allocate( pair_energy(( N**2 + N )/2) )
    pair_energy = 0.0d0

    if(allocated(spatial_pos)) deallocate( spatial_pos)
    allocate( spatial_pos(3, N) )
    do i = 1, N ! calculate spatial positions
      threevec = s(:,i)
      spatial_pos(:, i) = mat3x3_vec3(h0,threevec) ! careful - h0 is the normalised unit cell matrix
    end do
    spatial_pos = spatial_pos*(V**(1.0d0/3.0d0)) ! map to spatial dimensions (from h0)

    coeff_1 = 1.0d0 ! coeff_1 = sig**6 = 1.0
    coeff_2 = 4.0d0 ! coeff_2 = 4*eps*(sig**6)
    cutoff2 = cutoff**2

    jth: do j = 1, N

       ith: do i = j, N

          if (present(iat) .and. j /= iat .and. i /= iat) then
             cycle
          endif

          diff_ij = s(:,j) - s(:,i)

          if (ceiling(diff_ij(1)) .gt. 0) then
             nx_mod_low = 0 
             nx_mod_high = -1
          else
             nx_mod_low = 1
             nx_mod_high = 0
          end if
          if (ceiling(diff_ij(2)) .gt. 0) then
             ny_mod_low = 0 
             ny_mod_high = -1
          else
             ny_mod_low = 1
             ny_mod_high = 0
          end if
          if (ceiling(diff_ij(3)) .gt. 0) then
             nz_mod_low = 0 
             nz_mod_high = -1
          else
             nz_mod_low = 1
             nz_mod_high = 0
          end if

          spatial_diff_ij =  spatial_pos(:,j) - spatial_pos(:,i)    

          do iz = -n_shift(3) + nz_mod_low, n_shift(3) + nz_mod_high
             do iy = -n_shift(2) + ny_mod_low, n_shift(2) + ny_mod_high
                do ix = -n_shift(1) + nx_mod_low, n_shift(1) + nx_mod_high

                   do k = 1,3
                      latt_vec(k) = lattice_vector(k, ix, iy, iz)
                   end do

                   diff_ij_shift = spatial_diff_ij + latt_vec

                   dist_ij2 = 0.0_dp
                   do k =1,3
                      dist_ij2 = dist_ij2  + diff_ij_shift(k)**2
                   end do

                   if (dist_ij2 > cutoff2) cycle
                   if (dist_ij2 < 1.0e-50_dp) cycle

                   d_2 = 1.0_dp/dist_ij2
                   d_6 = d_2**3

                   ener_ij = coeff_2 * d_6 * ( coeff_1*d_6 - 1.0_dp )

                   pair_energy_pos_ij = i + (2*N-j)*(j-1)/2          

                   pair_energy( pair_energy_pos_ij ) = &
                        & pair_energy( pair_energy_pos_ij ) + ener_ij

                end do
             end do
          end do

       end do ith
    end do jth

    if(present(iat) ) then
       pair_energy_pos_ij = iat + (2*N-iat)*(iat-1)/2
       pair_energy(pair_energy_pos_ij) = pair_energy(pair_energy_pos_ij)*0.5_dp

       pair_energy_pos_ij =  iat + (2*N-iat)*(iat-1)/2
       do i = iat, N
          lj_energy_correct_longrange = lj_energy_correct_longrange + pair_energy(pair_energy_pos_ij)
          pair_energy_pos_ij = pair_energy_pos_ij + 1
       end do
       do j = 1, iat - 1
          pair_energy_pos_ij = iat + (2*N-j)*(j-1)/2   
          lj_energy_correct_longrange = lj_energy_correct_longrange + pair_energy(pair_energy_pos_ij)
       end do

    else
       do i = 1, N
          pair_energy_pos_ij = i + (2*N-i)*(i-1)/2
          pair_energy(pair_energy_pos_ij) = pair_energy(pair_energy_pos_ij)*0.5_dp
       end do

       lj_energy_correct_longrange = 0.0_dp
       do i = 1, int(0.5*N*(N+1))
          lj_energy_correct_longrange = lj_energy_correct_longrange + pair_energy(i)
       end do
    end if

  end function lj_energy_correct_longrange

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! VARIOUS NUMERICAL FUNCTIONS AND ROUTINES REQUIRED BY THE ABOVE !!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  pure function mat3x3_vec3(A,x)

    real(dp) :: mat3x3_vec3(3)
    real(dp), intent(in) :: A(3,3), x(3)

    mat3x3_vec3(1) = A(1,1)*x(1)
    mat3x3_vec3(2) = A(2,1)*x(1)
    mat3x3_vec3(3) = A(3,1)*x(1)

    mat3x3_vec3(1) = mat3x3_vec3(1) + A(1,2)*x(2)
    mat3x3_vec3(2) = mat3x3_vec3(2) + A(2,2)*x(2)
    mat3x3_vec3(3) = mat3x3_vec3(3) + A(3,2)*x(2)

    mat3x3_vec3(1) = mat3x3_vec3(1) + A(1,3)*x(3)
    mat3x3_vec3(2) = mat3x3_vec3(2) + A(2,3)*x(3)
    mat3x3_vec3(3) = mat3x3_vec3(3) + A(3,3)*x(3)

  end function mat3x3_vec3

  function lattice_height(lattice)

    real(dp), dimension(3,3), intent(in) :: lattice
    real(dp), dimension(3)  :: lattice_height

    real(dp), dimension(3)  :: vec
    real(dp) :: V

    V = lattice_Cell_Volume(lattice)

    vec = cross_product( lattice(:,1), lattice(:,2) )
    lattice_height(3) = V / vec3_length( vec )

    vec = cross_product( lattice(:,3), lattice(:,1) )
    lattice_height(2) = V / vec3_length( vec )

    vec = cross_product( lattice(:,2), lattice(:,3) )
    lattice_height(1) = V / vec3_length( vec )

  end function lattice_height

  pure function lattice_cell_volume(lattice)

    real(dp), intent(in)   :: lattice(3,3)
    real(dp)               :: lattice_Cell_Volume

    real(dp), dimension(3)  :: a, b, c, bc

    a = lattice(:,1)
    b = lattice(:,2)
    c = lattice(:,3)

    bc = cross_product(b,c)

    lattice_Cell_Volume = abs(inner_prod_3(a,bc))

  end function lattice_cell_volume

  pure function cross_product(vec1,vec2)

    real(dp), dimension(3), intent(in) :: vec1, vec2
    real(dp), dimension(3)             :: cross_product

    cross_product(1) = vec1(2) * vec2(3) - vec1(3) * vec2(2)
    cross_product(2) = vec1(3) * vec2(1) - vec1(1) * vec2(3)
    cross_product(3) = vec1(1) * vec2(2) - vec1(2) * vec2(1)

  end function cross_product

  pure function inner_prod_3(vec1,vec2)

    real(dp) :: inner_prod_3
    real(dp), dimension(3), intent(in) :: vec1,vec2

    inner_prod_3 = vec1(1)*vec2(1) + vec1(2)*vec2(2) + vec1(3)*vec2(3)

  end function inner_prod_3

  function vec3_length(vec)

    real(dp)              :: vec3_length
    real(dp), dimension(3), intent(in) :: vec

    vec3_length = sqrt( inner_prod_3( vec, vec ) )

  endfunction vec3_length

  !   ************************************************************
  !   * Calculate orthonormal 3-vectors v_out_1, v_out_2         *
  !   * co-planar to 3-vectors v_in_1, v_in_2                    *
  !   ************************************************************
  subroutine orthonormal_basis(v_in_1, v_in_2, v_out_1, v_out_2)

    implicit none

    real(dp), dimension(3), intent(in) :: v_in_1, v_in_2
    real(dp), dimension(3), intent(out) :: v_out_1, v_out_2

    real(dp) :: spare1
    integer :: i

    ! v_out_1 = v_in_1/|v_in_1|
    ! v_out_2 = v_in_2 - (v_in_2 . v_out_1) v_out_1
    ! v_out_2 = v_out_2/|v_out_2|

    spare1 = sqrt(inner_prod_3(v_in_1,v_in_1))

    do i = 1,3
      v_out_1(i) = v_in_1(i)/spare1
    end do

    spare1 = inner_prod_3(v_in_2,v_out_1)

    do i = 1,3
      v_out_2(i) = v_in_2(i) - (spare1*v_out_1(i))
    end do

    spare1 = sqrt(inner_prod_3(v_out_2,v_out_2))

    do i = 1,3
      v_out_2(i) = v_out_2(i) /spare1
    end do

  end subroutine orthonormal_basis

  subroutine polynomial_rand(expo,maximum,single,array)
    ! generates random numbers distributed as
    ! p(x) = x**expo : 0 <= x <= maximum

    double precision, intent(in) :: maximum, expo
    double precision, optional, intent(in out) :: single
    double precision, optional, intent(in out) :: array(:)

    logical :: present_flag = .false.
    double precision :: exp_to_use
   
    if(abs(expo+1.0d0)>1.0d-12) then
      exp_to_use = 1.0d0/(expo+1.0d0)
    else
      write(*,*) "Error, subroutine polynomial_rand: ", &
                & "expo=",expo," Values at or very close to -1.0d0 are not allowed.", &
                & "Skipping rest of routine."
      return
    end if

    if (present(single)) then 
      call random_number(single)
      single = single**exp_to_use
      single = single*maximum
      present_flag = .true.
    end if

    if (present(array)) then 
      call random_number(array)
      array = array**exp_to_use
      array = array*maximum
      present_flag = .true.
    end if

    if (.not. present_flag) then
      write(*,*) "Error, rands.F95, subroutine polynomial_rand: ", &
                & "Routine called without specifying either single or array. ", &
                & "No random numbers generated."
    end if

  end subroutine polynomial_rand

end module twin_peaks_routines
