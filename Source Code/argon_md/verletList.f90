subroutine update_verlet_list(TotAtom, r, Rcut, vL,Box)
    implicit none
    integer, parameter :: dp=kind(0.d0)
    integer, intent(in) :: TotAtom
    real(kind=dp), intent(in) :: r(TotAtom, 3)
    real(kind=dp), intent(in) :: Rcut, Box
    integer, intent(inout) :: vL(TotAtom,TotAtom)

    integer :: i, j, count
    real(kind=dp) :: dr(3), r2, rv

    rv = Rcut + 2.0_dp  ! Adjust the buffer distance for the Verlet list

    ! Loop through particles and update Verlet list
    do i = 1, TotAtom
        ! Reset Verlet list for particle i
        count = 0

        ! Check distances from other particles and update Verlet list
        do j = i+1, TotAtom
            
            dr = r(i, :) - r(j, :)
            dr = dr - Box * anint(dr / Box)
            r2 = dot_product(dr, dr)

            if (r2 <= rv**2) then
                count = count + 1
                vL(i, count+1) = j
                
            endif
        end do

        ! Store the count in the first column of the vL array
        vL(i, 1) = count
    end do
end subroutine update_verlet_list
