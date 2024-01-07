subroutine force_calc(TotAtom,Box,Rcut,r,Sig,Eps,Force,PE,vL)
 use general, only: dp,atom1,atom2
 implicit none 
 integer,intent(in) :: TotAtom
 real(kind=dp),intent(in) :: r(TotAtom,3)
 real(kind=dp) :: fac2,fac6,r2,df,fc(3),dr(3),Ecut,R2cut
 real(kind=dp),intent(out) :: Force(TotAtom,3),PE
 integer, intent(inout) :: vL(TotAtom,TotAtom)

 integer::i,j,neighbor_index
 real(kind=dp) :: Box,Rcut,Eps,Sig

 R2cut=Rcut*Rcut 

 ! calculating Ecut at Rcut 
 fac2=Sig*Sig/R2cut
 fac6=fac2*fac2*fac2
 Ecut=4.d0*Eps*fac6*(fac6-1)

 ! calculting force and energy 


 PE=0.d0
 Force=0.d0

do i = 1, TotAtom
        do j = 1, vL(i, 1)  ! Assuming vL(i, 1) stores the number of neighbors for atom i
            neighbor_index = vL(i, j + 1)  ! Adding 1 to skip the first element which stores the number of neighbors
            dr = r(i, :) - r(neighbor_index, :)
            ! Apply periodic boundary conditions
            dr = dr - Box * anint(dr / Box)
            r2 = dot_product(dr, dr)
            ! Check if the neighbor is within the cutoff radius
            if (r2 <= R2cut) then
                r2 = 1.0_dp / r2
                fac2 = r2 * Sig * Sig
                fac6 = fac2 * fac2 * fac2
                df = 48.0_dp * Eps * r2 * fac6 * (fac6 - 0.5_dp)
                fc = df * dr
                ! Accumulate forces
                Force(i, :) = Force(i, :) + fc
                Force(neighbor_index, :) = Force(neighbor_index, :) - fc
                ! Accumulate potential energy
                PE = PE + 4.0_dp * Eps * fac6 * (fac6 - 1.0_dp) - Ecut ! Shifted to zero at cutoff
            end if
        end do
    end do

!print

! writing forces 
! open(unit=200,file='md.force',action='write') 
!   write(200,*)  
! do atom1=1,TotAtom
!   write(200,"(3F15.8)") Force(atom1,:)
! enddo

end subroutine force_calc 
