 do atom1=1,TotAtom-1
   do atom2=atom1+1,TotAtom
     dr=r(atom1,:)-r(atom2,:)
     dr=dr-Box*anint(dr/Box)
     r2=dot_product(dr,dr)
     if(r2<=R2cut) then          ! r2cut  is square of rcut  
       r2=1/r2 
       fac2=r2*Sig*Sig 
       fac6=fac2*fac2*fac2 
       df=48.d0*Eps*r2*fac6*(fac6-0.5d0)
       fc(1)=df*dr(1)
       fc(2)=df*dr(2)
       fc(3)=df*dr(3)
       Force(atom1,:)=Force(atom1,:)+fc(:) 
       Force(atom2,:)=Force(atom2,:)-fc(:) 
       PE=PE+4.d0*Eps*fac6*(fac6-1)-Ecut        !shifted to zero at cutoff 
     endif 
   enddo
 enddo