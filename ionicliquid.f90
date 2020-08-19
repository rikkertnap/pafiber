! --------------------------------------------------------------|
! ionicliquid.f90:                                                    |
! constructs the vector function  needed by the                 |
! routine solver, which solves the SCMFT eqs for weak poly-     |
! electrolytes onto a tethered planar surface                   |
! --------------------------------------------------------------|


module ionicliquid

    use precision_definition

    implicit none

    contains

   
subroutine density_Im_simple()
    
    use field, only: xIm,rhoIm,xpa,xsol,psi,epsfcn
    use volume, only : nr
    use parameters, only : lb, vIm, zIm, chiIm,vsol
    use parameters, only: bornrad, expmu
    use dielectric_const, only : born

    integer :: i
    real(dp) :: lbr

    do i=1,nr  
        lbr = lb/epsfcn(i)     ! local Bjerrum length
        xIm(i)   = expmu%Im  *(xsol(i)**vIm)*exp(-born(lbr,bornrad%Im,zIm)-psi(i)*zIm)*exp(xpa(i)*chiIm)
        rhoIm(i) = xIm(i)/(vIm/vsol)
    enddo    

end subroutine

! subroutine density_Im_conf(rhopolA,rhopolACa,rhoIm_in)

!     use globals
!     use parameters
!     use volume
!     use chains
!     use field
!     use vectornorm
!     use VdW
!     use dielectric_const, only : dielectfcn, born

!     implicit none

!     !     .. scalar arguments

!     integer(8), intent(in) :: nn

!     !     .. array arguments

!     real(dp), intent(in) :: rhopolA(nsize),rhopolACa(nsize),rhoIm_in(nsize)
    
!     !     .. local variables
    
!     real(dp) :: nIm(ln_node,nsize,nsegtypes)  ! auxilairy variable for computing P(\alpha) 
!     real(dp) :: exppi(nsize,nsegtypes)        ! auxilairy variable for computing P(\alpha)
!     real(dp) :: phi(nsize)                     ! phi=1-xsol-sum_ix_i
!     real(dp) :: pro
!     real(dp) :: lbr,expborn,avgvol,Etotself,expsqrgrad, Eself
!     real(dp) :: expsqrgradpsi(nsize),expdeltaGAA(nsize,4),expEtotself(nsize)
!     integer  :: n, i,j,k,l,c,s,jmin,ln,t,k2,k3,k4       !  dummy indices
!     real(dp) :: norm
!     real(dp) :: integra_q
!     real(dp) :: xA(3),sumxA,sumxB, sgxA,qAD, constA ! disociation variables 
!     real(dp) :: xCO, xCN, bON, qON                 ! pairing interaction variables 
!     integer  :: count_sc
   
    
!     !     .. executable statements 
    
!     n=nsize

!     do t=1,nsegtypes
!         do i=1,n
!             rhoIm(i,t)=0.0_dp 
!         enddo    
!     enddo    
    
!     do l=1,n
!         q(l)=0.0_dp
!     enddo    

 
!     do i=1,n                  ! init volume fractions
!         xpol(i)    = 0.0_dp                                   ! volume fraction polymer
!         rhoqpol(i) = 0.0_dp                                   ! charge density AA monomoer

!         lbr = lb/epsfcn(i)     ! local Bjerrum length

!         Etotself = &           ! total self energy    
!             born(lbr,bornrad%AA  ,zpolAA(1))*rhopolA(i)   + & ! rhpolA(i)   = fdisA(i,1)*rhopolin(i,tA)
!             born(lbr,bornrad%AACa,zpolAA(4))*rhopolACa(i) + & ! rhopolACa(i)= fdisA(i,4)*rhopolin(i,tA)
!             born(lbr,bornrad%Na,zNa)*xNa(i)/(vNa*vsol)     + & 
!             born(lbr,bornrad%Cl,zCl)*xCl(i)/(vCl*vsol)     + &
!             born(lbr,bornrad%Rb,zRb)*xRb(i)/(vRb*vsol)     + & 
!             born(lbr,bornrad%Ca,zCa)*xCa(i)/(vCa*vsol)     + &
!             born(lbr,bornrad%Hplus,1 )*xHplus(i)/vsol      + &
!             born(lbr,bornrad%OHmin,-1)*xOHmin(i)/vsol
!             born(lbr,bornrad%Im,zIm)*rhoIm_in(i,tA)     + &

!         expEtotself(i) = Etotself*(Depsfcn(i)/epsfcn(i))  
!         !print*,"Etotself=",Etotself,"rhopolA=",rhopolA(i),"rhopolACa=",rhopolACa(i)
!         !print*,"lbr=",lbr,born(lbr,bornrad%pol,zpolAA(1)),born(lbr,bornrad%polCa,zpolAA(1))
!         !print*,"bornrad=",bornrad
!         ! pdf contribution varying dielectric constant
!         if(i/=1) then
!             expsqrgradpsi(i) = constqE* Depsfcn(i)*((psi(i+1)-psi(i-1))**2)             
!         else
!             expsqrgradpsi(i) = constqE* Depsfcn(1)*((psi(2)-psi(1))**2)   
!         endif

!         Eself=expsqrgradpsi(i)+expEtotself(i)

!     enddo

  
!     do t=1,nsegtypes
  
!         if(ismonomer_chargeable(t)) then    
!             do i=1,n  
!                 lbr=lb/epsfcn(i)
!                 expborn    = -born(lbr,bornrad%pol,-1)+ expEtotself(i)*vpolIm(t)*vsol 
!                 expsqrgrad = expsqrgradpsi(i)*vIm(1)      ! no mutipilcation with vsol because defintion constqE
!                 exppi(i,t) = (xsol(i)**vpolAA(1))*exp(psi(i)+expsqrgrad+expborn)   ! auxilary variable palpha
!             enddo  
!         else

!             do i=1,n

!                 expborn    = expEtotself(i)*vIm(t)*vsol 
!                 expsqrgrad = expsqrgradpsi(i)*vIm(t)
!                 exppi(i,t) = (xsol(i)**vpol(t))*exp(expsqrgrad+expborn)
!             enddo  
!         endif   
!     enddo               

!     ! Van der Waals   
!     ! if(geometry=="planar") then

!     !     do t=1, nsegtypes ! symmetrize 
!     !         do i=1,range
!     !             rhopolin(1-i,t)  =rhopolin(i,t)
!     !             rhopolin(nr+i,t) =0.0_dp
!     !         enddo
!     !     enddo    
    
!     !     do t=1,nsegtypes
!     !         if(isrhoselfconsistent(t)) call VdW_contribution_exp_planar(rhopolin,exppi(:,t),t)
!     !     enddo

!     ! else

!     !     do t=1,nsegtypes 
!     !         if(isrhoselfconsistent(t)) call VdW_contribution_exp_curved(rhopolin,exppi(:,t),t)
!     !     enddo
    
!     ! endif    
    
!     do t=1, nsegtypes
!         do i=1,ln_node
!             do j=1,nsize
!                 npol(i,j,t)=0.0_dp
!             enddo
!         enddo
!     enddo 

!     !  .. computation ionic liquid volume fraction 

!     do ln=1,nsize+ln_extra     ! loop starting layers  

!         l=layeroffset+ln

!         q(ln) = 0.0_dp         ! init q
     
!         do c=1,cuantas         ! loop over cuantas

!             pro=1.0_dp     

!             do s=1,nseg        ! loop over segments 
!                 k=indexchain(s,c,ln)     
!                 t=type_of_monomer(s)
! 				if (k<=nsize) then ! probability distribution fnc 
!                     pro = pro*exppi(k,t)
!                	else
!                     pro = pro*exppibulk(k,t)
!                	endif
!             enddo    
            
!             q(ln) = q(ln)+pro

!             do s=1,nseg
!                 k=indexchain(s,c,ln) 
!                 t=type_of_monomer(s)
!                 nIm(ln,k,t)=nIm(ln,k,t)+pro ! unnormed ionic density at k given that the 'begining'of chain is at l
!             enddo
         
!         enddo
      
!     enddo

!     !     .. construction of density ionic liquid     
!     !     .. volume fraction unnormed expmu need to be determined yet
      
!     do t=1,nsegtypes   
!         do i=1,n           
!             do ln=1,nsize+ln_extra 
!                 k=layeroffset+ln
!                 rhoIm(i,t)=rhoIm(i,t) + (deltaG(k)/deltaG(i))*nIm(ln,i,t)
!             enddo
!         enddo
!     enddo    

!     !     .. construction of fcn and volume fraction IM    
                                                  
!     do t=1, nsegtypes
!         do i=1,n
!             rhoIm(i,t) = expmuIm * rhoIm(i,t) ! density ionice of type t 
!             xpol(i)=xpol(i)+rhopol(i,t)*vpolIm(t)*vsol  
!         enddo    
!     enddo
        
!     do i=1,n 
!        rhoqIm(i)  = rhopol(i,tA)*zIm  ! charge density ionic liquid assume only one type of charged ionic liquid
!        phiIm(i)   = q(i)*expmuIm
!     enddo
      
! end subroutine density_il_conf

end module
