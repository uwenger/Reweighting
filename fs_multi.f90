!     ifort -r8 -o fs_multi fs_multi.f90 ranlxd_generator.f90 fs_multi_par.f90 -L ~/Libraries/NAG_Mark19/nag/fldau19da/ -lnag 
!     18/09/2025: to get rid of the warnings:
!     ifort -r8 -o fs_multi fs_multi.f90 fs_multi_par.f90 ranlxd_generator.f90 -L ~/Libraries/NAG_Mark19/nag/fldau19da/ -lnag -Wl,-ld_classic
!     15/06/2026: with nbeta > 9 there is a segmentation fault due to the stacksize hard limit of 63MB on macOS. To overcome use -heap-arrays to automatically push all automatic and temporary arrays onto the heap instead of stack:
!     ifort -r8 -o fs_multi fs_multi.f90 fs_multi_par.f90 ranlxd_generator.f90 -L ~/Libraries/NAG_Mark19/nag/fldau19da/ -lnag -Wl,-ld_classic  -heap-arrays


      program fs_multi

!     Ferrenberg-Swendson multi-histogram interpolation method.
!     (Note: this version uses binning of the data.)

      use fs_multi_par

      implicit none


      
!     number of  energy bins and measurements:
      integer :: nbin,nmeas(0:max_nbeta)
!     the bin energies and the beta value:
      real(kind=dp) :: bin_energy(0:max_nbin),beta_pr

      ! the values of the obsrvable bins and the probablilites:
      real(kind=dp) :: obsbin_val(0:obs_nbin),obs_prob(0:obs_nbin,max_nbeta)

      !     for the observable probability distribution:
      real(kind=dp) :: obs_en_prob(0:obs_nbin,0:max_nbin,0:max_nbeta)

      
!     the spectral density:
      real(kind=dp) :: en_dens(0:max_nbin,0:nbtrp)

!     the probability distributions:
      real(kind=dp) :: new_en_prob(0:max_nbin,0:nbtrp)

!     the effective observable at action S:
      real(kind=dp) :: eff_obs(0:max_nbin,0:nbtrp),eff_obs2(0:max_nbin,0:nbtrp)
!     the new observable and its error:
      real(kind=dp) :: obs(0:nbtrp),obs2(0:nbtrp)

!     the input file name:
      character*70 in_fname
!     the width of the bin:
      real(kind=dp) :: bin_width 
!     the beta values of the MC runs and the average beta_k(0):
      integer nbeta,ndbeta
      real(kind=dp) :: beta(0:max_nbeta),dbeta

!     for the bootstrap:
      real(kind=dp) :: btrp_av,btrp_err
      integer :: ibtrp

!     the estimated integrated autocorrelation time 
!     (used for calculating the error):
      real(kind=dp) :: tau_avg

!     the average
      real(kind=dp) :: act_av(0:max_nbeta)

!     which quantity to calculate: L -> abs. value of the Polyakov loop
!                                  S -> Pol. loop susceptibility
      character*1 quantity

!     the average beta value
      real(kind=dp) :: beta0
      
!     for estimating the max. of the susceptibility:
      real(kind=dp) :: chi_max,beta_max

!     the bias corrected critical beta:
      real(kind=dp) :: beta_corr

!     the lattice size:
      integer lsize,tsize
      real(kind=dp) :: vol

!     the energy for the ratio of weights:
      real(kind=dp) :: e0
!     and the number of deconfined phases:
      integer nc

!     auxiliary variables:
      real(kind=dp) :: e2, e4, aux, aux2, pmax1, pmax2
      integer i,j,nfail,k
      logical :: equal_prob_found = .False.
      real(kind=dp) :: beta_equal_prob

!     reading the data:
!--------------------------------------------------------------------
!     read in the input parameters:
      open(55,file='in_fs_multi.dat',form='formatted',status='old')
      read(55,*) in_fname
      print *,'Input file name:'
      print *,in_fname
      read(55,*) lsize,tsize
      print *,'Lattice size:',lsize,tsize
      read(55,*) bin_width
      print *,'Width of the bin:',bin_width
      read(55,*) nbeta
      print *,'Number of beta values:',nbeta
      beta(0)=0.
      do i=1,nbeta
         read(55,*) beta(i)
         beta(0)=beta(0)+beta(i)
         print '("Beta value",i3,":",f10.5)',i,beta(i)
      enddo
      beta(0)=beta(0)/float(nbeta)
      beta0=beta(0)
      print *,'Average beta value:',beta(0)
      read(55,*) dbeta
      print *,'delta beta:',dbeta
      read(55,*) ndbeta
      print *,'number of new beta values:',2*ndbeta+1
      print *
      read(55,*) beta(0)
      print *,'starting beta value:',beta(0)
      beta0=beta(0)
      print *
      read(55,*) quantity
      print *,'Quantity to calculate:',quantity
      if(quantity.eq.'R' .or. quantity.eq.'r') then
         read(55,*) e0
         read(55,*) nc
         print *,'E_0 for ratio of weights:',e0
         print *,'Number of deconfined phases:',nc
      endif
      close(55)

!     get the spectral densities wbar -> en_dens(S,0:nbtrp),
!     where in the last index label 0 denotes the original distribution
!     and 1..nbtrp distributions used for the bootstrap:
!----------------------------------------------------------------------
      print *
      print *,'Generating the energy distributions:'
      print *,'************************************'

      call gen_en_distr(nbin,nmeas,bin_energy,nbeta,beta,& 
           act_av,en_dens,eff_obs,eff_obs2,tau_avg)

      if(quantity.eq.'L' .or. quantity.eq.'l') then
         print *
         print *,'The abs. value of the Polyakov loop:'
         print *,'beta     |L|(beta)    bootstrap error'
      elseif(quantity.eq.'S' .or. quantity.eq.'s') then
         print *
         print *,'The Polyakov loop susceptibility:'
         print *,'beta     Chi(beta)    bootstrap error'
      elseif(quantity.eq.'C' .or. quantity.eq.'c') then
         print *
         print *,'The specific heat C = beta^2 V (<e^2> - <e>^2) with V=L^3 Lt:'
         print *,'beta     C(beta)    bootstrap error'
      elseif(quantity.eq.'B' .or. quantity.eq.'b') then
         print *
         print *,'The Binder cumulant B = 1 - <e^4> / (3 <e^2>^2):'
         print *,'beta     B(beta)    bootstrap error'
      elseif(quantity.eq.'R' .or. quantity.eq.'r') then
         print *
         print *,'The ratio of weight'
         print *,'         R = ln sum_{e<e0} P(e) / sum_{e>=e0} P(e):'
         print *,'beta     R(beta)    bootstrap error'
      elseif(quantity.eq.'E' .or. quantity.eq.'e') then
         print *
         print *,'The internal energy <e>:'
         print *,'beta     e(beta)    bootstrap error'
      elseif(quantity.eq.'H' .or. quantity.eq.'h') then
         print *
         print *,'The specific heat H = beta^2 V (<e^2> - <e>^2) with V=L^3 Lt:'
         print *,'beta     H(beta)    bootstrap error'
      endif

      open(47,file='obs_beta_pr.plo',form='formatted',status='unknown')
!     loop over new beta values:
      chi_max=0.
      do i=-ndbeta,ndbeta
         beta_pr=beta(0)+i*dbeta

!     get the energy probability distributions for the new beta value:
         call gen_new_distr(nbin,bin_energy,act_av,i*dbeta,en_dens,new_en_prob)

!     get the observable at the new beta value:
         call get_obs(nbin,new_en_prob,eff_obs,eff_obs2,obs,obs2)

         open(33,file='new_en_prob.plo',status='unknown',form='formatted')
         do j=0,nbin
            write(33,*) bin_energy(j)+act_av(0),new_en_prob(j,0)
         enddo
         close(33)

!     calculate the observable:
!------------------------------
         if(quantity.eq.'S' .or. quantity.eq.'s') then
!     calculate the Polyakov loop susceptibility
!     chi_pol=(<|P|^2> - <|P|>^2):
            do j=0,nbtrp
               obs(j)=lsize**3*(obs2(j)-obs(j)**2)
!               print *,'j,obs(j)=',j,obs(j)
            enddo
!     find a crude approximation of the max. of the susceptibility:
            if(obs(0).gt.chi_max) then
               beta_max=beta_pr
               chi_max=obs(0)
               beta_max=beta_max-beta(0)
!     print *,'beta(0)=',beta(0)
!     print *,'beta_max=',beta_max
            endif

            ! Check for the maxima and the minimum:
            call check_obs_distr_maxima(obs_en_prob,en_dens,nbin,obsbin_val,lsize,tsize,beta_pr,pmax1,pmax2)
            if(pmax1 < pmax2 .and. .not.equal_prob_found) then
               print *,'Equal prob found:'
               print *,beta_pr,pmax1, pmax2
               beta_equal_prob = beta_pr
               equal_prob_found = .True.
            endif
         endif  
         if(quantity.eq.'C' .or. quantity.eq.'c') then
!     calculate the specific heat
!     C=beta^2 lsize^3 tsize (<e^2> - <e>^2):
            do j=0,nbtrp
               obs(j)=beta_pr**2*lsize**3*tsize*(obs2(j)-obs(j)**2)
            enddo
!     find a crude approximation of the max. of the specific heat:
            if(obs(0).gt.chi_max) then
               beta_max=beta_pr
               chi_max=obs(0)
               beta_max=beta_max-beta(0)
     print *,'DEBUG: beta(0)=',beta(0)
     print *,'DEBUG: beta_max=',beta_max
            endif
         endif  
         if(quantity.eq.'B' .or. quantity.eq.'b') then
!     calculate the Binder cumulant
!     B = 1 - <e^4> / (3. * <e^2>^2):
!$$$            do j=0,nbtrp
!$$$               obs(j)=1.-obs2(j)/(3.*obs(j)**2)
!$$$            enddo

            do j=0,nbtrp
               e2=0.0
               e4=0.0
               do k=0,nbin
                     e4=e4+new_en_prob(k,j)*  &
                      ((bin_energy(k)+act_av(0))/ &
                      float(tsize*lsize**3))**4
                     e2=e2+new_en_prob(k,j)*  &
                      ((bin_energy(k)+act_av(0))/ &
                         float(tsize*lsize**3))**2
               enddo
                obs(j)=lsize**3*(e4/(e2**2)-1.)
            enddo
!     find a crude approximation of the max. of the Binder cumulant:
            if(obs(0).gt.chi_max) then
               beta_max=beta_pr
               chi_max=obs(0)
               beta_max=beta_max-beta(0)
!     print *,'beta(0)=',beta(0)
!     print *,'beta_max=',beta_max
            endif
         endif  
         if(quantity.eq.'R' .or. quantity.eq.'r') then
!     calculate the ratio of weight:
!     R = sum_{E<E_0} P(E) / sum_{E>=E_0} P(E) 
            do j=0,nbtrp
               e2=1.0e-8
               e4=1.0e-8
               do k=0,nbin
!                  if((bin_energy(k)+act_av(0)).lt.55450) then
!                  if(bin_energy(k).lt.0.0) then
                  if((bin_energy(k)+act_av(0))/lsize**3.lt.e0) then
                     e2=e2+new_en_prob(k,j)
                  else
                     e4=e4+new_en_prob(k,j)
                  endif
               enddo
               obs(j)=log(e2/e4)
            enddo
!     find a crude approximation of the solution to R = q :
!$$$            if((1.-(obs(0) - log(4.))**2).gt.chi_max) then
            if((1.-(obs(0) - log(float(nc)))**2).gt.chi_max) then
               beta_max=beta_pr
!$$$               chi_max=(1.-(obs(0) - log(4.))**2)
               chi_max=(1.-(obs(0) - log(float(nc)))**2)
               beta_max=beta_max-beta(0)
!            print *,'beta(0)=',beta(0)
!            print *,'beta_max=',beta_max
            endif
         endif
         if(quantity.eq.'E' .or. quantity.eq.'e') then
!     calculate the internal energy:
!     e = E / V = sum_{E} P(E) * E / V
            do j=0,nbtrp
               chi_max=0.0
               beta_max=0.0
               do k=0,nbin
                     chi_max=chi_max+new_en_prob(k,j)*bin_energy(k)
               enddo
               obs(j)=(chi_max+act_av(0))/float(lsize**3)
            enddo
         endif
         if(quantity.eq.'H' .or. quantity.eq.'h') then
!     calculate the specific heat again:
!     H = 
            do j=0,nbtrp
               aux=0.0
               aux2=0.0
               do k=0,nbin
!     <e> = <E/V>:
                     aux=aux+new_en_prob(k,j)* &
                      (bin_energy(k))/float(tsize*lsize**3)
!                     chi_max=chi_max+new_en_prob(k,j)*
!     &                 (bin_energy(k))
!     <e^2> = <(E/V)^2>:
                     aux2=aux2+new_en_prob(k,j)* &
                      ((bin_energy(k))/ &
                         float(tsize*lsize**3))**2
!                     beta_max=beta_max+new_en_prob(k,j)*
!     &                 (bin_energy(k))**2
               enddo
!                obs(j)=lsize**3*(beta_max-chi_max**2)*beta_pr**2
!                obs(j)=beta_pr**2*(beta_max-chi_max**2)*tsize/float(6)
                obs(j)=beta_pr**2*(aux2-aux**2)*lsize**3*tsize
!               obs(j)=(beta_max-chi_max**2)*beta_pr**2/
!     &              float(tsize**2*lsize**3)
            enddo
!     find a crude approximation of the max. of the specific heat:
            if(obs(0).gt.chi_max) then
               chi_max=obs(0)
               beta_max=beta_pr-beta(0)
            endif
            ! Check for the maxima and the minimum:
            call check_obs_distr_maxima(obs_en_prob,en_dens,nbin,obsbin_val,lsize,tsize,beta_pr,pmax1,pmax2)
            if(pmax1 < pmax2 .and. .not.equal_prob_found) then
               print *,'Equal prob found:'
               print *,beta_pr,pmax1, pmax2
               beta_equal_prob = beta_pr
               equal_prob_found = .True.
            endif
 
         endif

!     and its bootstrap error:
!------------------------------
         btrp_av=0.
         do j=1,nbtrp
            btrp_av=btrp_av+obs(j)
         enddo
         btrp_av=btrp_av/float(nbtrp)   !bootstrap average
         btrp_err=0.
         do j=1,nbtrp
            btrp_err=btrp_err+(obs(j)-btrp_av)**2
         enddo
!     for uncorrelated data:
         btrp_err=sqrt(btrp_err/float(nbtrp-1))
!     multiply the error by 2 tau due to autocorrelation:
!         btrp_err=2.*tau_avg*btrp_err/float(nbtrp)
!         btrp_err=sqrt(btrp_err)

!     write it out:
!--------------------------------------------
!!$         print *,beta_pr,obs(0),btrp_err
         if(ndbeta.eq.0) then
            write(47,*) beta_pr,obs(0),btrp_err
         else
            write(47,'(f10.5,3f18.6)') beta_pr,obs(0),obs(0)+btrp_err,obs(0)-btrp_err
         endif
      enddo
      close(47)
      if (.not.(quantity.eq.'L' .or. quantity.eq.'l')) then
         write(*,'("Crude estimate of beta_c=",f15.6)') beta_max+beta(0)
      endif
      
      if(quantity.eq.'S' .or. quantity.eq.'s') then
!     get the location of the maximum of the susceptibility:
         !-----------------------------------------------------------
         vol=(float(lsize)/float(tsize))**3
         call get_betac(beta_max,btrp_err,beta_corr,beta,nfail)
!     for correlated data multiply the error by 2 tau due 
!     to autocorrelation:
!         btrp_err=btrp_err**2
!         btrp_err=2.*tau_avg*btrp_err*(nbtrp-1)/float(nbtrp)
         print *
         print *,'Location of the maximum of the susceptibility:'
         print *,'**********************************************'
         if(nfail.ne.0) then
            print '("Attention: nfail=",i3)',nfail
         endif
         print '(2f12.7," +/- ",f12.7," beta_c(S_max)")',1.0_dp/vol,beta_max,btrp_err
         print '(2f12.7," +/- ",f12.7,"              ")',1.0_dp/vol,beta_corr,btrp_err
         print *
!     get the energy probability distributions for the new beta value:
         dbeta = beta_max-beta(0)
         call gen_new_distr(nbin,bin_energy,act_av,dbeta,en_dens,new_en_prob)
!     get the observable at the new beta value:
         call get_obs(nbin,new_en_prob,eff_obs,eff_obs2,obs,obs2)
!     calculate the Polyakov loop susceptibility
!     chi_pol=(<|P|^2> - <|P|>^2):
            do j=0,nbtrp
               obs(j)=lsize**3*(obs2(j)-obs(j)**2)
            enddo
            btrp_av=0.
            do j=1,nbtrp
               btrp_av=btrp_av+obs(j)
            enddo
            btrp_av=btrp_av/float(nbtrp) !bootstrap average
            btrp_err=0.
            do j=1,nbtrp
               btrp_err=btrp_err+(obs(j)-btrp_av)**2
            enddo
            btrp_err=sqrt(btrp_err/float(nbtrp-1))
            

         print *,'Susceptibility maximum S_max and S_max/V:'
         print *,'***********************************************'
         print '(2f12.7," +/- ",f12.7," S_max(beta_c)")',1.0_dp/vol,obs(0),btrp_err
         
         print '(2f12.7," +/- ",f12.7," S_max_over_V")',1.0_dp/vol,obs(0)/vol,btrp_err/vol
         print *
         ! Calculate the surface tension:
!!$         call print_obs_distr(obs_en_prob,en_dens,nbin,obsbin_val,lsize,tsize)

         ! Calculate the surface tension at equal probability height:
         dbeta = beta_equal_prob-beta(0)
         call gen_new_distr(nbin,bin_energy,act_av,dbeta,en_dens,new_en_prob)
         call print_obs_distr(obs_en_prob,en_dens,nbin,obsbin_val,lsize,tsize)

         
      endif
      if(quantity.eq.'B' .or. quantity.eq.'b') then
!     get the location of the maximum of the susceptibility:
!-----------------------------------------------------------
         call get_betac(beta_max,btrp_err,beta_corr,beta,nfail)
!     for correlated data multiply the error by 2 tau due 
!     to autocorrelation:
!         btrp_err=btrp_err**2
!         btrp_err=2.*tau_avg*btrp_err*(nbtrp-1)/float(nbtrp)
         print *
         print *,'Location of the maximum of the Binder cumulant:'
         print *,'***********************************************'
         if(nfail.ne.0) then
            print '("Attention: nfail=",i3)',nfail
         endif
         print '(f12.7," +/- ",f12.7)',beta_max,btrp_err
         print '(f12.7," +/- ",f12.7)',beta_corr,btrp_err
         print *
      endif
      if(quantity.eq.'C' .or. quantity.eq.'c' .or. quantity.eq.'H' &
           .or. quantity.eq.'h') then
         vol=(float(lsize)/float(tsize))**3
!     get the location of the maximum of the susceptibility:
!-----------------------------------------------------------
         call get_betac(beta_max,btrp_err,beta_corr,beta,nfail)
!     for correlated data multiply the error by 2 tau due 
!     to autocorrelation:
!         btrp_err=btrp_err**2
!         btrp_err=2.*tau_avg*btrp_err*(nbtrp-1)/float(nbtrp)
         print *
         print *,'Location of the maximum of the specific heat:'
         print *,'***********************************************'
         if(nfail.ne.0) then
            print '("Attention: nfail=",i3)',nfail
         endif
         print '(2f12.7," +/- ",f12.7," beta_c(C_max)")',1.0_dp/vol,beta_max,btrp_err
         print '(2f12.7," +/- ",f12.7," bias corrected")',1.0_dp/vol,beta_corr,btrp_err
         print *
!     get the energy probability distributions for the new beta value:
         dbeta = beta_max-beta(0)
         call gen_new_distr(nbin,bin_energy,act_av,dbeta,&
             en_dens,new_en_prob)

            do j=0,nbtrp
               aux=0.0
               aux2=0.0
               do k=0,nbin
!     <e> = <E/V>:
                  aux=aux+new_en_prob(k,j)* &
                      (bin_energy(k))/float(tsize*lsize**3)
!     <e^2> = <(E/V)^2>:
                  aux2=aux2+new_en_prob(k,j)* &
                      ((bin_energy(k))/ &
                      float(tsize*lsize**3))**2
               enddo
               obs(j)=beta_max**2*(aux2-aux**2)*lsize**3*tsize
            enddo
            btrp_av=0.
            do j=1,nbtrp
               btrp_av=btrp_av+obs(j)
            enddo
            btrp_av=btrp_av/float(nbtrp) !bootstrap average
            btrp_err=0.
            do j=1,nbtrp
               btrp_err=btrp_err+(obs(j)-btrp_av)**2
            enddo
            btrp_err=sqrt(btrp_err/float(nbtrp-1))
            
            
         print *,'Specific heat maximum C_max and C_max / V:'
         print *,'***********************************************'
         print '(2f12.7," +/- ",f12.7," C_max")',1.0_dp/vol,obs(0),btrp_err
         vol=(float(lsize)/float(tsize))**3
         print '(2f12.7," +/- ",f12.7," C_max/V")',1.0_dp/vol,obs(0)/vol,btrp_err/vol
         print *

         ! C_max/(beta**2*L**3*Lt) = L->oo => 1/4 (<s_c> - <s_d>)^2 = 1/4 L_h^2
         print *,'beta_max=',beta_max
         do j=0,nbtrp
            obs(j)=4.0*obs(j)/(beta_max**2*lsize**3*tsize)
            obs(j)=sqrt(obs(j))
         enddo
         btrp_av=0.
         do j=1,nbtrp
            btrp_av=btrp_av+obs(j)
         enddo
         btrp_av=btrp_av/float(nbtrp) !bootstrap average
         btrp_err=0.
         do j=1,nbtrp
            btrp_err=btrp_err+(obs(j)-btrp_av)**2
         enddo
         btrp_err=sqrt(btrp_err/float(nbtrp-1))
            
         print *,'Latent heat (from C_max):'
         print *,'***********************************************'
         print '(2f12.7," +/- ",f12.7," L_h")',1.0_dp/vol,obs(0),btrp_err
         print '(2f12.7," +/- ",f12.7," L_h/Tc^4")',1.0_dp/vol,tsize**4*obs(0),tsize**4*btrp_err

         ! Calculate the surface tension at the specific heat peak:
         call print_obs_distr(obs_en_prob,en_dens,nbin,obsbin_val,lsize,tsize)

!!$         ! Calculate the surface tension at equal probability height:
!!$         dbeta = beta_equal_prob-beta(0)
!!$         call gen_new_distr(nbin,bin_energy,act_av,dbeta,en_dens,new_en_prob)
!!$         call print_obs_distr(obs_en_prob,en_dens,nbin,obsbin_val,lsize,tsize)

 
      endif
      if(quantity.eq.'R' .or. quantity.eq.'r') then
!     get the location of the crossing of the ratio with the number of phases:
!-----------------------------------------------------------
         call get_betac(beta_max,btrp_err,beta_corr,beta,nfail)
         print *
!$$$         print *,'Solution of R = ln N with N = 4:'
         print '(" Solution of R = ln N with N =",i2,":")',nc
         print '(" ***********************************************")'
         if(nfail.ne.0) then
            print '("Attention: nfail=",i3)',nfail
         endif
         print '(f12.7," +/- ",f12.7)',beta_max,btrp_err
         print '(f12.7," +/- ",f12.7)',beta_corr,btrp_err
         print *


!     get the energy probability distributions for the new beta value:
         dbeta = beta_max-beta(0)
         call gen_new_distr(nbin,bin_energy,act_av,dbeta,en_dens,new_en_prob)

            do j=0,nbtrp
               e2=1.0e-8
               e4=1.0e-8
               chi_max=1.0e-8
               beta_max= 1.0e-8
               do k=0,nbin
                  aux=(bin_energy(k)+act_av(0))/float(lsize**3)
                  aux2=bin_energy(k)/float(lsize**3*tsize)
                  if(aux.lt.e0) then
                     e2=e2+new_en_prob(k,j)
                     chi_max = chi_max + aux2*new_en_prob(k,j)
                  else
                     e4=e4+new_en_prob(k,j)
                     beta_max = beta_max + aux2*new_en_prob(k,j)
                  endif
               enddo
               obs(j)=beta_max/e4 - chi_max/e2
            enddo
            btrp_av=0.
            do j=1,nbtrp
               btrp_av=btrp_av+obs(j)
            enddo
            btrp_av=btrp_av/float(nbtrp) !bootstrap average
            btrp_err=0.
            do j=1,nbtrp
               btrp_err=btrp_err+(obs(j)-btrp_av)**2
            enddo
            btrp_err=sqrt(btrp_err/float(nbtrp-1))
            
            
         print *,'Latent heat Delta e:'
         print *,'***********************************************'
         print '(f12.7," +/- ",f12.7)',obs(0),btrp_err
         print *
         

      endif
         
      print *,'Program ended!'

      stop
  
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!--------------------------------------------------------------------------
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
contains

      subroutine get_betac(bc,err_bc,bc_corr,beta,nfail)

!     This subroutine calculates the location of the maximum of the 
!     reweighted susceptibility using a NAG-library minimization routine.
!     Note, that 'bc' is actually the difference to the 'average' beta
!     value 'beta(0)'.
      
      use fs_multi_par

      implicit none

!     bootstrap sample number:
      integer :: nb

!     variables for the NAG-routine (in double precision):
      real(kind=dp) :: chi_max,chi_max_der,e1,e2,low,up,dbc,bound
      integer max_cal,ifail
!      external chi_beta

!     the bootstrap maximum locations:
      real(kind=dp) :: beta_c(0:nbtrp)

!     the critical beta and its error:
      real(kind=dp), intent(out) :: bc,err_bc
!     the bias corrected critical beta:
      real(kind=dp), intent(out) :: bc_corr
      
!     the beta values of the MC runs and the average beta_k(0):
      real(kind=dp), intent(in) :: beta(0:max_nbeta)

      !     Number of fails to locate beta_c:
      integer, intent(out) :: nfail
      
!     auxiliary variables:
      real(kind=dp) :: btrp_av,btrp_err
      integer :: nbound

      nfail=0
      nbound=0

!      print *,'DEBUG: beta(0), bc, dble(bc)=',beta(0), bc, dble(bc)

      
!     loop over the bootstrap samples:
      do nb=0,nbtrp
         ibtrp=nb ! store this globally to make it available in subroutine chi_beta
         !     convert to double precision:
         dbc=dble(bc)
!     calculate the maximum of the susceptibility:
!-------------------------------------------------
!     the relative and absolute accuracy:
         e1=1.d-8
         e2=1.d-8
         !     lower and upper bound of search interval:
         bound=0.01
         low=dbc-bound
         up=dbc+bound
         !     maximal number of calls:
         max_cal=2000
         ifail=-1
         !      print *,'Calling E04BBF:'
         call E04BBF(chi_beta,e1,e2,low,up,max_cal,dbc,chi_max,chi_max_der,ifail)
         chi_max=-chi_max
         chi_max_der=-chi_max_der

         !     check if the solution is at the boundary of the search interval:
         if(dbc.le.(bc-bound) .or. dbc.ge.(bc+bound)) nbound=nbound+1
         
         beta_c(nb)=real(dbc)+beta(0)
!         print *,'DEBUG: nb, beta_c(nb)=', nb, beta_c(nb)
         if(ifail.ne.0) then
            print *,'**********************'
            print *,'Maximization failed!'
            print *,'nb=',nb
            print *,'ifail=',ifail
            print *,'**********************'
            print *,'low,up=',low,up
            print *,'max_cal=',max_cal
            print *,'chi_max_der=',chi_max_der
            nfail=nfail+1
            beta_c(nb)=0.
         endif
         if((ifail.ne.0) .and. (nb==0)) then
            print *,'**********************'
            print *,'Maximization failed!'
            print *,'on original data set!'
            print *,'**********************'
            print *,'ifail=',ifail
            print *,'low,up=',low,up
            print *,'max_cal=',max_cal
            print *,'chi_max_der=',chi_max_der
         endif
!      print *,'max_cal,beta_c=',max_cal,beta_c(nb)

!     end of loop over the bootstrap samples:
      enddo

!     print out error message if the bound of the search interval 
!     was reached for some bootstrap samples:
      if(nbound.ne.0) then
         print *
         print *,'****************************************'
         print *,'Maximization reached search boundary'
         print *,'in number of samples, nbound=',nbound
         print *,'****************************************'
      endif


!     calculate the error from the bootstrap samples:
      btrp_av=0.
      do nb=1,nbtrp
         btrp_av=btrp_av+beta_c(nb)
      enddo
      btrp_av=btrp_av/float(nbtrp-nfail)   !bootstrap average
      btrp_err=0.
      do nb=1,nbtrp
         if(beta_c(nb).ne.0.) then
            btrp_err=btrp_err+(beta_c(nb)-btrp_av)**2
         endif
      enddo
!     for uncorrelated data:
      btrp_err=sqrt(btrp_err/float(nbtrp-nfail-1))

!     for checking the error:
      open(10,file='betac.plo',form='formatted',status='unknown')
      do nb=0,nbtrp
         write(10,'(i6,3f20.10)') nb,beta_c(nb),bc+beta(0)+bound,&
             bc+beta(0)-bound
      enddo
      write(10,*) '&'
      write(10,'(i6,3f20.10)')  0,btrp_av,bc+beta(0),btrp_av
      write(10,'(i6,3f20.10)')  nbtrp,btrp_av,bc+beta(0),btrp_av
      write(10,*) '&'
      write(10,'(i6,3f20.10)')  0,beta_c(0),beta_c(0)+btrp_err,beta_c(0)-btrp_err
      write(10,'(i6,3f20.10)')  nbtrp,beta_c(0),beta_c(0)+btrp_err,beta_c(0)-btrp_err


      close(10)


      bc=beta_c(0)
      err_bc=btrp_err
!     the bias corrected critical beta:
!     bias = bc - btrp_av ==> bc_corr = bc + bias
      bc_corr=2*bc-btrp_av


      return
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!--------------------------------------------------------------------------
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine chi_beta(dbeta,chi,chi_der)

!     This function calculates the susceptibility chi and its derivative at 
!     beta=beta(0)+dbeta. It is called by the minimization routine E04BBF, 
!     thus chi -> -chi, chi_der -> -chi_der.

      use fs_multi_par

      implicit none

      real(kind=dp), intent(in) :: dbeta
      real(kind=dp), intent(out) :: chi,chi_der

      ! The following are globally defined:
      !------------------------------------
!     number of bins:
!      integer nbin
!     energy of the bins:
!      real(kind=dp) :: bin_energy(0:max_nbin)
!     the distributions:
!      real(kind=dp) :: en_dens(0:max_nbin,0:nbtrp),new_en_prob(0:max_nbin,0:nbtrp)

!     the effective observable at action S:
!      real(kind=dp) :: eff_obs(0:max_nbin,0:nbtrp),eff_obs2(0:max_nbin,0:nbtrp)

!     the energy averages:
      real(kind=dp) :: act_av(0:max_nbeta)

!     the energy for the ratio of weights:
      real(kind=dp) :: e0
      integer :: nc

!     bootstrap sample number:
      integer :: nb

!!$!     the quantity to calculate:
!!$      character*1 :: quantity

!!$!     the average beta value
!!$      real(kind=dp) :: beta

!!$!     the lattice size:
!!$      integer :: lsize,tsize

!     auxiliary variables:
      integer :: i
      real(kind=dp) :: aux,L,L2,L2S,LS,S


      nb = ibtrp
!     calculate the normalization factor Sum_S Wbar(S) exp(-dbeta*dS):
!      print *,'Calculating the normalization:'
      aux=0.0_dp
!      print *,'nbin=',nbin
!      print *,'dbeta=',dbeta
      do i=0,nbin
         aux=aux+en_dens(i,nb)*exp(-dbeta*bin_energy(i))
!         print *,'DEBUG: en_dens, bin_energy=',i,en_dens(i,nb), bin_energy(i)
      enddo
!      print *,'DEBUG: aux=',aux


      
!     calculate the new probability:
!      print *,'Calculating new prop.:'
      do i=0,nbin
         new_en_prob(i,nb)=en_dens(i,nb)*exp(-dbeta*bin_energy(i))/aux
      enddo

      L=0.
      L2=0.
      LS=0.
      L2S=0.
      S=0.
!      print *,'Quantity=',quantity
!     sum over the energy bins:
!!$      if (quantity.eq.'H' .or. quantity.eq.'h') then
!!$         do i=0,nbin
!!$            L=L+bin_energy(i)*new_en_prob(i,nb)
!!$            L2=L2+bin_energy(i)**2*new_en_prob(i,nb)
!!$            LS=LS+bin_energy(i)*bin_energy(i)*new_en_prob(i,nb)
!!$            L2S=L2S+bin_energy(i)**3*new_en_prob(i,nb)
!!$            S=S+bin_energy(i)*new_en_prob(i,nb)
!!$         enddo
!!$!     calculate the susceptibility:
!!$         chi=-(L2-L**2)*(beta+dbeta)**2/float(tsize**2*lsize**3)
!!$!     and its derivative:
!!$!         chi_der=-(-L2S+L2*S-2*L*(L*S-LS))*(beta+dbeta)**2/
!!$!     &        float(tsize**2*lsize**3)
!!$!         chi_der=-((L2-L**2)*2*(beta+dbeta)+
!!$!     &        (-L2S+2*L*LS)*(beta+dbeta)**2)/float(tsize**2*lsize**3)
!!$         chi_der=2*chi*(S+1.0/(beta+dbeta))- &
!!$             (L2*S-L2S)*(beta+dbeta)**2/float(tsize**2*lsize**3)
!!$
!!$
!!$!      if(nb==0) print '("dbeta,chi,chi_der=",3e20.10)',
!!$!     &       beta+dbeta,-chi,-chi_der
      if (quantity.eq.'H' .or. quantity.eq.'h' .or. quantity.eq.'C' .or. quantity.eq.'c') then
      !==================================================
         do i=0,nbin
            L=L+bin_energy(i)/float(tsize*lsize**3)*new_en_prob(i,nb)                        ! <s>
            L2=L2+(bin_energy(i)/float(tsize*lsize**3))**2*new_en_prob(i,nb)                 ! <s^2>
            LS=LS+bin_energy(i)/float(tsize*lsize**3)*bin_energy(i)*new_en_prob(i,nb)        ! <s S>
            L2S=L2S+(bin_energy(i)/float(tsize*lsize**3))**2*bin_energy(i)*new_en_prob(i,nb) ! <s^2 S>
            S=S+bin_energy(i)*new_en_prob(i,nb)                                              ! <S>
         enddo
!     calculate the susceptibility:
         chi = (L2-L**2)*(beta(0)+dbeta)**2*lsize**3*tsize  
!     and its derivative:
!         chi_der=-(-L2S+L2*S-2*L*(L*S-LS))*(beta(0)+dbeta)**2* &
!             float(tsize*lsize**3)
!         chi_der=-((L2-L**2)*2*(beta(0)+dbeta)+ &
!              (-L2S+2*L*LS)*(beta(0)+dbeta)**2)*float(tsize*lsize**3)
         chi_der= chi*(S+2.0/(beta(0)+dbeta)) - (beta(0)+dbeta)**2*float(lsize**3)*tsize * &
              (L2S - 2.0*L*LS + L**2*S) 
!!$         chi_der=  chi*2.0/(beta(0)+dbeta) + (beta(0)+dbeta)**2*float(lsize**3)*tsize * &
!!$              (L2*S - L2S - 2.0*(L**2*S-L*LS)) 

         if(nb==0) print '("beta,chi,chi_der=",3e20.10,2i4)', beta(0)+dbeta,chi,chi_der,lsize,tsize
         ! minus sign because we need to turn the maximum into a minimum
         chi = - chi
         chi_der = -chi_der
      elseif (quantity.eq.'B' .or. quantity.eq.'b') then
      !==================================================
         do i=0,nbin
            aux = (bin_energy(i)+act_av(0))/float(tsize*lsize**3)
            L=L+aux*new_en_prob(i,nb)
            L2=L2+aux**2*new_en_prob(i,nb)
            LS=LS+aux**4*new_en_prob(i,nb)
            L2S=L2S+aux**3*new_en_prob(i,nb)
            S=S+aux**5*new_en_prob(i,nb)
         enddo
!     calculate the susceptibility:
         chi=-(LS/(L2**2) - 1.)*float(lsize**3)
!     and its derivative:
         chi_der=(S/(L2**2) - LS*L2S/(L2**3))*float(lsize**3)**2
!      print *,'dbeta,chi,chi_der=',dbeta,chi,chi_der

      elseif(quantity.eq.'R' .or. quantity.eq.'r') then
      !==================================================
         do i=0,nbin
            LS=(bin_energy(i)+act_av(0))/lsize**3
            if(LS.lt.e0) then
               L=L+new_en_prob(i,nb)
               L2=L2+LS*new_en_prob(i,nb)
            else
               S=S+new_en_prob(i,nb)
               L2S=L2S+LS*new_en_prob(i,nb)
            endif
         enddo
!$$$         chi=(log(L/S) - log(4.))**2
!$$$         chi_der=-2.*(log(L/S) - log(4.))*lsize**3*(L2/L - L2S/S)
         chi=(log(L/S) - log(float(nc)))**2
         chi_der=-2.*(log(L/S) - log(float(nc)))*lsize**3*(L2/L - L2S/S)
      else  ! Polyakov loop susceptibility 'S'
      !==================================================
         do i=0,nbin
            L=L+eff_obs(i,nb)*new_en_prob(i,nb)
            L2=L2+eff_obs2(i,nb)*new_en_prob(i,nb)
            LS=LS+eff_obs(i,nb)*bin_energy(i)*new_en_prob(i,nb)
            L2S=L2S+eff_obs2(i,nb)*bin_energy(i)*new_en_prob(i,nb)
            S=S+bin_energy(i)*new_en_prob(i,nb)
         enddo
!     calculate the susceptibility:
         chi=-(L2-L**2)
!     and its derivative:
         chi_der=-(-L2S+L2*S-2*L*(L*S-LS))
!      if(nb==0) print '("dbeta,chi,chi_der=",3e20.10)',
         !     &       beta+dbeta,-lsize**3*chi,-chi_der
!         print *,'nb=',nb
         if(nb==0) print '("beta,chi,chi_der=",4e20.10,i4)', beta(0),beta(0)+dbeta,-chi,-chi_der,lsize

      endif

!      print *,'DEBUG in subroutine chi_beta: nb=', nb, dbeta, chi, chi_der

!      print *,'dbeta,chi,chi_der=',dbeta,chi,chi_der

      return
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!--------------------------------------------------------------------------
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine get_obs(nbin,new_en_prob,eff_obs,eff_obs2,obs,obs2)

      use fs_multi_par

      implicit none

!     number of bins:
      integer nbin

!     the new distribution:
      real(kind=dp) :: new_en_prob(0:max_nbin,0:nbtrp)

!     the effective observable at action S:
      real(kind=dp) :: eff_obs(0:max_nbin,0:nbtrp),eff_obs2(0:max_nbin,0:nbtrp)

!     the observable and its error:
      real(kind=dp) :: obs(0:nbtrp),obs2(0:nbtrp)

!     auxiliary variables:
      integer i,j

!     loop over the bootstrap samples:
      do j=0,nbtrp
      obs(j)=0.
      obs2(j)=0.
!     sum over the energy bins:
      do i=0,nbin
         obs(j)=obs(j)+eff_obs(i,j)*new_en_prob(i,j)
         obs2(j)=obs2(j)+eff_obs2(i,j)*new_en_prob(i,j)
!         obs2(j)=obs2(j)+eff_obs(i,j)**2*new_en_prob(i,j)
      enddo
      enddo


      return
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!--------------------------------------------------------------------------
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine gen_new_distr(nbin,bin_energy,act_av,dbeta,en_dens,new_en_prob)

      use fs_multi_par

      implicit none

!     number of bins:
      integer nbin

!     energy of the bins:
      real(kind=dp) :: bin_energy(0:nbin)
!     the average of the action:
      real(kind=dp) :: act_av(0:max_nbeta)

!     the difference beta value: 
      real(kind=dp) :: dbeta

!     the distributions:
      real(kind=dp) :: en_dens(0:max_nbin,0:nbtrp),new_en_prob(0:max_nbin,0:nbtrp)

!     auxiliary variables:
      integer i,nb
      real(kind=dp) :: aux,act

!     loop over the bootstrap samples:
      do nb=0,nbtrp

!     calculate the normalization factor Sum_S Wbar(S) exp(-dbeta*dS):
      aux=0.
      do i=0,nbin
         aux=aux+en_dens(i,nb)*exp(-dbeta*bin_energy(i))
      enddo

!     calculate the new probability:
      do i=0,nbin
         new_en_prob(i,nb)=en_dens(i,nb)*exp(-dbeta*bin_energy(i))/aux
      enddo

!     end of loop over the bootstrap samples:
      enddo

!     write out the new distribution:
!      print *,'Writing out the new distribution:'
      open(55,file='en_new_distr.plo',form='formatted',status='unknown')
      do i=0,nbin
         act=bin_energy(i)+act_av(0)
!     divide by the bin width for comparison:         
!         write(55,*) act, new_en_prob(i)/bin_width
         write(55,*) act,new_en_prob(i,0)
      enddo
      close(55)

      return
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!--------------------------------------------------------------------------
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine gen_en_distr(nbin,nmeas,bin_energy,nbeta,beta,act_av,en_dens,obs_bar,obs2_bar,tau_avg)

!     This routine calculates the energy probability distributions of 
!     MC runs with given beta values.
!     It needs a input file formatted like:
!
!     beta, value of the action, observable
!
!     The input file name is read in from the input data file 
!     'in_gen_en_distr.f', as well as the width of the energy bins 
!     'bin_width'.
!
!     Output is a file containing the energy distribution as a normalized 
!     histogram:
!     
!     average energy of the bin, probability

!     for ranlux generator:
      use ranlxd_generator
      use fs_multi_par

      implicit none

!     for ranlux generator:
      integer :: level


!     the input file name:
      character*70 in_fname

!     the width of the bin:
      real(kind=dp) :: bin_width 
!     energy of the bins:
      real(kind=dp) :: bin_energy(0:max_nbin)

!     number of beta values:
      integer nbeta

!     the free energies (guess values):
      real(kind=dp) :: fe(1:nbeta,0:nbtrp)

!     the energy averages:
      real(kind=dp) :: act_av(0:max_nbeta)

!     the beta values of the MC runs:
      real(kind=dp) :: beta(0:max_nbeta),beta_read
!     and the beta differences:
      real(kind=dp) :: dbeta(1:nbeta)

!     the action and the observable:
      real(kind=dp) :: act,obs,min_act(0:nbeta),max_act(0:nbeta)
      real(kind=dp) :: energy(1:max_meas,1:nbeta),obsrvb(1:max_meas,1:nbeta)

!     the distributions:
      real(kind=dp) :: en_prob(0:max_nbin,1:max_nbeta,0:nbtrp)

!     the spectral densities:
      real(kind=dp) :: en_dens(0:max_nbin,0:nbtrp)

!     the effective observable at action S:
      real(kind=dp) :: eff_obs(0:max_nbin,1:max_nbeta,0:nbtrp), &
          eff_obs2(0:max_nbin,1:max_nbeta,0:nbtrp)
      real(kind=dp) :: obs_bar(0:max_nbin,0:nbtrp),obs2_bar(0:max_nbin,0:nbtrp)

!     the errors on the distributions from the autocorrelation times:
      real(kind=dp) :: gtau(1:nbeta),tau(1:nbeta)

!     the estimated averaged integrated autocorrelation time 
!     (used for calculating the error):
      real(kind=dp) :: tau_avg

!     counts the number of measurements(nmeas(0) being the total number
!     of measurements):
      integer nmeas(0:max_nbeta)

!     bin number and number of bins:
      integer bin_number,nbin

!     for binning the observable:
      integer obsbin_number
      real(kind=dp) :: obsbin_width,obs_prob(0:obs_nbin,max_nbeta),min_obs,max_obs!   ,obsbin_val(0:obs_nbin)

!!$!     for the observable probability distribution:
!!$      real(kind=dp) :: obs_en_prob(0:obs_nbin,0:max_nbin,0:max_nbeta)


!     the spatial and temporal lattice size:
      integer lsize,tsize

!     the energy for the ratio of weights:
      real(kind=dp) :: e0
!     for the calculation of the ratio of weights:
      integer ir0, ir1
 
!     random variables:
      integer irvar,iseed
      real(kind=dp) :: rvar!,ranf

      real(kind=dp), dimension(1:1) :: r


!     which quantity to calculate: L -> abs. value of the Polyakov loop
!                                  S -> Pol. loop susceptibility
      character*1 quantity

!     auxiliary variables:
      integer i,k,j,nb,nblock,ndbeta
      real(kind=dp) :: av_obs(1:max_nbeta,0:nbtrp),av_obs2(1:max_nbeta,0:nbtrp),aux,aux2,norm,pk


!     read in the data file name again:
!-------------------------------------------------------------------
      open(55,file='in_fs_multi.dat',form='formatted',status='old')
      read(55,*) in_fname
      print *,'Data file name:'
      print *,in_fname
      read(55,*) lsize, tsize
      read(55,*) bin_width
      print *,'Width of the bin:',bin_width
      read(55,*) nbeta
      print *,'Number of beta values:',nbeta
      do i=1,nbeta
         read(55,*) aux
         print '("Beta value",i3,":",f10.5)',i,aux
      enddo
      read(55,*) aux
      print *,'delta beta:',aux
      read(55,*) i
      print *,'number of new beta values:',2*i+1
      read(55,*) aux
      print *,'starting beta value:',aux
      read(55,*) quantity
      print *,'Quantity to calculate:',quantity
      if(quantity.eq.'R' .or. quantity.eq.'r') then
         read(55,*) e0
      endif

      close(55)

!      stop

!     read in the data, count # of measurements, action average, 
!     min. and max. action value:
!-------------------------------------------------------------------
      act_av(0)=0.
      nmeas(0)=0
      min_act(0)=1.e+10
      max_act(0)=0.
      min_obs=1.e+10
      max_obs=0.
      do k=1,nbeta
         act_av(k)=0.
         av_obs(k,0)=0.
         av_obs2(k,0)=0.
         nmeas(k)=0
         min_act(k)=1.e+10
         max_act(k)=0.
         ir0=0
         ir1=0
         open(55,file=in_fname,form='formatted',status='old')
!         open(55,file=in_fname,form='unformatted',status='old')
         do i=1,max_nbeta*max_meas
!         do i=1,180000          !nbeta*max_meas
!            read(55,*,end=435) beta_read,act,aux,aux,obs  
!            read(55,*,end=435) beta_read,act,aux,aux,obs  
!            read(55,*,end=435) beta_read,act,aux,aux,obs  
!            read(55,*,end=435) beta_read,act,aux,aux,obs  
!!$            read(55,*,end=435) beta_read,act,aux,aux,obs
!!$!---------------------------------------------------
!!$!     this is for the APE action:
!!$            act=act/3._dp !divide the action by xi*nc
!$$$!---------------------------------------------------
!$$$!     this is for the Manton action (Philippe's version):
!$$$            read(55,*,end=435) beta_read,act,aux,obs
!$$$            act=(1.-act)*lsize**3*tsize*6 
!$$$            obs=sqrt(aux**2+obs**2)
!$$$!---------------------------------------------------
!$$$!     this is for the Wilson action:
!$$$            read(55,*,end=435) beta_read,act,aux,aux,obs
!$$$            act=(1.-act)*lsize**3*tsize*6 
!$$$            obs=sqrt(aux**2+obs**2)
!$$$!---------------------------------------------------
!$$$!     this is for the large N_c Wilson action and the QCDNCF90 Manton action:
!$$$            read(55,*,end=435) beta_read,act,aux,
!$$$     &           aux2,obs
!$$$!---------------------------------------------------
!$$$!            print '(f8.4,f20.8,3e20.7)',beta_read,act,aux,aux2,obs
!$$$            if (Quantity.eq.'C' .or. quantity.eq.'c' .or. 
!$$$     &           quantity.eq.'H' .or. quantity.eq.'h') then
!$$$               obs = act/float(lsize**3*tsize) !energy per site
!$$$            elseif(quantity.eq.'E' .or. quantity.eq.'e') then
!$$$               obs = act/float(lsize**3) !internal energy
!$$$            elseif(quantity.eq.'B' .or. quantity.eq.'b') then
!$$$               obs = (act/float(lsize**3))**2
!$$$            elseif(quantity.eq.'R' .or. quantity.eq.'r') then
!$$$               obs = act/float(lsize**3)
!$$$            endif
!$$$!---------------------------------------------------
!$$$!     this is for the large N_c Wilson action, specific heat:
!$$$            read(55,*,end=435) beta_read,act,aux,aux,obs
!$$$            obs = act/float(lsize**3*tsize)
!---------------------------------------------------
!$$$!     this is for the large N_c Wilson action, Binder cumulant:
!$$$            read(55,*,end=435) beta_read,act,aux,aux,obs
!$$$            obs = obs**2
!---------------------------------------------------
!$$$!     this is for the large N_c Wilson action for the top. charge:
!$$$            read(55,*,end=435) beta_read,act,obs,aux,aux
!$$$            obs = abs(obs)
!---------------------------------------------------
!!$!     this is for the canonical data from Philippe, prepared with &
!!$!     combine_data.f90:
!!$            read(55,*,end=435) beta_read,act,obs,aux
!!$!            act=(1.-act)*lsize**3*tsize*6 
!!$            if (Quantity.eq.'C' .or. quantity.eq.'c' .or. &
!!$                quantity.eq.'H' .or. quantity.eq.'h') then
!!$               obs = act/float(lsize**3*tsize) !energy per site
!!$            elseif(quantity.eq.'E' .or. quantity.eq.'e') then
!!$               obs = act/float(lsize**3) !internal energy
!!$            elseif(quantity.eq.'B' .or. quantity.eq.'b') then
!!$               obs = (act/float(lsize**3))**2
!!$            elseif(quantity.eq.'R' .or. quantity.eq.'r') then
!!$               obs = act/float(lsize**3)
!!$            elseif(quantity.eq.'P' .or. quantity.eq.'p') then
!!$               obs = aux
!!$            endif
!---------------------------------------------------
!     this is for Kieran's finite-T FP data:
!            read(55,*,end=435) beta_read,act,aux,aux,obs
            ! Version for data after 15 June 2026:
            read(55,*,end=435) beta_read,act,obs,aux,aux
            ! NOTE: the FP action value from the simulations has \beta/3.0 factored out, so
            ! we need to factor in 1/3.0:
            act = act/3.0_dp
!------------------------------------------------------
!            act = (1.0-act/(lsize**3*tsize*6.0))*lsize**3*tsize*6.0
            !            print *,beta_read,act,aux,aux,obs
!            read(55,*,end=435) beta_read,act,obs,aux,aux
            if (Quantity.eq.'C' .or. quantity.eq.'c' .or. &
                quantity.eq.'H' .or. quantity.eq.'h') then
               obs = act/float(lsize**3*tsize) !energy per site
            elseif(quantity.eq.'E' .or. quantity.eq.'e') then
               obs = act/float(lsize**3) !internal energy
            elseif(quantity.eq.'B' .or. quantity.eq.'b') then
               obs = (act/float(lsize**3))**2
            elseif(quantity.eq.'R' .or. quantity.eq.'r') then
               obs = act/float(lsize**3)
            elseif(quantity.eq.'P' .or. quantity.eq.'p') then
               !obs = obs
            endif
!!$!---------------------------------------------------


            if(abs(beta_read-beta(k)).lt.1.e-15) then
               nmeas(k)=nmeas(k)+1
               act_av(0)=act_av(0)+act
               act_av(k)=act_av(k)+act
               max_act(k)=max(max_act(k),act)
               min_act(k)=min(min_act(k),act)
               max_obs=max(max_obs,obs)
               min_obs=min(min_obs,obs)
               if(quantity.eq.'R' .or. quantity.eq.'r') then
                  if (obs.lt.e0) then
                     av_obs(k,0)=av_obs(k,0)+1.
                  else
                     av_obs2(k,0)=av_obs2(k,0)+1.
                  endif
                  if(act .lt. e0) then
                     ir0 = ir0 + 1
                  else
                     ir1 = ir1 + 1
                  endif
               else
                  av_obs(k,0)=av_obs(k,0)+obs
                  av_obs2(k,0)=av_obs2(k,0)+obs**2
               endif
!     the measured quantities:
               energy(nmeas(k),k)=act
               obsrvb(nmeas(k),k)=obs
            endif
         enddo
         print *,'End of file reached, number of lines in data file'
         print *,'is larger than',nbeta*max_meas
!     stop
 435  continue
      close(55)
         

      max_act(0)=max(max_act(0),max_act(k))
      min_act(0)=min(min_act(0),min_act(k))

      print *
      print '("Number of beta=",f10.5," measurements: ",i6)',beta(k),nmeas(k)
      nmeas(0)=nmeas(0)+nmeas(k)
      
      print '("Min. action value:",f20.8)',min_act(k)
      print '("Max. action value:",f20.8)',max_act(k)
      act_av(k)=act_av(k)/float(nmeas(k))
      print '("Average action value:",f15.6)',act_av(k)
      av_obs(k,0)=av_obs(k,0)/float(nmeas(k))
      av_obs2(k,0)=av_obs2(k,0)/float(nmeas(k))
      print '("Average of the observable:  ",f15.6)',av_obs(k,0)
      print '("Average of the observable^2:",f15.6)',av_obs2(k,0)
      if(quantity.eq.'S' .or. quantity.eq.'s') then
         print '("Susceptibility:",f15.6)',lsize**3*(av_obs2(k,0)-av_obs(k,0)**2)
      elseif(quantity.eq.'C' .or. quantity.eq.'c' .or. quantity.eq.'H' .or. quantity.eq.'h' ) then
         print '("Specific heat:",f15.6)',beta(k)**2*lsize**3*tsize*(av_obs2(k,0)-av_obs(k,0)**2)
      elseif(quantity.eq.'B' .or. quantity.eq.'b') then
         print '("Binder cumulant:",f15.6)',lsize**3*(av_obs2(k,0)/(av_obs(k,0)**2)-1.)
      elseif(quantity.eq.'R' .or. quantity.eq.'r') then
!     check:
         if (ir0+ir1 .ne. nmeas(k)) then
            print *,"Warning: ir0+ir1 .ne. nmeas(k):"
            print *,"ir0, ir1, nmeas(k)=",ir0, ir1, nmeas(k)
         endif
!         print '("Ratio of weights,",f15.6,", e0=",f15.6)',
!     &       log(float(ir0)/float(ir1)) , e0
         print '("Ratio of weights,",f15.6,", e0=",f15.6)',log(av_obs(k,0)/av_obs2(k,0)) , e0
         
      endif
      print *,'*******************************************'
      print *

!     end of loop over beta values:
      enddo
      if(quantity.eq.'P' .or. quantity.eq.'p') then
         quantity='S'
      endif

      act_av(0)=act_av(0)/float(nmeas(0))
      print '("Action average of all measurements:",f15.6)',act_av(0)
      print *
      print  '("Min. obs. value:",f15.6)',min_obs
      print '("Max. obs. value:",f15.6)',max_obs
      print *

!     calculate the number of bins:
      nbin=int((max_act(0)-min_act(0))/bin_width)+1
!     calculate the bin width for the observable:
      obsbin_width=(max_obs-min_obs)/float(obs_nbin)
      if(nbin.gt.max_nbin) then
         print *,'Number of bins is too big, increase bin width!'
         print *,'or increase max_nbin!'
         print *,'nbin, max_nbin=',nbin, max_nbin
         stop
      endif
      print '("Number of energy bins:",i5)',nbin
!     generate table:    bin number <-> energy
      act=min_act(0)-bin_width/2.
      do i=0,nbin
         act=act+bin_width
!         bin_energy(i)=act
         bin_energy(i)=act-act_av(0)
      enddo
!     generate table:    bin number <-> obs. value
      obs=min_obs-obsbin_width/2.
      do i=0,obs_nbin
         obs=obs+obsbin_width
         obsbin_val(i)=obs
      enddo

!     binning the data: 
!--------------------------------------------------------------------
!     loop over the beta values:
      do k=1,nbeta
      print *
      print *,'Binning the data:'
      print *,'from beta=',beta(k)
!     set some arrays to zero:
      do i=0,nbin
         en_prob(i,k,0)=0.
         eff_obs(i,k,0)=0.
         eff_obs2(i,k,0)=0.
         do j=0,obs_nbin
            obs_en_prob(j,i,k)=0.
         enddo
      enddo
      do i=0,obs_nbin
         obs_prob(i,k)=0.
      enddo

      do i=1,nmeas(k)
         act=energy(i,k)
         obs=obsrvb(i,k)
!     calculate the bin number:
         bin_number=int((act-min_act(0))/bin_width)
         en_prob(bin_number,k,0)=en_prob(bin_number,k,0)+1.
!     calculate the obs. bin number:
         obsbin_number=int((obs-min_obs)/obsbin_width)
         obs_prob(obsbin_number,k)=obs_prob(obsbin_number,k)+1./nmeas(k)
!     fill in the energy-obs. probability:
         obs_en_prob(obsbin_number,bin_number,k)= &
             obs_en_prob(obsbin_number,bin_number,k)+1./nmeas(k)
!         obs_prob(obsbin_number,k)=
!     &        obs_prob(obsbin_number,k)+exp(-beta(k)*(act-act_av(k)))
!     calculate the effective observable:
         eff_obs(bin_number,k,0)=eff_obs(bin_number,k,0)+obs
         eff_obs2(bin_number,k,0)=eff_obs2(bin_number,k,0)+obs**2
      enddo

!     normalize the observable:
      do i=0,nbin
         if(en_prob(i,k,0).gt.0.) then
            eff_obs(i,k,0)=eff_obs(i,k,0)/en_prob(i,k,0)
            eff_obs2(i,k,0)=eff_obs2(i,k,0)/en_prob(i,k,0)
         endif
      enddo

!     normalize the obs. prob. distr.:
      aux=0.
      do i=0,obs_nbin
         aux=aux+obs_prob(i,k)
      enddo
      do i=0,obs_nbin
         obs_prob(i,k)=obs_prob(i,k)/aux
      enddo
!     end of loop over the beta values:
      enddo

!     write out the normalized distributions:
!      print *,'Writing out the distribution:'
      open(55,file='en_distr.plo',form='formatted',status='unknown')
      write(55,'("action ",20f15.4)') (beta(k),k=1,nbeta)
      do i=0,nbin
         act=bin_energy(i)+act_av(0)
!     divide by the bin width for comparison:         
!         write(55,*) act, en_prob(i)/bin_width
         write(55,'(20f15.5)') act/(lsize**3*tsize),(en_prob(i,k,0)/float(nmeas(k)),k=1,nbeta)
      enddo
      close(55)

! -------------------------------
! initialize random generator
      print *
      open(20,file='iseed.dat',status='old')
      read(20,*) iseed
      close(20)
      iseed=iseed+2
      open(20,file='iseed.dat',status='old')
      write(20,*) iseed
      close(20)
      print '("Initialize random generator:",i9)',iseed
!     ranlux generator:
      level=1
      call rlxd_init(level,iseed)
      ! -------------------------------
!!$      do k=1,nbeta
!!$         call ranlxd(r)
!!$         print *,'ranf()',1.0-ranf(), r(1)
!!$      enddo

!     write out the normalized distributions:
!      print *,'Writing out the distribution:'
      open(55,file='en_btrp_distr.plo',form='formatted',status='unknown')

!     start the bootstrap:
        print *
        print '("Generating",i5," bootstrap samples:")',nbtrp
        print '("Blocksize for bootstrap:",i6)',blocksize
        print *
      do nb=1,nbtrp
!     loop over the beta values:
         do k=1,nbeta
!         print '("Sample",i5)',k
!         print *,'Binning the data:'
         do i=0,nbin
            en_prob(i,k,nb)=0.
            eff_obs(i,k,nb)=0.
            eff_obs2(i,k,nb)=0.
         enddo
!     generate a bootstrap sample by picking one of the blocks:
!     the number of blocks:
         nblock=nmeas(k)/blocksize
!         print *,'DEBUG: k,nblock=', k, nblock
         av_obs(k,nb)=0.
         av_obs2(k,nb)=0.
         do j=1,nblock
            !            rvar=ranf()
            call ranlxd(r)
            rvar = 1.0 - r(1)
            irvar=int(rvar*nblock)+1
!            print *,'DEBUG: rvar,irvar=',rvar,irvar
            do i=1,blocksize
!     and do the binning: 
               act=energy((irvar-1)*blocksize+i,k)
               obs=obsrvb((irvar-1)*blocksize+i,k)
!     calculate the bin number:
               bin_number=int((act-min_act(0))/bin_width)
               en_prob(bin_number,k,nb)=en_prob(bin_number,k,nb)+1.
!     calculate the effective observable:
               eff_obs(bin_number,k,nb)=eff_obs(bin_number,k,nb)+obs
               eff_obs2(bin_number,k,nb)=eff_obs2(bin_number,k,nb)+obs**2 
               if(quantity.eq.'R' .or. quantity.eq.'r') then
                  if (obs.lt.e0) then
                     av_obs(k,nb)=av_obs(k,nb)+1.
                  else
                     av_obs2(k,nb)=av_obs2(k,nb)+1.
                  endif
               else
                  av_obs(k,nb)=av_obs(k,nb)+obs
                  av_obs2(k,nb)=av_obs2(k,nb)+obs**2
               endif
            enddo
         enddo
!     normalize the distribution and the observable:
!         print *,'Normalizing the distribution:'
         do i=0,nbin
            if(en_prob(i,k,nb).gt.0.) then
               eff_obs(i,k,nb)=eff_obs(i,k,nb)/en_prob(i,k,nb)
               eff_obs2(i,k,nb)=eff_obs2(i,k,nb)/en_prob(i,k,nb)
            endif
         enddo
!     normalize the averages:
         av_obs(k,nb)=av_obs(k,nb)/float(blocksize*nblock)
         av_obs2(k,nb)=av_obs2(k,nb)/float(blocksize*nblock)

!     end of loop over beta values:
      enddo

!!$      do i=0,nbin
!!$         act=bin_energy(i)+act_av(0)
!!$         !     divide by the bin width for comparison:         
!!$!!!$         write(55,*) act, en_prob(i)/bin_width
!!$         write(55,'(4f15.8)') act, &
!!$              (en_prob(i,k,nb)/float(nblock*blocksize),k=1,nbeta)
!!$      enddo
      
!     end of bootstrap loop:
      enddo

      close(55)

!     calculate the average and the bootstrap error for each beta
!     value seperately (and write it out):
!----------------------------------------------------------------
      open(56,file='obs+err.plo',form='formatted',status='unknown')
      if(quantity.eq.'L' .or. quantity.eq.'l') then
      do k=1,nbeta
!     calculate the bootstrap average:
         aux=0.
         do nb=1,nbtrp
            aux=aux+av_obs(k,nb)
         enddo
         aux=aux/float(nbtrp)
!     calculate the error:
         norm=0.
         do nb=1,nbtrp
            norm=norm+(av_obs(k,nb)-aux)**2
         enddo
         norm=norm/float(nbtrp-1)
         norm=sqrt(norm)
!     and write it out:
         write(56,*) beta(k),av_obs(k,0),norm
      enddo
      endif
      if(quantity.eq.'S' .or. quantity.eq.'s') then

      do k=1,nbeta
!     calculate the bootstrap average:
         aux=0.
         do nb=1,nbtrp
            aux=aux+lsize**3*(av_obs2(k,nb)-av_obs(k,nb)**2)
!            print *,'nb, av_obs(.,nb), av_obs2(.,nb)=', nb, av_obs(k,nb), av_obs2(k,nb)
         enddo
         aux=aux/float(nbtrp)
!     calculate the error:
         norm=0.
         do nb=1,nbtrp
            norm=norm+(lsize**3*(av_obs2(k,nb)-av_obs(k,nb)**2)-aux)**2
         enddo
         norm=norm/float(nbtrp-1)
         norm=sqrt(norm)
!     and write it out:
         write(56,*) beta(k),lsize**3*(av_obs2(k,0)-av_obs(k,0)**2),norm
      enddo
      endif
      if(quantity.eq.'C' .or. quantity.eq.'c' .or. quantity.eq.'H' .or. quantity.eq.'h') then
      do k=1,nbeta
!     calculate the bootstrap average:
         aux=0.
         do nb=1,nbtrp
            aux=aux+beta(k)**2*lsize**3*tsize*(av_obs2(k,nb)-av_obs(k,nb)**2)
         enddo
         aux=aux/float(nbtrp)
!     calculate the error:
         norm=0.
         do nb=1,nbtrp
            norm=norm+ &
                (beta(k)**2*lsize**3*tsize*(av_obs2(k,nb)-av_obs(k,nb)**2)-aux)**2
         enddo
         norm=norm/float(nbtrp-1)
         norm=sqrt(norm)
!     and write it out:
         write(56,*) beta(k),beta(k)**2*lsize**3*tsize*(av_obs2(k,0)-av_obs(k,0)**2),norm
      enddo
      endif
      if(quantity.eq.'B' .or. quantity.eq.'b') then
      do k=1,nbeta
!     calculate the bootstrap average:
         aux=0.
         do nb=1,nbtrp
            aux=aux + lsize**3*(av_obs2(k,nb)/(av_obs(k,nb)**2) - 1.)
         enddo
         aux=aux/float(nbtrp)
!     calculate the error:
         norm=0.
         do nb=1,nbtrp
            norm=norm+ &
                (lsize**3*(av_obs2(k,nb)/(av_obs(k,nb)**2) - 1.)-aux)**2
         enddo
         norm=norm/float(nbtrp-1)
         norm=sqrt(norm)
!     and write it out:
         write(56,*) beta(k),lsize**3*(av_obs2(k,0)/(av_obs(k,0)**2) - 1.),norm
      enddo
      endif
      if(quantity.eq.'E' .or. quantity.eq.'e') then
      do k=1,nbeta
!     calculate the bootstrap average:
         aux=0.
         do nb=1,nbtrp
            aux=aux+av_obs(k,nb)
         enddo
         aux=aux/float(nbtrp)
!     calculate the error:
         norm=0.
         do nb=1,nbtrp
            norm=norm+(av_obs(k,nb)-aux)**2
         enddo
         norm=norm/float(nbtrp-1)
         norm=sqrt(norm)
!     and write it out:
         write(56,*) beta(k),av_obs(k,0),norm
      enddo
      endif
      if(quantity.eq.'R' .or. quantity.eq.'r') then
      do k=1,nbeta
!     calculate the bootstrap average:
         aux=0.
         do nb=1,nbtrp
            aux=aux+log(av_obs(k,nb)/av_obs2(k,nb))
         enddo
         aux=aux/float(nbtrp)
!     calculate the error:
         norm=0.
         do nb=1,nbtrp
            norm=norm+(log(av_obs(k,nb)/av_obs2(k,nb))-aux)**2
         enddo
         norm=norm/float(nbtrp-1)
         norm=sqrt(norm)
!     and write it out:
         write(56,*) beta(k),log(av_obs(k,0)/av_obs2(k,0)),norm
      enddo
      endif



      close(56)


!     first get the autocorrelation times:
!     for the time being set them all equal: 
!     they should maybe also be included in the bootstrap analysis!
      call get_autocor(nbeta,nmeas,energy,obsrvb,tau)
      do k=1,nbeta
         gtau(k)=1.+2*tau(k)
!         gtau(k)=1.
      enddo

!$$$!     calculate the averaged integratd autocorrelation time 
!$$$!     (used for calculating the error):
!$$$      tau_avg=0.
!$$$      norm=0.
!$$$      do k=1,nbeta
!$$$         pk=(1.+2.*tau(k))/float(nmeas(k))
!$$$         norm=norm+pk
!$$$         tau_avg=tau_avg+tau(k)*pk
!$$$      enddo
!$$$      tau_avg=tau_avg/norm
!$$$      print *
!$$$      print '("Avg. int. autocorrelation time:",f12.5)',tau_avg

!     calculate the max. integrated autocorrelation time 
!     (used for calculating the error):
      tau_avg=0.
      do k=1,nbeta
         tau_avg=max(tau_avg,tau(k))
      enddo
      print *
      print '("Max. int. autocorrelation time:",f12.5)',tau_avg

!     calculate the beta differences (dbeta=beta-beta_av):
      do k=1,nbeta
         dbeta(k)=beta(k)-beta(0)
      enddo

!     and set the initial free energies:
      do nb=0,nbtrp
         do k=1,nbeta
            fe(k,nb)=dbeta(k)*(act_av(k)-act_av(0))
!            if(nb==0) print *,'k,dbeta(k),fe(k)=',k,dbeta(k),fe(k,0)
         enddo
      enddo

!     solve now for the free energies 'fe' (treated as free parameters)
!----------------------------------------------------------------------
!     in order to get the unique spectral density function 
!     Wbar(S) ( = en_dens(0:nbin,0:nbtrp)):
      print *,'Solve now for the free energies'
      call solve_fe(nbin,nbeta,nmeas,bin_energy,fe,dbeta,en_prob,gtau,en_dens)


!     calculate the effective observable for every bootstrap sample:
!-------------------------------------------------------------------
      open(54,file='obs_btrp_distr.plo',form='formatted',status='unknown')

!     loop over the bootstrap samples:
      do nb=0,nbtrp
!     loop over the energies:
         do i=0,nbin
            aux=0.
            do k=1,nbeta
               aux=aux+exp(-dbeta(k)*bin_energy(i)+fe(k,nb))*nmeas(k)/gtau(k)
            enddo
            obs_bar(i,nb)=0.
            obs2_bar(i,nb)=0.
            norm=0.
            do k=1,nbeta
               pk=exp(-dbeta(k)*bin_energy(i)+fe(k,nb))*nmeas(k)/gtau(k)
               pk=pk/aux
               norm=norm+pk
               obs_bar(i,nb)=obs_bar(i,nb)+pk*eff_obs(i,k,nb)
               obs2_bar(i,nb)=obs2_bar(i,nb)+pk*eff_obs2(i,k,nb)
            enddo
!            print *,'i,Norm of p_k(S=i)=',i,norm
         enddo

!     end loop over bootstrap samples:
      enddo

!     calculate the energy-observable distribution:
      nb=0
      do i=0,nbin
         aux=0.
         do k=1,nbeta
            aux=aux+exp(-dbeta(k)*bin_energy(i)+fe(k,nb))*nmeas(k)/gtau(k)
         enddo
         do j=0,obs_nbin
            obs_en_prob(j,i,0)=0.
            do k=1,nbeta
               pk=exp(-dbeta(k)*bin_energy(i)+fe(k,nb))*nmeas(k)/gtau(k)
               pk=pk/aux
               obs_en_prob(j,i,0)=obs_en_prob(j,i,0)+pk*obs_en_prob(j,i,k)
            enddo
         enddo
      enddo
!     print out the reweighted observable distribution
!     w(P)=sum_S w(P,S) exp(-dbeta*S):
!     write out the normalized obs. distributions:
      norm=0.
      do i=0,obs_nbin
         aux=0.
         do j=0,nbin
!     Note: dbeta(0)=0 -> exp(-dbeta*S)=0.:
            aux=aux+obs_en_prob(i,j,0)*en_dens(j,0)
         enddo
         norm=norm+aux
      enddo
      open(56,file='obs_distr.plo',form='formatted',status='unknown')
      do i=0,obs_nbin
         aux=0.
         do j=0,nbin
!     Note: dbeta(0)=0 -> exp(-dbeta*S)=0.:
            aux=aux+obs_en_prob(i,j,0)*en_dens(j,0)/norm
         enddo
         obs=obsbin_val(i)
         write(56,'(10f10.5)') obs,(obs_prob(i,k),k=1,nbeta),aux
      enddo
      close(56)




!!$      call print_obs_distr(obs_en_prob,en_dens,nbin,obsbin_val,lsize,tsize)


!$$$      do i=0,nbin
!$$$         act=bin_energy(i)+act_av(0)
!$$$!         write(54,'(4f15.8)') act,obs_bar(i,0),obs2_bar(i,0)
!$$$         write(54,'(4f15.8)') act,obs_bar(i,0)
!$$$      enddo

!$$$      close(54)

      return
      end
      
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!--------------------------------------------------------------------------
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine check_obs_distr_maxima(obs_en_prob,en_dens,nbin,obsbin_val,lsize,tsize,beta_pr,pmax1,pmax2)

        
        ! This routine calculates the probability distribution of the observable and determines the maximum

        use fs_multi_par
        
        implicit none
        
        !     the spectral densities:
        real(kind=dp), intent(in) :: en_dens(0:max_nbin,0:nbtrp)
        
        !     for the observable probability distribution:
        real(kind=dp), intent(in) :: obs_en_prob(0:obs_nbin,0:max_nbin,0:max_nbeta)
        
        !     the bin values of the observable:
        real(kind=dp), intent(in) :: obsbin_val(0:obs_nbin)

        !     number of energy bins:
        integer, intent(in) :: nbin

        !     lattice size in spatial and temporal direction:
        integer, intent(in) :: lsize,tsize

        ! The current beta value:
        real(kind=dp), intent(in) :: beta_pr

        ! The two maxima:
        real(kind=dp), intent(out) :: pmax1, pmax2

        !     auxiliary for calculating bootstrap error on obs. distribution:
        real(kind=dp) :: paux(0:obs_nbin,0:nbtrp)

        !     auxiliaries for finding the extremas:
        integer imin(0:nbtrp), imax1(0:nbtrp),imax2(0:nbtrp)
        real(kind=dp) :: pmin

        !     auxiliaries:
        integer nb,i,j
        real(kind=dp) :: norm,aux,err,fact,add,aux0

        !     loop over observable bins:
        nb=0
        norm=0.
        do i=0,obs_nbin
           aux=0.
           !     loop over energies:
           do j=0,nbin
              !               aux=aux+obs_en_prob(i,j,0)*en_dens(j,nb)
              aux=aux+obs_en_prob(i,j,0)*new_en_prob(j,nb)
           enddo
           paux(i,nb)=aux
           norm=norm+aux
        enddo
        !     norm the distributions:
        do i=0,obs_nbin
           paux(i,nb)=paux(i,nb)/norm
        enddo

        !     now search for the two peaks and the minimum in order to
        !      determine ln(p_min/p_max):
!-------------------------------------------------------------
        !     first do it for extrema position fixed from nb=0:
        nb=0
        pmax1=0.
        pmax2=0.
        pmin=1.
        !      do i=1,15
        do i=1,obs_nbin/3
           if(paux(i,nb) .ge. pmax1) then
              pmax1=paux(i,nb)
              imax1(nb)=i
           endif
        enddo
        do i=obs_nbin/2,obs_nbin
           if(paux(i,nb) .ge. pmax2) then
              pmax2=paux(i,nb)
              imax2(nb)=i
           endif
        enddo
        !      do i=10,35
        !      do i=15,25
        do i=imax1(nb),imax2(nb)
           if(paux(i,nb) .le. pmin) then
              pmin=paux(i,nb)
              imin(nb)=i
           endif
        enddo

        print '("i_max1, p_max1, i_max2, p_max2=",f12.7,i4,f16.11,i4,f16.11,i4,f16.11)', beta_pr,imax1(0), pmax1, imax2(0), pmax2, imin(0), pmin
        
        
      end subroutine check_obs_distr_maxima

      
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!--------------------------------------------------------------------------
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine print_obs_distr(obs_en_prob,en_dens,nbin,obsbin_val,lsize,tsize)


!     This routine calculates and prints out the probability distribution
!     of the observable together with the error.
!     Note:
!     w(P)=sum_S w(P,S) exp(-dbeta*S)
!     dbeta(0)=0 -> exp(-dbeta*S)=1.


      use fs_multi_par

      implicit none

!     the spectral densities:
      real(kind=dp), intent(in) :: en_dens(0:max_nbin,0:nbtrp)

!     for the observable probability distribution:
      real(kind=dp), intent(in) :: obs_en_prob(0:obs_nbin,0:max_nbin,0:max_nbeta)

!     the bin values of the observable:
      real(kind=dp), intent(in) :: obsbin_val(0:obs_nbin)

!     number of energy bins:
      integer, intent(in) :: nbin

!     lattice size in spatial and temporal direction:
      integer, intent(in) :: lsize,tsize

!     auxiliary for calculating bootstrap error on obs. distribution:
      real(kind=dp) :: paux(0:obs_nbin,0:nbtrp)

!     auxiliaries for finding the extremas:
      integer imin(0:nbtrp), imax1(0:nbtrp),imax2(0:nbtrp)
      real(kind=dp) :: pmin,pmax1,pmax2

!     auxiliaries:
      integer nb,i,j
      real(kind=dp) :: norm,aux,err,fact,add,aux0

!     loop over bootstrap samples:
      do nb=0,nbtrp
         norm=0.
!     loop over observable bins:
         do i=0,obs_nbin
            aux=0.
!     loop over energies:
            do j=0,nbin
!               aux=aux+obs_en_prob(i,j,0)*en_dens(j,nb)
               aux=aux+obs_en_prob(i,j,0)*new_en_prob(j,nb)
            enddo
            paux(i,nb)=aux
            norm=norm+aux
         enddo
!     norm the distributions:
         do i=0,obs_nbin
            paux(i,nb)=paux(i,nb)/norm
         enddo
      enddo

!!$      ! Comparison of the density of states and the new energy probablility distribution:
!!$      nb=0
!!$      print *,'nb=0: en_dens(j,nb),new_en_prob(j,nb)'
!!$      do j=0,nbin
!!$         print *,en_dens(j,nb),new_en_prob(j,nb)
!!$      enddo

      
!     open file for writing results: 
      open(54,file='obs_btrp_distr.plo',form='formatted', status='unknown')

!     calculate the bootstrap error:
!-----------------------------------
!     loop over observable bins:
      do i=0,obs_nbin
!     the bootstrap average:
         aux=0.
         do nb=1,nbtrp
            aux=aux+paux(i,nb)
         enddo
         aux=aux/float(nbtrp)
!     the bootstrap error:
         err=0.
         do nb=1,nbtrp
            err=err + (aux - paux(i,nb))**2
         enddo
         err=err/float(nbtrp-1)
         err=sqrt(err)
!     print it out:
         write(54,'(3f15.10)') obsbin_val(i),paux(i,0),err
      enddo

      close(54)

      
!     now search for the two peaks and the minimum in order to
!      determine ln(p_min/p_max):
!-------------------------------------------------------------
!     first do it for extrema position fixed from nb=0:
      nb=0
      pmax1=0.
      pmax2=0.
      pmin=1.
      !      do i=1,15
      do i=1,obs_nbin/3
         if(paux(i,nb) .ge. pmax1) then
            pmax1=paux(i,nb)
            imax1(nb)=i
         endif
      enddo
      do i=obs_nbin/2,obs_nbin
         if(paux(i,nb) .ge. pmax2) then
            pmax2=paux(i,nb)
            imax2(nb)=i
         endif
      enddo
!      do i=10,35
!      do i=15,25
      do i=imax1(nb),imax2(nb)
      if(paux(i,nb) .le. pmin) then
            pmin=paux(i,nb)
            imin(nb)=i
         endif
      enddo


!     print out the p_max1,2 which should be equal:
!-------------------------------------------------- 
      print *
      print '(" The surface tension:")'
      print '(" ************************************")'
      i=imax1(0)
      aux=0.
      do nb=1,nbtrp
         aux=aux+paux(i,nb)
      enddo
      aux=aux/float(nbtrp)
      err=0.
      do nb=1,nbtrp
         err=err + (aux - paux(i,nb))**2
      enddo
      err=err/float(nbtrp-1)
      err=sqrt(err)
      print '(" p_max1=",f12.7," +/- ",f12.7)',paux(i,0),err
      i=imax2(0)
      aux=0.
      do nb=1,nbtrp
         aux=aux+paux(i,nb)
      enddo
      aux=aux/float(nbtrp)
      err=0.
      do nb=1,nbtrp
         err=err + (aux - paux(i,nb))**2
      enddo
      err=err/float(nbtrp-1)
      err=sqrt(err)
      print '(" p_max2=",f12.7," +/- ",f12.7)',paux(i,0),err
      i=imin(0)
      aux=0.
      do nb=1,nbtrp
         aux=aux+paux(i,nb)
      enddo
      aux=aux/float(nbtrp)
      err=0.
      do nb=1,nbtrp
         err=err + (aux - paux(i,nb))**2
      enddo
      err=err/float(nbtrp-1)
      err=sqrt(err)
      print '(" p_min =",f12.7," +/- ",f12.7)',paux(i,0),err


      ! The correct expression should be (cf. Lucini et al., arXiv:0502003):
      !    \hat \sigma = -1/2 Lt^2/Ls^2 { \ln p_min/p_max - 1/2 \ln Ls - c}
      ! where \hat \sigma = \sigma /T_c^3 is dimensionless.
      ! We define
      !    \hat \Sigma = -1/2 Lt^2/Ls^2 \ln p_min/p_max + 1/2 Lt^2/Ls^2 1/2 \ln Ls
      ! which has the correct thermodynamic limit
      !    \hat \Sigma -(Lt/Ls -> oo)-> \hat \sigma
      ! and FV corrections \propto Lt^2/L_s^2:
      !---------------------------------------------------------------------------
      fact=0.5*(float(tsize)/float(lsize))**2
      add=0.5*fact*log(float(lsize))
!     now calculate the surface tension from these fixed locations:
      open(52,file='surface_tension.plo',form='formatted',status='unknown')
      ! use the first maximum:
      !-----------------------
      j=imin(0)
      i=imax1(0)
      aux=0.
      do nb=1,nbtrp
         aux=aux+log(paux(j,nb)/paux(i,nb))
      enddo
      aux=aux/float(nbtrp)
      err=0.
      do nb=1,nbtrp
         err=err + (aux - log(paux(j,nb)/paux(i,nb)))**2
      enddo
      err=err/float(nbtrp-1)
      err=sqrt(err)
      print '("ibin_min=",i8,", ibin_max1=",i8,", \hat \Sigma = ",f12.7," +\- ",f12.7,", (bias corr. =",f12.7,")")', &
           j,i,-fact*log(paux(j,0)/paux(i,0))+add,fact*err,-fact*aux+add
      write(52,'(2i8,3f12.7)') j,i,-fact*log(paux(j,0)/paux(i,0))+add,fact*err,-fact*aux+add
      ! use the second maximum:
      !------------------------
      j=imin(0)
      i=imax2(0)
      do nb=1,nbtrp
         aux=aux+log(paux(j,nb)/paux(i,nb))
      enddo
      aux=aux/float(nbtrp)
      err=0.
      do nb=1,nbtrp
         err=err + (aux - log(paux(j,nb)/paux(i,nb)))**2
      enddo
      err=err/float(nbtrp-1)
      err=sqrt(err)
      print '("ibin_min=",i8,", ibin_max2=",i8,", \hat \Sigma = ",f12.7," +\- ",f12.7,", (bias corr. =",f12.7,")")', &
           j,i,-fact*log(paux(j,0)/paux(i,0))+add,fact*err,-fact*aux+add
      write(52,'(2i8,3f12.7)') j,i,-fact*log(paux(j,0)/paux(i,0))+add,fact*err,-fact*aux+add

      ! use the average of the two maxima:
      !-----------------------------------
!$$$      j=imin(0)
!$$$      i=imax2(0)
!$$$      do nb=1,nbtrp
!$$$         aux=aux+log(2.*paux(j,nb)/(paux(i,nb)+paux(imax1(0),nb)))
!$$$      enddo
!$$$      aux=aux/float(nbtrp)
!$$$      err=0.
!$$$      do nb=1,nbtrp
!$$$         err=err + (aux - 
!$$$     &        log(2.*paux(j,nb)/(paux(i,nb)+paux(imax1(0),nb))))**2
!$$$      enddo
!$$$      err=err/float(nbtrp-1)
!$$$      err=sqrt(err)
!$$$      norm = log(2.*paux(j,0)/(paux(i,0)+paux(imax1(0),0)))
!$$$      write(52,'(2i4,i8,3f12.7)') imax1(0),i,j,norm,err,aux
      ! use the geometric mean of the two maxima:
      !------------------------------------------
      j=imin(0)
      i=imax2(0)
      nb=0
      aux0 = log(paux(j,nb)/(paux(i,nb)**0.5*paux(imax1(0),nb)**0.5))
      do nb=1,nbtrp
         aux=aux+log(paux(j,nb)/(paux(i,nb)**0.5*paux(imax1(0),nb)**0.5))
      enddo
      aux=aux/float(nbtrp)
      err=0.
      do nb=1,nbtrp
         err=err + (aux - &
             log(paux(j,nb)/(paux(i,nb)**0.5* paux(imax1(0),nb)**0.5)))**2
      enddo
      err=err/float(nbtrp-1)
      err=sqrt(err)
      norm = log(paux(j,0)/(paux(i,0)**0.5*paux(imax1(0),0)**0.5))
      norm = -fact*norm+add
      err = fact*err
      print '("ibin_min=",i8,", ibin_max2=",i8,", \hat \Sigma = ",f12.7," +\- ",f12.7,", (bias corr. =",f12.7,")")', &
           j,i,-fact*aux0+add,err,-fact*aux+add
      
      write(52,'(2i8,3f12.7)') j,i,-fact*aux0+add,fact*err,-fact*aux+add
!      write(52,'(2i4,i8,3f12.7)') imax1(0),i,j,norm,err,-fact*aux+add

      !     Print out the surface tension:
      !-----------------------------------
      print *
      print '(" sigma/T_c^3 = ",3f12.7)',(float(tsize)/float(lsize))**2,norm,err
      print *



!     now calculate it for each bootstrap sample:
      do nb=0,nbtrp
         pmax1=0.
         pmax2=0.
         pmin=1.
         do i=1,obs_nbin/3
            if(paux(i,nb) .ge. pmax1) then
               pmax1=paux(i,nb)
               imax1(nb)=i
            endif
         enddo
         do i=obs_nbin/2,obs_nbin
            if(paux(i,nb) .ge. pmax2) then
               pmax2=paux(i,nb)
               imax2(nb)=i
            endif
         enddo
!         do i=10,35
         do i=imax1(nb),imax2(nb)
            if(paux(i,nb) .le. pmin) then
               pmin=paux(i,nb)
               imin(nb)=i
            endif
         enddo
      enddo

      aux=0.
      do nb=0,nbtrp
!     the bootstrap average:
         aux=aux+log(paux(imin(nb),nb)/paux(imax1(nb),nb))
      enddo
      aux=aux/float(nbtrp)
      err=0.
      do nb=1,nbtrp
!     the bootstrap error:         
         err=err + (aux - log(paux(imin(nb),nb)/paux(imax1(nb),nb)))**2
      enddo
      err=err/float(nbtrp-1)
      err=sqrt(err)
      write(52,'(2i8,3f12.7)') imin(0),imax1(0),&
          -fact*log(paux(imin(0),0)/paux(imax1(0),0))+add,fact*err,-fact*aux+add


      aux=0.
      do nb=0,nbtrp
!     the bootstrap average:
         aux=aux+log(paux(imin(nb),nb)/paux(imax2(nb),nb))
      enddo
      aux=aux/float(nbtrp)
      err=0.
      do nb=1,nbtrp
!     the bootstrap error:         
         err=err + (aux - log(paux(imin(nb),nb)/paux(imax2(nb),nb)))**2
      enddo
      err=err/float(nbtrp-1)
      err=sqrt(err)
      write(52,'(2i8,3f12.7)') imin(0),imax2(0),&
          -fact*log(paux(imin(0),0)/paux(imax2(0),0))+add,fact*err,-fact*aux+add


      close(52)

    end subroutine print_obs_distr


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!--------------------------------------------------------------------------
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine solve_fe(nbin,nbeta,nmeas,bin_energy,fe0,dbeta,en_prob,gtau,en_dens)

!     This routine solves for the free energies self-consistently by using 
!     the derivatives of the new values of f_n as functions of the old 
!     values in a iteration process:
!
!     Define
!     f_k = -Log(Zbar(beta_k,f_i)) -> F_k = f_k + Log(Zbar(beta_k,f_i))
!
!     and calculate the new value at iteration (i+1) through:
!     f_k^{(i+1)} = f_k^{(i)} - M_{ki}^{-1} F_i,
!
!     where M_{ki} = \partial F_k / \partial f_i.
!
!     The selfconsistently determined spectral density function 'wbar'
!     is returned.

      use fs_multi_par

      implicit none

!     number of bins:
      integer :: nbin
!     number of beta values:
      integer :: nbeta
!     number of measurements:
      integer :: nmeas(0:nbeta)

!     energy of the bins:
      real(kind=dp) :: bin_energy(0:nbin)

!     the beta differences:
      real(kind=dp) :: dbeta(1:nbeta)

!     the distributions (histograms):
      real(kind=dp) :: en_prob(0:max_nbin,1:max_nbeta,0:nbtrp)

!     the error of the probability distributions:
      real(kind=dp) :: gtau(1:nbeta)

!     the spectral densities:
      real(kind=dp) :: en_dens(0:max_nbin,0:nbtrp)

!     the free energies:
      real(kind=dp) :: fe0(1:nbeta,0:nbtrp),fe(1:nbeta)

!     the partition functions:
      real(kind=dp) :: zbar(1:nbeta)

!     F_k = f_k + Log(Zbar(beta_k,f_i))
      real(kind=dp) :: ff(1:nbeta)

!     the curvature matrix:
      real(kind=dp) :: mki(1:nbeta,1:nbeta)

!     mf(k) == M_(ki)^{-1} F_i:
      real(kind=dp) :: mf(1:nbeta)

!     some parameters concerning the iterative procedure:
!--------------------------------------------------------
!     max. number of iterations:
      integer iter
!     absolute accuracy:
      real(kind=dp) :: acc
      parameter(iter=1000,acc=1.e-5)

!     some auxiliary parameters:
!-------------------------------
      integer nb,k,i,ii,ifail
      real(kind=dp) :: nom,denom,aux,wkspce(1:nbeta),avg
!     the spectral density:
      real(kind=dp) :: wbar(0:nbin)


!     start loop over the bootstrap samples:
      do nb=0,nbtrp
!         print *,'Bootstrap sample:',nb

         !     set the initial free energies:
         do k=1,nbeta
            fe(k)=fe0(k,nb)
         enddo

!     start the iteration:
      do ii=1,iter
!$$$         print *,'iteration:',ii

!     Calculate the function F_k (=ff(k)):
!-----------------------------------
!     first get the density of states Wbar(S,f_i) for the given {f_i}:
!         print *,'get the density of states Wbar(S,{f_i}):'
!         print *,'fe=',fe
!!$!     loop over the energies:
!!$         do i=0,nbin
!!$!     calculate the numerator 'nom' and the denominator 'denom':
!!$            nom=0.
!!$            denom=0.
!!$            do k=1,nbeta
!!$               nom=nom+en_prob(i,k,nb)/gtau(k)
!!$               aux=exp(-dbeta(k)*bin_energy(i)+fe(k))
!!$               denom=denom+nmeas(k)*aux/gtau(k)
!!$            enddo
!!$!     the spectral density function:
!!$            wbar(i)=nom/denom
!!$!     end of loop over the energies:
!!$         enddo
!!$
!!$!     second get the function F_k:
!!$      print *,'get the function F_k({f_i}):'
!!$      do k=1,nbeta
!!$!     loop over the energies:
!!$         aux=0.
!!$         do i=0,nbin
!!$            aux=aux+wbar(i)*exp(-dbeta(k)*bin_energy(i))
!!$         enddo
!!$         ff(k)=fe(k)+log(aux)
!!$      enddo
!!$
!!$!     third get the curvature M_{ki}:
!!$      print *,'get the curvature M_{ki}:'
!!$      call get_mki(nbin,nbeta,nmeas,bin_energy,dbeta,gtau,fe,wbar,mki)
!!$
!!$!     invert the matrix, i.e. M_(ki)^{-1} F_i:
!!$      print *,'invert the matrix:'
!!$      ifail=-1
!!$      call F04ARF(mki,nbeta,ff,nbeta,mf,wkspce,ifail)
!!$      if(ifail.ne.0) then
!!$         print *,'ifail=',ifail
!!$         print *,'ifail=1 -> The matrix is singular!'
!!$      endif
!!$
!!$!     and calculate the new free energy values:
!!$      print *,'get the new free energy values:'
!!$      aux=0.
!!$      do k=1,nbeta
!!$         fe(k)=fe(k)-mf(k)
!!$         aux=max(aux,abs(mf(k)))
!!$      enddo
!!$
!!$      print '("Iteration",i6,": max. correction=",f15.9)',ii,aux
!!$      print '("fe:",20f15.6)',fe
!!$      if(aux.lt.acc) then
!!$         print *,'desired accuracy reached!'
!!$         exit
!!$      endif
!!$
!!$!     end of iteration loop:
!!$      enddo
!!$      if(ii>=iter) then
!!$         print '("No solution found after",i6," iterations!")',iter
!!$      endif
 
!     checking the solution:
!------------------------------------------------------
!      print *,'Checking the solution:'
!      print *,'**********************'
!      print '("fe:",20f15.6)',fe
!      print *,'Calculating the partition functions:'
!     first get the density of states Wbar(S,f_i) for the given {f_i}:
!     loop over the energies:

         do i=0,nbin
!$$$            print *,'bin=',i
!$$$            print *,'bin_energy(i)=',bin_energy(i)
!     calculate the numerator 'nom' and the denominator 'denom':
            nom=0.
            denom=0.
            do k=1,nbeta
               nom=nom+en_prob(i,k,nb)/gtau(k)
               aux=exp(-dbeta(k)*bin_energy(i)+fe(k))
               denom=denom+nmeas(k)*aux/gtau(k)
!$$$               print *,'exp(-dbeta(k)*bin_energy(i)+fe(k))=',aux
            enddo
!     the spectral density function:
            wbar(i)=nom/denom
!     end of loop over the energies:
         enddo

!$$$         print *,'Calling part_funct'

      call part_funct(nbin,nbeta,bin_energy,dbeta,wbar,zbar)
      avg=0.
      do k=1,nbeta
         ff(k)=-log(zbar(k))
         avg=avg+ff(k)
      enddo
      avg=avg/float(nbeta)
!      print '("-log(zbar):",20f15.6)',ff
!--------------------------------------------------------
      aux=0.
      do k=1,nbeta
         ff(k)=ff(k)-avg
         mf(k)=fe(k)-ff(k)
         aux=max(aux,abs(mf(k)))
         fe(k)=ff(k)
      enddo
!!$      print '("Iteration",i6,": max. correction=",5f15.9)',ii,(mf(k),k=1,nbeta)
!!$      print '("fe:",20f15.6)',fe
!!$      if(aux.lt.acc) then
!!$         print *,'desired accuracy reached!'
!!$         goto 679
!!$      endif
!     end of loop over iteration:
      enddo


 679  continue

!     write over the spectral function:
      do i=0,nbin
         en_dens(i,nb)=wbar(i)
      enddo
!     and the free energies:
      do k=1,nbeta
         fe0(k,nb)=fe(k)
      enddo


!     end loop over the bootstrap samples:
      enddo


      return
    end subroutine solve_fe

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!--------------------------------------------------------------------------
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!$
!!$      subroutine get_mki(nbin,nbeta,nmeas,bin_energy,dbeta,gtau,fe,wbar,mki)
!!$
!!$      use fs_multi_par
!!$
!!$      implicit none
!!$
!!$!     number of beta values, energy bins and measurements:
!!$      integer nbeta,nbin,nmeas(1:nbeta)
!!$!     the bin energies and the beta differences:
!!$      real(kind=dp) :: bin_energy(0:nbin),dbeta(1:nbeta)
!!$!     the probability distribution errors:
!!$      real(kind=dp) :: gtau(1:nbeta)
!!$
!!$!     the density of states:
!!$      real(kind=dp) :: wbar(0:nbin)
!!$!     the partition functions:
!!$      real(kind=dp) :: zbar(1:nbeta)
!!$
!!$!     the free energy values:
!!$      real(kind=dp) :: fe(1:nbeta)
!!$
!!$!     the curvature matrix:
!!$      real(kind=dp) :: mki(1:nbeta,1:nbeta)
!!$
!!$!     auxiliary variables:
!!$      integer i,j,k
!!$      real(kind=dp) :: aux,nom,denom,norm(0:nbin)
!!$ 
!!$!     calculate the curvature M_{ki}:
!!$!------------------------------------
!!$!     first get the partition functions for every beta value:
!!$      print *,'get the partition functions:'
!!$      call part_funct(nbin,nbeta,bin_energy,dbeta,wbar,zbar)
!!$
!!$!     precalculate sum_{k=1}^K n_k/g_k exp^{-beta_k*S+f_k}:
!!$!     loop over all the energy bins:
!!$!      print *,'sum_{k=1}^K n_k/g_k exp^{-beta_k*S+f_k}:'
!!$!      print *,'g=',g
!!$      do i=0,nbin
!!$         denom=0.
!!$         do k=1,nbeta
!!$!            print *,'i,j=',i,j
!!$            aux=exp(-dbeta(k)*bin_energy(i)+fe(k))
!!$            denom=denom+aux*nmeas(k)/gtau(k)
!!$         enddo
!!$         norm(i)=denom
!!$!         print *,'i,norm(i)=',i,norm(i)
!!$      enddo
!!$
!!$!     for every k calculate...
!!$      do k=1,nbeta
!!$!     ...the derivative with respect to f_i:
!!$         do i=1,nbeta
!!$!            print *,'k,i=',k,i
!!$!     loop over all the energy bins:
!!$            aux=0.
!!$            do j=0,nbin
!!$!               print *,'loop over energy bins, j=',j
!!$               nom=exp(-dbeta(i)*bin_energy(j)+fe(i))
!!$               nom=nom*nmeas(i)/gtau(i)
!!$               aux=aux+exp(-dbeta(k)*bin_energy(j))*wbar(j)*nom/norm(j)
!!$            enddo
!!$            mki(k,i)=-aux/zbar(k)
!!$            if(i.eq.k) then
!!$               mki(k,i)=mki(k,i)+1.
!!$            endif
!!$         enddo
!!$      enddo
!!$
!!$      return
!!$      end
!!$
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!--------------------------------------------------------------------------
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      subroutine part_funct(nbin,nbeta,bin_energy,dbeta,wbar,zbar)

!     This routine calulates the partition functions for all beta values
!     with a given set of free energy values {f_i} (implicitly contained
!     in wbar).

      use fs_multi_par

      implicit none

!     number of energy bins and beta values:
      integer nbin,nbeta

!     the density of states:
      real(kind=dp) :: wbar(0:nbin)

!     the energy of the bins:
      real(kind=dp) :: bin_energy(0:nbin)
!     the beta differences:
      real(kind=dp) :: dbeta(1:nbeta)

!     the partition functions:
      real(kind=dp) :: zbar(1:nbeta)


!     auxiliary variables:
      integer i,k
      real(kind=dp) :: aux

!     loop over the beta values:
      do k=1,nbeta
!     loop over the energies:
         aux=0.
         do i=0,nbin
            aux=aux+wbar(i)*exp(-dbeta(k)*bin_energy(i))
         enddo
         zbar(k)=aux
      enddo


      return
    end subroutine part_funct

      
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!--------------------------------------------------------------------------
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine get_autocor(nbeta,nmeas,energy,obs,tau)  

!     This subroutine calculates the autocorrelation time for the different
!     MC runs from (the measured energies) and the observables.

      use fs_multi_par

      implicit none

!     number of beta values:
      integer nbeta

!     counts the number of measurements(nmeas(0) being the total number
!     of measurements):
      integer nmeas(0:nbeta)

!     the autocorrelation times:
      real(kind=dp) :: tau(1:nbeta)

!     the action and the observable:
      real(kind=dp) :: energy(1:max_meas,1:nbeta),obs(1:max_meas,1:nbeta)

!     the autocorrelation functions:
      real(kind=dp) :: acf(1:max_meas)

!     auxiliary variables:
      integer k,i,j
      real(kind=dp) :: auxtau,vev,vev2,denom,cor,vev_i,vev_j,auxi,auxj
      real(kind=dp) :: tauerr

      integer max_j

      open(67,file='acf.plo',form='formatted',status='unknown')
      print *
      print *,'Calculating the autocorrelation times:'

!     loop over the beta values:
      do k=1,nbeta
         print *,'k=',k
!     calculate the denominator of the autocorrelation function:
         vev=0.
         vev2=0.
         do i=1,nmeas(k)
!            vev=vev+obs(i,k)
            auxi=dble(energy(i,k))
            vev=vev+auxi
!            vev2=vev2+obs(i,k)**2
            vev2=vev2+auxi**2
         enddo
         vev=vev/float(nmeas(k))
         vev2=vev2/float(nmeas(k))
         denom=vev2-vev**2

!     loop over the MC time differences:
         max_j=min(nmeas(k)-1,5000)
         do j=1,max_j             !nmeas(k)-1
            cor=0.
            vev_i=0.
            vev_j=0.
            do i=1,nmeas(k)-j
!               cor=cor+obs(i,k)*obs(i+j,k)
!               vev_i=vev_i+obs(i,k)
!               vev_j=vev_j+obs(i+j,k)
               auxi=dble(energy(i,k))
               auxj=dble(energy(i+j,k))
               cor=cor+auxi*auxj
               vev_i=vev_i+auxi
               vev_j=vev_j+auxj
            enddo
            cor=cor/float(nmeas(k)-j)
            vev_i=vev_i/float(nmeas(k)-j)
            vev_j=vev_j/float(nmeas(k)-j)
!     construct the autocorrelation function:
!            acf(j)=(cor-vev**2)/denom
            acf(j)=(cor-vev_i*vev_j)/denom
         enddo
!     write it out:
         do j=1,max_j
            write(67,*) j,acf(j)
         enddo

!     calculate the integrated autocorrelation time:
         auxtau=0.5
!         do j=1,max_j           !nmeas(k)-1
         max_j=min(nmeas(k)-1,3000)
         do j=1,max_j
            if(acf(j).lt.0.) goto 491
            auxtau=auxtau+acf(j)*(1.-j/float(max_j))
         enddo

 491     continue

         tau(k)=auxtau



         auxtau=0.5
         do j=1,max_j
          ! this is the correct one, but the correction factor in parantheses
          ! can and is usually dropped:
          ! auxtau=auxtau+acf(j)*(1.-j/float(max_j))
          ! so we simply use:
            auxtau=auxtau+acf(j)

          ! this implements the self-consistent truncation 
          ! j_max >= 6 \tau_int(j_max):
            if(j >= 6.0*auxtau .or. j>=max_j) then
               tauerr=auxtau*sqrt(2.0*(2*j+1)/float(nmeas(k)))
               print *
               print '("------------------------------------------- ")'
               print '("Int autocorr time from self-consistent ")'
               print '("truncation at k_max =",i5,":")',j
               print *
               print '("  tau_int = ",f5.1," +/- ",f5.1)',auxtau,tauerr
               print *
               print '("------------------------------------------- ")'
               exit
            endif
         enddo
         tau(k)=auxtau


!     end of loop over beta values:
      enddo
      
      print '("tau:",12f12.5)',tau

      close(67)

      return
    end subroutine get_autocor


!=================================================================
      real(kind=dp) function ranf()
!=================================================================
!     for ranlux generator:

      use ranlxd_generator

      real(kind=dp), dimension(1:1) :: r
      
      call ranlxd(r)
      ranf=1.0_dp-r(1)
!      write(*,'("DEBUG: r(1)=",f24.16)') r(1)
!      stop
      return
    end function ranf

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!--------------------------------------------------------------------------
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  end program fs_multi
 
