module fs_multi_par

!     maximal number of bins:
      integer :: max_nbin
      parameter(max_nbin=500)

      !     number of bins for observable:
      integer :: obs_nbin
      parameter(obs_nbin=40)

!     number of bootstrap samples:
      integer :: nbtrp
      parameter(nbtrp=1000)

!     bootstrap sample size:
      integer :: blocksize
      parameter(blocksize=100)

!     maximal number of measurements:
      integer :: max_meas
      parameter(max_meas=500001)

!     maximal number of beta values:
      integer :: max_nbeta
      parameter(max_nbeta=20)

!     double precision
      integer,parameter :: dp=selected_real_kind(14)
      
end module fs_multi_par




