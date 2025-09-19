module fs_multi_par

!     maximal number of bins:
      integer :: max_nbin
      parameter(max_nbin=500)

!     number of bootstrap samples:
      integer :: nbtrp
      parameter(nbtrp=1000)

!     bootstrap sample size:
      integer :: blocksize
      parameter(blocksize=10)

!     maximal number of measurements:
      integer :: max_meas
      parameter(max_meas=500001)

!     maximal number of beta values:
      integer :: max_nbeta
      parameter(max_nbeta=10)

!     double precision
      integer,parameter :: dp=selected_real_kind(14)
      
end module fs_multi_par




