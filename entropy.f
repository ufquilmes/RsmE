c
c
      program acomo
      implicit none
      integer i,j,k,l,n,n1,n2,n3,n4,nfrmax,ndat,ndat1,nume
      integer nsum2,nsum3,nfreq
      real*8 kb,temp,hpb,clight,voga
      parameter (nfrmax=200000,nfreq=5514)
      parameter (kb=1.3806488d-23) ! joule/K
      parameter (temp=310.0d0) ! K
      parameter (clight=2.9979245800d10) ! cm/s
      parameter (hpb=1.054571817d-34) !joule·s
      parameter (voga=6.02214076d23)      
      real*8 a1,val,asum(nfrmax,nfreq),sasum(nfrmax)
      real*8 ap(nfrmax),fsum
      real*8 ava(nfrmax,nfreq)
      real*8 avb(nfrmax,nfreq)
      real*8 amodu(nfreq)
      real*8 afreq(nfreq),entropy
      character*15 nada
      dimension nume(nfreq)
      character*4 exh


      open(75,file='frequencies.dat')


      entropy=0.0d0
     
      !read output from cpptraj with frequency of the modes obtained from the mwcovar covariance matrix 
       do k = 1, nfreq
        read(75,12) nume(k), afreq(k)  !viene en cm-1 en mwcovar
        afreq(k)=afreq(k)*clight        ! de cm-1 a s-1
        entropy=entropy+((hpb*afreq(k)/(kb*temp))/(dexp(hpb*afreq(k)/
     &  (kb*temp))-1))-dlog(1-dexp(-hpb*afreq(k)/(kb*temp)))
       enddo
       entropy=entropy*kb*voga !en J/Kmol
       write(6,22) entropy,entropy/4.184d0,temp*entropy/4.184d3      !entropia en J y cal, TS kcal/mol
       
      close (75)


   12 format(2x,i3,2(2x,f14.5))
   22 format(3(2x,f14.7))


      return
      end
