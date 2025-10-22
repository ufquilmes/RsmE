************************************************************************
      program umbrella
************************************************************************
      implicit double precision (a-h,o-z)
      parameter (nbinf=3)     
      parameter (ncajas=30)   !parametro para la funcion del error  
      parameter (ncenmax=680, nbinmax=nbinf*ncenmax, nmcmax=1000000)
      parameter (lwork= 10*nbinmax) !parametro para subrutina dham
      parameter (temp=303.0d0)
      parameter (gas_c=1.987d-3)
      parameter (crit_conv=1.0d-6)
      parameter (akt=temp*gas_c)  !RT  in kcal/mol
      parameter (betha=1.0d0/akt)  !1/RT  in kcal/mol
      dimension xk(ncenmax),xi(nbinmax),xp(ncenmax,nmcmax) !xp son los datos
      dimension xpmink(ncenmax),xpmaxk(nbinmax) !min y max de cada ventana
      dimension ncont(nbinmax,ncenmax),ncontk(nbinmax)
      dimension ntmat(nbinmax,nbinmax),ammat(nbinmax,nbinmax) !para dham
      dimension amark(nbinmax,nbinmax) !para dham
      dimension vl(nbinmax,nbinmax),vr(nbinmax,nbinmax) !dham diag
      dimension wr(nbinmax),wi(nbinmax) !dham diag
      dimension work(lwork) !dham diag
      dimension nj1(nmcmax),nj0(nmcmax)
      dimension pn1(nmcmax),pn0(nmcmax)
      dimension funct(nmcmax),du10(nmcmax) !error good practice
      dimension centro(nmcmax),u1(nmcmax),u0(nmcmax) !error en good practice
      dimension xf_min(ncenmax),xf_max(ncenmax),bar_id(ncenmax)
      dimension up(nbinmax,ncenmax)
      dimension prb(nbinmax),prb0(nbinmax) !prob wham
      dimension hark(ncenmax),fact(ncenmax)
      dimension filename(ncenmax)
      dimension ntot(nbinmax),ntotk(ncenmax)
      character*11 filename
      character*4 exh
      character*1 jobvl,jobvr !dham diag
      INTEGER givename
      CHARACTER*200 str,nombre
      LOGICAL :: wham,dham,bar,errores,angulo


      !si lo de abajo es 0 no hay input y toma default . no hacer dham y si errores
      givename=iargc()

      if (givename.eq.0) then

         wham=.TRUE.
         niter=500000
         dham=.FALSE.
         bar=.TRUE.
         nbariter=150
         errores=.FALSE.
         angulo=.FALSE.

      else

         call getarg(1,nombre)
         open(1,file=nombre,STATUS='OLD',ACCESS='SEQUENTIAL')
         read(1,*) str
         read(1,*) wham
         read(1,*) str
         read(1,*) niter
         read(1,*) str
         read(1,*) dham
         read(1,*) str
         read(1,*) bar
         read(1,*) str
         read(1,*) nbariter
         read(1,*) str
         read(1,*) errores
         read(1,*) str
         read(1,*) angulo
         close(1)

      endif

      ndump=int(niter/5) !cada cuanto escribe el temp el wham

      !read input file with files, equilibrium position and constant
      open(999,file="wham.inp")
      api=acos(-1.0d0)  

      do k=1,ncenmax
        read(999,*,end=300) filename(k),xk(k),hark(k)
        if (angulo) hark(k)=hark(k)/(180.0d0/api)**2    
      enddo
  300 continue
      ncenter= k-1
      close(999) 

      !stop if #windows is bigger than maxvalue
      if (ncenter.eq.ncenmax) then
        write(6,*) "Muchas ventanas pa: modifica el codigo"
        stop
      endif

      !loop read sobre las k simulaciones
      do k=1,ncenter
        open(998,file=filename(k))
        do j=1,nmcmax
          read(998,*,end=500) ndumb,xp(k,j)
        enddo
  500   continue
        close(998)
        ntotk(k)= j-1
        if (ntotk(k).eq.nmcmax) then
          write(6,*) "Muchos datos en la ventana=",k
          stop
        endif
      enddo

      !busca el minimo y el máximo observado en xp
      !y el minimo y el maximo en cada ventana

      do k=1,ncenter
        xpmink(k)= 1.0d3   
        xpmaxk(k)= -xpmink(k)   
        do j=1,ntotk(k)
          if (xp(k,j).lt.xpmink(k)) then
            xpmink(k)=xp(k,j)
          endif
          if (xp(k,j).gt.xpmaxk(k)) then
            xpmaxk(k)=xp(k,j)
          endif
        enddo
      enddo
      
      xlowest=xpmink(1)
      xhighest=xpmaxk(ncenter)

      !warning por si un boludo tiene el min/max en otro lugar 
      do k=2,ncenter-1
        if (xpmink(k).lt.xlowest.or.xpmaxk(k).gt.xhighest) then
        write(6,*) "donde tenes el minimo o maximo pa?"
        STOP
        else
        endif
      enddo


!###############################################################################
!#################calculo de error conjunto con BAR previo a wham y dham  ######
!###############################################################################
        if (.not.errores) then  
            if (.not.bar) go to 114  
        else 
        endif         
  
        open(699,file='bar.out')
        !contiene la diferencia de energia entre distribuciones y la acumulada 
        !para la curva !!!OJO que la acumulada va en kcal/mol y la otra en KbT units
         
        bar_acu=0.0d0
        !la funcion es de a pares por eso empieza en 2. 2-1,3-2,4-3,etc
        nnom=0
        ncorrec_tail=6

        do n=2,ncenter

         bar_id(n)=0.0d0
         nnom=nnom+1 
         if (ntotk(n-1).ne.ntotk(n)) then
          write(6,*) "Calculo del error: esto no anda"
          write(6,*) "Habria que dividir en dos el loop"
          STOP
         endif

         if (.not.errores) go to 112  


         open(99,file='func.'//filename(n))


         do i= 1+ncorrec_tail,ncajas-ncorrec_tail

          nj1(i)=0
          nj0(i)=0
          xf_min(i)= xk(n-1)+(i-1)*(xk(n)-xk(n-1))/ncajas
          xf_max(i)= xk(n-1)+i*(xk(n)-xk(n-1))/ncajas

          do j=1,ntotk(n)
           if (xp(n-1,j).gt.xf_min(i).and.xp(n-1,j).le.xf_max(i)) then
           nj0(i)=nj0(i)+1
           endif
           if (xp(n,j).gt.xf_min(i).and.xp(n,j).le.xf_max(i)) then
           nj1(i)=nj1(i)+1
           endif
          enddo
          pn1(i)= dfloat(nj1(i))/dfloat(ntotk(n))
          pn0(i)= dfloat(nj0(i))/dfloat(ntotk(n))

          centro(i)=xf_min(i)+0.5d0*(xf_max(i)-xf_min(i))
          u1(i)=(hark(n)*0.5d0)*(centro(i)-xk(n))**2
          u0(i)=(hark(n-1)*0.5d0)*(centro(i)-xk(n-1))**2
          du10(i)=u1(i)-u0(i)
          funct(i)= log(pn1(i)/pn0(i))+betha*du10(i)
          write(99,25) du10(i),funct(i)

         enddo


  112 continue
        !por si hay dos constantes diferente o diferente #puntos, no va a andar    
         if (hark(n-1).ne.hark(n)) then
          write(6,*) "esto no anda. pone la misma k"
          STOP
         endif
       

          if (nnom.eq.1) df=0.0d0
          !Eq 18 Kim 2012 JPC
          !M es cero si ntot0=ntot1

          do it=1,nbariter

           sumo0=0.0d0
           sumo1=0.0d0
           do j=1,ntotk(n)
            ubar1=(hark(n-1)*0.5d0)*(xp(n-1,j)-xk(n))**2
            ubar0=(hark(n-1)*0.5d0)*(xp(n-1,j)-xk(n-1))**2
            du=ubar1-ubar0
            sumo0=sumo0+1.0d0/(1.0d0+dexp(betha*(-df+du)))
            ubar1=(hark(n)*0.5d0)*(xp(n,j)-xk(n))**2
            ubar0=(hark(n)*0.5d0)*(xp(n,j)-xk(n-1))**2
            du=ubar1-ubar0
            sumo1=sumo1+1.0d0/(1.0d0+dexp(betha*(df-du)))
           enddo
           df=(1.0d0/(betha))*dlog(sumo1/sumo0)+df
            if (mod(it,nbariter).eq.0) then
             bar_acu=bar_acu+(betha*df)*akt
             bar_id(n)=bar_acu
             if (.not.bar) then
             else
             write (699,25) 0.5d0*(xk(n)+xk(n-1)),bar_acu,betha*df
             endif   
             if (.not.errores) then
             else
             write (99,*) "#",xk(n),bar_acu,betha*df   
             close(99)  
             endif   
            endif

          enddo 

  113 continue

        enddo

        close(699) 

  114 continue
       
!###################################################################################
!###################################################################################
!#################        fin calculo de error y BAR            ####################
!###################################################################################
!###################################################################################

!###########comun a wham y dham############
!###########comun a wham y dham############

      !number of bins as times the centers. For wham
      nbin=nbinf*ncenter

c intervalo entre centros de cajas
      delbin= (xhighest-xlowest)/dfloat(nbin)
c pone los centros de las cajas donde calculara el PMF
      xbmin= xlowest+0.5d0*delbin
      do i=1,nbin
        xi(i)= xbmin+delbin*dfloat(i-1)
      enddo
c calculo de los u_i^(k)  !center es el numero de simulaciones
      do k=1,ncenter
        do i=1,nbin
          up(i,k)= 0.5d0*hark(k)*(xi(i)-xk(k))**2
        enddo
      enddo


c calculamos poblaciones y matriz de transicion
!ntmat=0 hace cero toda la matriz 
 
      ncont= 0
      ntmat= 0
      do k=1,ncenter
        xx= xp(k,1)
        ibin0= int((xx-xlowest)/delbin)+1
        ncont(ibin0,k)= ncont(ibin0,k)+1
        do j=2,ntotk(k)
          xx= xp(k,j)
          ibin= int((xx-xlowest)/delbin)+1
          ncont(ibin,k)=ncont(ibin,k)+1
          ntmat(ibin,ibin0)=ntmat(ibin,ibin0)+1
          ibin0= ibin
        enddo
      enddo

!########## fin comun a dham y wham #########

      if (.not.dham) go to 122  

!######################################DHAMMMMMMMMMMMMMMMMMMMMMMMM##############
!######################################DHAMMMMMMMMMMMMMMMMMMMMMMMM##############
!######################################DHAMMMMMMMMMMMMMMMMMMMMMMMM##############
!######################################DHAMMMMMMMMMMMMMMMMMMMMMMMM##############
!######################################DHAMMMMMMMMMMMMMMMMMMMMMMMM##############
!######################################DHAMMMMMMMMMMMMMMMMMMMMMMMM##############
!######################################DHAMMMMMMMMMMMMMMMMMMMMMMMM##############


c calcula matriz de markov sin normalizar
      do i=1,nbin
        do j=1,nbin
          deno= 0.0d0
          do k=1,ncenter
            if (ncont(i,k).ge.1) then
              expo=0.5d0*(up(j,k)-up(i,k))/akt
              deno=deno+ dfloat(ncont(i,k))*dexp(-expo)
            endif
          enddo
          if (ntmat(j,i).eq.0) then
            ammat(j,i)= 0.0d0
          else
            ammat(j,i)= dfloat(ntmat(j,i))/deno
          endif
        enddo
      enddo

c calcula la matriz de markov normalizada
      do i=1,nbin
        suma= 0.0d0
        do j=1,nbin
          suma=suma+ammat(j,i)
        enddo
        do j=1,nbin
          if (suma.gt.0.0d0) then
            amark(j,i) =ammat(j,i)/suma
          else
            amark(j,i)= ammat(j,i)
            write(6,*) "suma de ammat=0"
            write(6,*) "no deberia pasar por aca"
            write(6,*) "hay un error"
            STOP
          endif
        enddo
      enddo

c calculo de dham
      jobvl= "n"
      jobvr= "v"
      lda= nbinmax
      ldvl= nbinmax
      ldvr= nbinmax
      info= 0
      call dgeev(jobvl,jobvr,nbin,amark,lda,wr,wi,vl,ldvl,vr,ldvr,
     #           work,lwork,info)

c buscamos autovector con autovalor mas cercano a 1.0d0
      eigv= 100.0d0
      do i=1,nbin
        test= dabs(wr(i)-1.0d0) 
        if (test.lt.eigv) then
          eigv= test
          imin= i
        endif
      enddo
      if (vr(1,imin).lt.0.0d0) vr=-vr

      write(6,*) "Autovector",imin,"Autovalor",wr(imin)

c buscamos el minimo de energia libre del pozo de la izquierda
      nlim= nbin/2
      pmax= 0.0d0
      do i=1,nlim
        if (vr(i,imin).gt.pmax) then
          pmax= vr(i,imin)
        endif
      enddo
      const= -akt*dlog(pmax)
      open(15,file="dham.out")
      do i=1,nbin
        write(15,99) xi(i),-akt*log(vr(i,imin))-const,vr(i,imin),i
      enddo
      close(15)


!############################FIN#######DHAMMMMMMMMMMMMMMMMMMMMMMMM##############
!############################FIN#######DHAMMMMMMMMMMMMMMMMMMMMMMMM##############
!############################FIN#######DHAMMMMMMMMMMMMMMMMMMMMMMMM##############
!############################FIN#######DHAMMMMMMMMMMMMMMMMMMMMMMMM##############
!############################FIN#######DHAMMMMMMMMMMMMMMMMMMMMMMMM##############


  122  continue  


!######################################WHAMMMMMMMMMMMMMMMMMMMMMMMM##############
!######################################WHAMMMMMMMMMMMMMMMMMMMMMMMM##############
!######################################WHAMMMMMMMMMMMMMMMMMMMMMMMM##############
!######################################WHAMMMMMMMMMMMMMMMMMMMMMMMM##############
!######################################WHAMMMMMMMMMMMMMMMMMMMMMMMM##############
!######################################WHAMMMMMMMMMMMMMMMMMMMMMMMM##############
!!!!!!!! energia libre por WHAM
!cantidad por bin, ntot

      if (.not.wham) go to 123  

      do i=1,nbin
        ntot(i)=0
        do k=1,ncenter
          ntot(i)=ntot(i)+ncont(i,k)
        enddo
      enddo


c loop sobre iteraciones del wham
c inicializacion cabeza de perro

      do i=1,nbin
      prb0(i)= 1.0d0/dfloat(nbin)
      enddo

      do iter=1,niter
        do k=1,ncenter
          suma= 0.0d0
          do i=1,nbin
            suma= suma+prb0(i)*dexp(-up(i,k)/akt)
          enddo
          fact(k)= 1.0d0/suma
        enddo
 
        do i=1,nbin
          deno= 0.0d0
          do k=1,ncenter
            deno= deno+ntotk(k)*fact(k)*dexp(-up(i,k)/akt)
          enddo
          prb(i)= dfloat(ntot(i))/deno
        enddo

        resto= 0.0d0
        suma= 0.0d0
        do i=1,nbin
          resto= resto+(dlog(prb(i))-dlog(prb0(i)))**2
          suma= suma+prb(i)
        enddo
        resto=dsqrt(resto)
        if (resto.le.crit_conv) goto 700 !finish iteration process if fulfill


 
c         bajar resultados temporarios
          if (mod(iter,ndump).eq.0) then
 
           nsum=int(int(iter)/ndump)
           exh=char(48+nsum)
 
           open(86,file='temp_wham.'//exh)
 
           pmax= 0.0d0
           do i=1,nbin
             if (prb(i).gt.pmax) then
               pmax= prb(i)
             endif
           enddo
           const= -akt*dlog(pmax)
 
           open(86,file='temp_wham.'//exh)
           do j=1,nbin
             write(86,99) xi(j),-akt*dlog(prb(j))-const,prb(j),j
           enddo
           write(86,86) "#Iter=",iter,"Resto=",resto
           close(86)

          endif
          !!!!fin de baja de datos temporarios.
   
        !actualizamos las probabilidades
        do i=1,nbin
         prb0(i)= prb(i)
        enddo

      enddo


  700 continue


      
      pmax= 0.0d0
      do i=1,nbin
        if (prb(i).gt.pmax) then
          pmax= prb(i)
        endif
      enddo
      const= -akt*dlog(pmax)

      open(16,file="wham.out")
      do j=1,nbin
        write(16,99) xi(j),-akt*dlog(prb(j))-const,prb(j),j
      enddo
      write(16,87) "#Minimo=",xlowest, "#Maximo=",xhighest
      write(16,87) "#Resto= ",resto
      close(16)
      

  123  continue  

!#################################################################
!################### otro calculo de error: test de kuber   ######
!#################################################################
      if (.not.errores) go to 124 

!     calculo de las entropias relativas
      open(500,file="entropias1.dat")
      do k=1,ncenter
        suma= 0.0d0
        do i=1,nbin
          pki= fact(k)*dexp(-up(i,k)/akt)*prb(i)
          fact1= dfloat(ncont(i,k))/dfloat(ntotk(k))
          if (fact1.gt.0.0d0) then
            fact2=dlog(fact1/pki)
            suma=suma+fact1*fact2
          endif
        enddo
        write(500,99) xk(k),suma
      enddo
      close(500)

      open(600,file="entropias2.dat")
      do k=1,ncenter
        suma= 0.0d0
        do i=1,nbin
          pki= fact(k)*dexp(-up(i,k)/akt)*prb(i)
          fact1= dfloat(ncont(i,k))/dfloat(ntotk(k))
          if (fact1.gt.0.0d0) then
            fact2= dlog(pki/fact1)
            suma= suma+ pki*fact2
          endif
        enddo
        write(600,99) xk(k),suma
      enddo
      close(600)

  124  continue  


!################################### fin calculo del error ######
!################################################################ 

  25  format(6(2x,f14.7))
  86  format(a6,2x,i10,2x,a6,2x,f16.14)
  87  format(2(a8,2x,f16.14))
  99  format(1x,f8.3,1x,f16.10,1x,f18.14,1x,i3)


      stop
      end
