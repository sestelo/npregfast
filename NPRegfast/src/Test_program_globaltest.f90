
program hola
implicit none
!integer,parameter::kbin=100,ncmax=5,nboot=10
integer i,p,nf,fact(100),iostat,iend,n,ikernel,iopt,model,kernel,nh,ipredict,r
integer,allocatable::F(:)
integer,allocatable::nc(:)
double precision,allocatable:: X(:),Y(:),w(:),h(:),C2(:,:),Xb(:),&
pb(:,:,:),li(:,:,:),ls(:,:,:),Dif(:,:,:,:),Difi(:,:,:,:),Difs(:,:,:,:),&
difc(:,:,:),difci(:,:,:),difcs(:,:,:),pboot(:,:,:,:),&
c(:,:),Cboot(:,:,:),a(:),ainf(:),asup(:),b(:),binf(:),bsup(:),&
predict(:,:,:),predictu(:,:,:),predictl(:,:,:),pcmin(:),pcmax(:)
double precision T,pvalor,D,Ci,Cs
integer kbin,ncmax,nboot
character*30 fichero






ncmax=5



!write (*,*) 'nombre del archivo de datos'
!read (*,*) fichero
!fichero="alo.txt"
!fichero="decus_phill_1.txt"
!fichero="prueba.txt"
!fichero="pediatria3.txt"
fichero="Datos_año.prn"
n=0
open (1,file=fichero)
	read (1,*,iostat=iend)
1   continue
    read (1,*,iostat=iend)
	if (iend.eq.0) then
		n=n+1
		goto 1
	end if		
close(1)



allocate (F(n),X(n),Y(n),W(n),h(100))

W=1
open (1,file=fichero)
	read (1,*)
	do i=1,n	
		read(1,*) F(i),X(i),Y(i)
	end do
close (1)





call Factores(F,n,fact,nf)

allocate (nc(nf),xb(kbin),pb(kbin,3,nf),pcmin(nf),pcmax(nf))

C2=-1
iopt=1
ikernel=1
nc=-1

if(ikernel.eq.1) p=1

!h=98 !80




!model=1 noparametrico
!model=2 alometrico
!model=0 contraste









nh=30
kernel=1


! Contraste
!**************************
! r=0  estimacion
! r=1 primera derivada
! r=2 segunda derivada

p=3
kbin=100
h=-1.0
r=1
pcmax=600
pcmin=0
nboot=100
call globaltest(F,X,Y,W,n,h,nh,p,kbin,fact,nf,kernel,nboot,r,T,pvalor)

!call localtest(F,X,Y,W,n,h,nh,p,kbin,fact,nf,kernel,nboot,pcmax,pcmin,r,D,Ci,Cs)
end
