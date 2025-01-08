      subroutine SRBFestimation(observationfl,surfhgtgrdfl,checkpointfl,para)
      implicit none
	character*800::observationfl,surfhgtgrdfl,checkpointfl,outfl,nstr(80)
	character*80000::line,line0,str,astr
      character(len=15)::znm(900000)
	integer sn,kln,astat(20),i,j,k,NF,NF4,Kt,nn,mm,obsn,lvl,order,ks,kp,kpnt,k0
      real*8 rec(8000),BLH(3),rln(3),rlnk(3),gr,NFD(5),dln(2),GMr,rhd(4),rr,sigmax,sigma0,prst(2)
	real*8 para(20),GRS(6),hd(6),GM,ae,pi,RAD,wgh,r0,dpth,nta,br,sta(4)!br-Bjerhammar
	real*8 val,unit(5),dlat,dr,tmp,mr,st(4),st0(4),bf(8),blat,sinf,tt,cosa,sina,lmt!第一行平行圈格网中心地心纬度°
	integer nlon,nlat,minN,maxN,nk,hrw,obsrw,ob2rw,wghrw,kobs,fknd,krbf
	integer ki,kj,nd,mk,ni,nj,mthd,edgn
	integer::status=0
	real*8,allocatable::mpn(:,:),mdp(:,:),obs(:,:),dl(:),lon(:,:),RBF4(:,:),RBFn(:,:)
      real*8,allocatable::BPB(:,:),BPL(:),BB(:),B2(:),xx(:),chs(:),sr(:),rlatlon(:,:),chd(:,:)
	real*8,allocatable::mpn4(:,:),mdp4(:,:),hgt(:,:),rst(:,:),rst2(:,:),RBF(:,:)!目标场元
	integer,allocatable::nln(:),nrd(:,:),node(:),enode(:),gpnt(:,:)!格网、未知数序号,每个观测量有效节点序号
!---------------------------------------------------------------------
      GRS(1)= 3.986004415d14; GRS(2)=6378137.d0; GRS(3)=1.0826359d-3
      GRS(4) = 7.292115d-5; GRS(5)=1.d0/298.25641153d0;ae=GRS(2)
      GM=GRS(1)*1.d-7;ae=GRS(2)
      pi=datan(1.d0)*4.d0;RAD=pi/180.d0;mr=36.d2/RAD
      unit(1)=1.d5;unit(2)=1.d0;unit(3)=1.d5;unit(4)=1.d9;unit(5)=mr
      minN=nint(para(1));maxN=nint(para(2));order=nint(para(3));krbf=nint(para(4))
      lvl=nint(para(5));kobs=nint(para(8));fknd=nint(para(9));mthd=nint(para(13))
      obsrw=nint(para(10));ob2rw=nint(para(11));wghrw=nint(para(12))
      dr=para(6)*1.d3/ae/RAD!球面角距°
      dpth=para(7)*1.d3
      kpnt=1
      if(minN<order)minN=order
      !打开计算面大地高格网文件
      open(unit=8,file=surfhgtgrdfl,status="old",iostat=status)
      if(status/=0) goto 902
      read(8,'(a)') line0
      call PickReclong(line0,kln,rec,sn)
      if(sn<6)then
         close(8);goto 902
      endif
      hd(1:6)=rec(1:6)
	nlat=nint((hd(4)-hd(3))/hd(6))
	nlon=nint((hd(2)-hd(1))/hd(5))
	hd(5)=(hd(2)-hd(1))/dble(nlon)
	hd(6)=(hd(4)-hd(3))/dble(nlat)
 	allocate(hgt(nlat,nlon), stat=astat(1))
 	allocate(rst(nlat,nlon), stat=astat(2))
 	allocate(rst2(nlat,nlon), stat=astat(3))
	if (sum(astat(1:3)) /= 0) then
          close(8);goto 902
      endif
 	do i=1,nlat
	   read(8,*,end=905)(hgt(i,j),j=1,nlon)
      enddo
905   close(8)
      k0=0; k=0
      open(unit=8,file=observationfl,status="old",iostat=status)
      if(status/=0) goto 901
      read(8,'(a)') line !头文件一行
      do while(.not.eof(8))
        read(8,'(a)') line
        call PickRecstr(line,kln,nstr,sn)
        if(sn>3)then
          k=k+1;znm(k)=nstr(1)
        endif
        call PickReclong(line,kln,rec,sn)
        if(sn>3)then
            if(rec(2)<hd(2).and.rec(2)>hd(1).and.rec(3)<hd(4).and.rec(3)>hd(3))k0=k0+1
        endif
      enddo
      close(8);obsn=k!观测数
      if(k0<6)goto 901
      NF=nint(dr*3600)!影响半径等分,间隔1″,NF+1→[0,dr]
      NF4=nint(dr*600)!间隔6″,NF4+1→[0,dr]
 	allocate(mpn(maxN-minN+1,NF+1), stat=astat(1))
 	allocate(mdp(maxN-minN+1,NF+1), stat=astat(2))
 	allocate(mpn4(maxN-minN+1,NF4+1), stat=astat(3))
 	allocate(mdp4(maxN-minN+1,NF4+1), stat=astat(4))
 	allocate(obs(4*obsn,6), stat=astat(5))!地心距，地心纬度，经度，权值，观测量1(，观测量2)
	if (sum(astat(1:5)) /= 0) goto 901
      obs(:,6)=0.d0;r0=0.d0;k=0!重新读取观测文件，转换为观测点球面坐标，计算平均地心距r0
      open(unit=8,file=observationfl,status="old",iostat=status)
      read(8,'(a)') line
      do while(.not.eof(8))
        read(8,'(a)') line
        call PickReclong(line,kln,rec,sn)
        if(sn<obsrw)goto 306
        BLH(1)=rec(3);BLH(2)=rec(2);BLH(3)=rec(4);wgh=1.d0
        if(wghrw>0.and.wghrw<sn+1)wgh=rec(wghrw);if(wgh<-1.d-5)wgh=1.d0
        call BLH_RLAT(GRS,BLH,rln)
        k=k+1;r0=r0+rln(1)/ae
        obs(k,1:3)=rln(1:3);obs(k,4)=wgh;obs(k,5)=rec(obsrw)
        if(ob2rw>4.and.ob2rw<sn.and.kobs==4)then
          k=k+1;obs(k,1:3)=rln(1:3);obs(k,4)=wgh;obs(k,5)=rec(ob2rw);obs(k,6)=9.d4
        endif
306     continue
      enddo
      close(8);r0=r0/dble(obsn)*ae!平均地心距
      obsn=k
   !计算格网节点/未知数个数Kt
      !nln(nn)-平行圈方向格网数nn=maxi-mini+1
      !sr(nn)平行圈方向格网面积与赤道格网面积之差的百分比
      !dl(nn)平行圈方向经度间隔°Kt格网总点数-节点数、未知数
      !lon(nn,mm)格网中心经度,mm为平行圈方向最多格网数
      dlat=180.d0/lvl;nd=nint(dr/dlat+0.5d0)!dlat格网间隔,积分半径对应的格网数
      rhd(1:4)=hd(1:4);BLH(2)=(hd(1)+hd(2))/2.d0;BLH(3)=0.d0!!!!!!!目标格网范围用球坐标表示
      BLH(1)=hd(3);call BLH_RLAT(GRS,BLH,rln);rhd(3)=rln(2)
      BLH(1)=hd(4);call BLH_RLAT(GRS,BLH,rln);rhd(4)=rln(2)!!!!!!!!!
      nn=nint((rhd(4)-rhd(3))/dlat+0.5);mm=nint((rhd(2)-rhd(1))/dlat+0.5)!mm平行圈方向最大格网数
      allocate(nln(nn),sr(nn),dl(nn),chs(obsn),nrd(nn,mm),gpnt(nn,mm),lon(nn,mm),rlatlon(2*(nn+mm),2),enode(2*(nn+mm)))
      call ReuterGrid(rhd,lvl,Kt,blat,nn,mm,nln,sr,dl,nrd,lon)!Kt节点数/未知数个数
      gpnt=0!计算格网中测点数，修正Reuter格网节点数Kt,序号nrd
      call Edgnode(enode,rlatlon,lvl,edgn,lon,blat,nln,gpnt,nn,mm)
      ks=edgn+obsn;allocate(chd(ks+obsn,3))
 	chd(1:ks,1)=r0;chd(1:edgn,2:3)=rlatlon(1:edgn,1:2)
      chd(edgn+1:edgn+obsn,1:3)=obs(1:obsn,1:3)
      call AdjReuterGrd(chd(1:ks,1:3),ks,Kt,blat,lvl,nn,mm,nln,dl,lon,nrd,gpnt)
 	allocate(BPB(Kt,Kt), stat=astat(1))
	if (sum(astat(1:1)) /= 0) goto 903
      allocate(RBF(NF+1,5),BPL(Kt),BB(Kt),B2(Kt),xx(Kt),node(Kt),RBF4(NF4+1,4),RBFn(maxN-minN+1,4))
   !计算minN~maxN阶NF+1组勒让德函数[minN~maxN]×[0,dr]
      call LegPn01(mpn,mdp,minN,maxN,NF,dr)
      call LegPn01(mpn4,mdp4,minN,maxN,NF4,dr)
   !由初始补偿深度dpth,计算SRFB曲线
      br=r0-dpth;nta=br/r0;rlnk(1)=br!初始补偿深度dpth和宽度参数nta
      call SRBF5all(RBF,order,krbf,mpn,mdp,minN,maxN,NF,nta)
      BPL=0.d0;BPB=0.d0!构造观测方程和法方程
      do k=1,obsn!0扰动重力，1高程异常，2空间异常，3扰动重力梯度，4垂线偏差
        rln(1:3)=obs(k,1:3);wgh=obs(k,4);val=obs(k,5)/unit(kobs+1);rr=rln(1)
        if(kobs==0.or.kobs==2)GMr=GM/rr/rr
        if(kobs==3)GMr=GM/rr/rr/rr
        if(kobs==1.or.kobs==4)then
           call RLAT_BLH(GRS,rln,BLH);call GNormalfd(BLH,NFD,GRS);gr=NFD(2)
           if(kobs==1)GMr=GM/rr/gr
           if(kobs==4)GMr=GM/rr/rr/gr
        endif
        call RtGridij(rln,ki,kj,blat,lvl,nn,mm,nln,dl,lon)!计算rln在格网中的序号ki>0,kj>0
        mk=0;BB=0.d0;node=0!BB-RBF系数,node(Kt)观测量作用半径范围内有效节点序号
        do i=ki-nd,ki+nd!增加实现节点序号与格网行列号对应，即节点序号二维整数数组
          if(i<1.or.i>nn)goto 1001!lon(i,j),blat第一行平行圈格网中心地心纬度°
          rlnk(2)=blat+(i-1.d0)*dlat
          do j=kj-nd,kj+nd
            if(j<1.or.j>nln(i)) goto 1002 
            if(nrd(i,j)<1) goto 1002 !nrd(i,j)未知数序号
            rlnk(3)=lon(i,j);call drln(rln,rlnk,dln)!距离m与夹角°
            if(dln(2)>dr)goto 1002
            call RBFvalue(RBF(:,kobs+1),NF,dr,dln(2),tmp)!由球面角距°dln(2)内插RBF值tmp
            mk=mk+1;node(mk)=nrd(i,j)
            BB(nrd(i,j))=GMr*tmp*dexp((minN-1.0)*dlog(r0/rr))
            if(kobs==4)then
              tt=dcos(dln(2)*RAD);sinf=dsqrt(1.d0-tt**2);if(sinf<1.d-12)sinf=1.d-12
              if(obs(k,6)<9.d4)then
	          cosa=(dcos(rln(2)*RAD)*dsin(rlnk(2)*RAD)-dsin(rln(2)*RAD)*dcos(rlnk(2)*RAD)*dcos(rlnk(3)*RAD-rln(3)*RAD))/sinf
                BB(nrd(i,j))= BB(nrd(i,j))*cosa
              else
	          sina=dcos(rlnk(2)*RAD)*dsin(rlnk(3)*RAD-rln(3)*RAD)/sinf
                BB(nrd(i,j))= BB(nrd(i,j))*sina
              endif
            endif
1002        continue
          enddo
1001      continue    
        enddo
        do i=1,mk
           ki=node(i)
           BPL(ki)=BPL(ki)+BB(ki)*val*wgh
           do j=1,i
             kj=node(j)
             BPB(ki,kj)=BPB(ki,kj)+BB(ki)*BB(kj)*wgh
	     enddo
        enddo
      enddo
      call Stat1d(obs(1:obsn,5),obsn,st0)
      tmp=0.d0
	do i=1,Kt
	   do j=1,i-1
	      BPB(j,i)=BPB(i,j)
         enddo
         tmp=tmp+BPB(i,i)**2/dble(Kt)
      enddo
      tmp=dsqrt(tmp)
      !以Reuter格网四周节点未知数为零组成观测方程，抑制边缘效应。
      !节点序号数组enode
      do i=1,edgn!edgn-Reuter格网四周节点数
         ki=enode(i); BPB(ki,ki)=BPB(ki,ki)+tmp/dsqrt(dble(obsn))
      enddo
	do i=1,Kt
         BPB(i,i)=BPB(i,i)+tmp*1.d-4
      enddo
5001  xx=0.d0
      !mthd=1 LU分解,2 Cholesky分解,3 最小二乘QR分解,4最小范数奇异值分解,5-岭估计
      if(mthd<5)call Equsolve(BPB,xx,Kt,BPL,mthd,bf)
      if(mthd==5)call RidgeEstimate(BPB,xx,Kt,BPL)
      chs=0.d0!计算残差并统计
      do k=1,obsn
        rln(1:3)=obs(k,1:3);rr=rln(1)
        if(kobs==0.or.kobs==2)GMr=GM/rr/rr
        if(kobs==3)GMr=GM/rr/rr/rr
        if(kobs==1.or.kobs==4)then
           call RLAT_BLH(GRS,rln,BLH);call GNormalfd(BLH,NFD,GRS);gr=NFD(2)
           if(kobs==1)GMr=GM/rr/gr
           if(kobs==4)GMr=GM/rr/rr/gr
        endif
        call RtGridij(rln,ki,kj,blat,lvl,nn,mm,nln,dl,lon)!计算rln在格网中的序号ki>0,kj>0
        mk=0;BB=0.d0;node=0!s-RBF系数,node(Kt)观测量作用半径范围内有效节点序号
        do i=ki-nd,ki+nd!增加实现节点序号与格网行列号对应，即节点序号二维整数数组，0表示!!!!!!!!!!!!!!!!!!
          if(i<1.or.i>nn)goto 2001!lon(i,j),blat第一行平行圈格网中心地心纬度°
          rlnk(2)=blat+(i-1.d0)*dlat
          do j=kj-nd,kj+nd!!!!!!!!!!!!!!!!!!!
            if(j<1.or.j>nln(i)) goto 2002 
            if(nrd(i,j)<1) goto 2002 !nrd(i,j)未知数序号
            rlnk(3)=lon(i,j);call drln(rln,rlnk,dln)!距离m与夹角°
            if(dln(2)>dr)goto 2002
            call RBFvalue(RBF(:,kobs+1),NF,dr,dln(2),tmp)!由球面角距°dln(2)内插RBF值val
            mk=mk+1;node(mk)=nrd(i,j)
            BB(nrd(i,j))=GMr*tmp*dexp((minN-1.0)*dlog(r0/rr))
            if(kobs==4)then
              tt=dcos(dln(2)*RAD);sinf=dsqrt(1.d0-tt**2);if(sinf<1.d-12)sinf=1.d-12
              if(obs(k,6)<9.d4)then
	          cosa=(dcos(rln(2)*RAD)*dsin(rlnk(2)*RAD)-dsin(rln(2)*RAD)*dcos(rlnk(2)*RAD)*dcos(rlnk(3)*RAD-rln(3)*RAD))/sinf
                BB(nrd(i,j))= BB(nrd(i,j))*cosa
              else
	          sina=dcos(rlnk(2)*RAD)*dsin(rlnk(3)*RAD-rln(3)*RAD)/sinf
                BB(nrd(i,j))= BB(nrd(i,j))*sina
              endif
            endif
2002        continue
          enddo
2001      continue    
        enddo
        val=0.d0
        do i=1,mk
          ki=node(i)
          val=val+BB(ki)*xx(ki)
        enddo
        chs(k)=obs(k,5)-val*unit(kobs+1)
      enddo
      call Stat1d(chs,obsn,st)
      open(unit=10,file="residuals.txt",status="replace")
	write(10,'(2F10.4,2F11.4,a,2F10.4,40F11.4)')(st0(i),i=1,4)," residuals:",(st(i),i=1,4)
      kp=1;if(kobs==4)kp=2
 	do i=1,obsn,kp
         call RLAT_BLH(GRS,obs(i,1:3),BLH)
	   if(abs(kobs-4)>0)write(10,'(a8,2F12.6,2F11.3,8F12.4)')trim(znm(i)),BLH(2),BLH(1),BLH(3),obs(i,4),obs(i,5),chs(i)
	   if(kobs==4)then
             write(10,'(a8,2F12.6,2F11.3,8F12.4)')trim(znm(i/2+1)),BLH(2),BLH(1),BLH(3),obs(i,4),
     *      obs(i,5),obs(i+1,5),chs(i),chs(i+1)
         endif
      enddo
      close(10)
      rst=0.d0;rst2=0.d0
 	do ni=1,nlat
        BLH(1)=hd(3)+(real(ni)-0.5d0)*hd(6)
        do nj=1,nlon
	    BLH(2)=hd(1)+(real(nj)-0.5d0)*hd(5);BLH(3)=hgt(ni,nj)
          call BLH_RLAT(GRS,BLH,rln);rr=rln(1)
          if(fknd==0.or.fknd==2)GMr=GM/rr/rr
          if(fknd==3)GMr=GM/rr/rr/rr
          if(fknd==1.or.fknd==4)then
             call GNormalfd(BLH,NFD,GRS);gr=NFD(2)
             if(fknd==1)GMr=GM/rr/gr
             if(fknd==4)GMr=GM/rr/rr/gr
          endif
          call RtGridij(rln,ki,kj,blat,lvl,nn,mm,nln,dl,lon)!计算rln在格网中的序号ki>0,kj>0
          mk=0;BB=0.d0;B2=0.d0;node=0!s-RBF系数,node(Kt)观测量作用半径范围内有效节点序号
          do i=ki-nd,ki+nd!增加实现节点序号与格网行列号对应，即节点序号二维整数数组，0表示
            if(i<1.or.i>nn)goto 3001!lon(i,j),blat第一行平行圈格网中心地心纬度°
            rlnk(2)=blat+(i-1.d0)*dlat
            do j=kj-nd,kj+nd
              if(j<1.or.j>nln(i)) goto 3002 
              if(nrd(i,j)<1) goto 3002 !nrd(i,j)未知数序号
              rlnk(3)=lon(i,j);call drln(rln,rlnk,dln)!距离m与夹角°
              if(dln(2)>dr)goto 3002
              call RBFvalue(RBF(:,fknd+1),NF,dr,dln(2),tmp)!由球面角距°dln(2)内插RBF值val
              mk=mk+1;node(mk)=nrd(i,j)
              BB(nrd(i,j))=GMr*tmp*dexp((minN-1.0)*dlog(r0/rr))
              if(fknd==4)then
                tt=dcos(dln(2)*RAD);sinf=dsqrt(1.d0-tt**2);if(sinf<1.d-12)sinf=1.d-12
	          cosa=(dcos(rln(2)*RAD)*dsin(rlnk(2)*RAD)-dsin(rln(2)*RAD)*dcos(rlnk(2)*RAD)*dcos(rlnk(3)*RAD-rln(3)*RAD))/sinf
	          sina=dcos(rlnk(2)*RAD)*dsin(rlnk(3)*RAD-rln(3)*RAD)/sinf
                BB(nrd(i,j))= BB(nrd(i,j))*cosa
                B2(nrd(i,j))= BB(nrd(i,j))*sina
              endif
3002          continue
            enddo
3001        continue    
          enddo
          val=0.d0;tmp=0.d0
          do i=1,mk
            ki=node(i)
            val=val+BB(ki)*xx(ki)
            if(fknd==4)tmp=tmp+B2(ki)*xx(ki)
          enddo
          rst(ni,nj)=val*unit(fknd+1)
          if(fknd==4)rst2(ni,nj)=tmp*unit(fknd+1)
        enddo
      enddo
      !计算并输出检核点处目标场元估值
        open(unit=8,file=checkpointfl,status="old",iostat=status)
        if(status/=0)goto 401
        open(unit=10,file="checkreslt.txt",status="replace")
        read(8,'(a)') line
        write(10,'(a)')trim(AdjustL(line))
        do while(.not.eof(8))
          read(8,'(a)') line
          call PickReclong(line,kln,rec,sn)
          if(sn<4)goto 406
          BLH(1)=rec(3);BLH(2)=rec(2);BLH(3)=rec(4)
          call BLH_RLAT(GRS,BLH,rln);rr=rln(1)
          if(fknd==0.or.fknd==2)GMr=GM/rr/rr
          if(fknd==3)GMr=GM/rr/rr/rr
          if(fknd==1.or.fknd==4)then
             call GNormalfd(BLH,NFD,GRS);gr=NFD(2)
             if(fknd==1)GMr=GM/rr/gr
             if(fknd==4)GMr=GM/rr/rr/gr
          endif
          call RtGridij(rln,ki,kj,blat,lvl,nn,mm,nln,dl,lon)!计算rln在格网中的序号ki>0,kj>0
          mk=0;BB=0.d0;B2=0.d0;node=0!s-RBF系数,node(Kt)观测量作用半径范围内有效节点序号
          do i=ki-nd,ki+nd!增加实现节点序号与格网行列号对应，即节点序号二维整数数组，0表示
            if(i<1.or.i>nn)goto 4001!lon(i,j),blat第一行平行圈格网中心地心纬度°
            rlnk(2)=blat+(i-1.d0)*dlat
            do j=kj-nd,kj+nd
              if(j<1.or.j>nln(i)) goto 4002 
              if(nrd(i,j)<1) goto 4002 !nrd(i,j)未知数序号
              rlnk(3)=lon(i,j);call drln(rln,rlnk,dln)!距离m与夹角°
              if(dln(2)>dr)goto 4002
              call RBFvalue(RBF(:,fknd+1),NF,dr,dln(2),tmp)!由球面角距°dln(2)内插RBF值val
              mk=mk+1;node(mk)=nrd(i,j)
              BB(nrd(i,j))=GMr*tmp*dexp((minN-1.0)*dlog(r0/rr))
              if(fknd==4)then
                tt=dcos(dln(2)*RAD);sinf=dsqrt(1.d0-tt**2);if(sinf<1.d-12)sinf=1.d-12
	          cosa=(dcos(rln(2)*RAD)*dsin(rlnk(2)*RAD)-dsin(rln(2)*RAD)*dcos(rlnk(2)*RAD)*dcos(rlnk(3)*RAD-rln(3)*RAD))/sinf
	          sina=dcos(rlnk(2)*RAD)*dsin(rlnk(3)*RAD-rln(3)*RAD)/sinf
                BB(nrd(i,j))= BB(nrd(i,j))*cosa
                B2(nrd(i,j))= BB(nrd(i,j))*sina
              endif
4002          continue
            enddo
4001        continue    
          enddo
          val=0.d0;tmp=0.d0
          do i=1,mk
            ki=node(i)
            val=val+BB(ki)*xx(ki)
            if(fknd==4)tmp=tmp+B2(ki)*xx(ki)
          enddo
          prst(1)=val*unit(fknd+1)
          if(fknd==4)then
             prst(2)=tmp*unit(fknd+1)
	       write(10,'(a,4F12.4)')trim(line),prst(1),prst(2)
          else
	       write(10,'(a,4F12.4)')trim(line),prst(1)
          endif
406       continue
        enddo
        close(8);close(10)
401     continue
      open(unit=10,file="SRBFestimate.dat",status="replace")
	write(10,'(a)')line0
      do i=1,nlat
	  write(10,'(15F12.4)')(rst(i,j),j=1,nlon)
      enddo
      if(fknd==3)then
        do i=1,nlat
	    write(10,'(15F12.4)')(rst2(i,j),j=1,nlon)
        enddo
      endif
      close(10)
      call SRBF4one(RBF4,RBFn,order,krbf,mpn4,mdp4,minN,maxN,NF4,nta)
      open(unit=10,file="SRBFspc.txt",status="replace")
	write(10,'(2I3,2I6,F8.2,a)')krbf,order,minN,maxN,dpth*1.d-3,
     *  " gravity disturbance, height anomaly, gravity gradient, vertical deflection" 
      tmp=dr/dble(NF4)*RAD*ae*1.d-3
 	do i=1,NF4
	   write(10,'(F12.3,8F13.5)')-(NF4-i+1.0)*tmp,(RBF4(NF4-i+2,j),j=1,4)
      enddo
 	do i=1,NF4+1
	   write(10,'(F12.3,8F13.5)')(i-1.0)*tmp,(RBF4(i,j),j=1,4)
      enddo
      close(10)
      open(unit=10,file="SRBFdgr.txt",status="replace")
	write(10,'(2I3,2I6,F8.2,a)')krbf,order,minN,maxN,dpth*1.d-3,
     * " gravity disturbance, height anomaly, gravity gradient, vertical deflection" 
 	do i=minN,maxN
	   write(10,'(I8,8F13.5)')i,(RBFn(i-minN+1,j),j=1,4)
      enddo
      close(10)
      open(unit=10,file="SRBFcenter.txt",status="replace")
	write(10,'(4I6,F8.3)')lvl,Kt,nn,nln(1),dlat*60.d0
      rln(1)=rr;k=0
      do i=1,nn!增加实现节点序号与格网行列号对应，即节点序号二维整数数组
          rln(2)=blat+(i-0.5d0)*dlat
          do j=1,nln(i)
            if(gpnt(i,j)>0)then
              k=k+1;rln(3)=lon(i,j);call RLAT_BLH(GRS,rln,BLH)
	        write(10,'(I6,2F12.6,2F8.3)')k,BLH(2),BLH(1),sr(i),dl(i)*60.d0
            endif
          enddo
      enddo
      close(10)
904   deallocate(BPB)
903   deallocate(nln,sr,dl,nrd,lon,node,enode,RBF,gpnt,RBF4,RBFn,rlatlon)
      deallocate(mpn,mdp,mpn4,mdp4,obs,BB,B2,BPL,xx,chs,chd)
901   deallocate(hgt,rst,rst2)
902   continue
101   format(a,40F12.4)
      end
!
!******************************************************************
!
      subroutine drln(rln1,rln2,dln)
      !由两点球坐标计算距离与夹角°
      implicit none
      real*8::rln1(3),rln2(3),dln(2),XYZ1(3),XYZ2(3),L2
      integer::i,j,n,m,k
      real*8::onei,Bn,pi,RAD,CnmCalc
!---------------------------------------------------------------
      RAD=datan(1.d0)/45.d0
      call RLAT_XYZ(rln1,XYZ1)
      call RLAT_XYZ(rln2,XYZ2)
      L2=(XYZ1(1)-XYZ2(1))**2+(XYZ1(2)-XYZ2(2))**2+(XYZ1(3)-XYZ2(3))**2
      dln(1)=dsqrt(L2)
      dln(2)=dabs(dacos((rln1(1)**2+rln2(1)**2-L2)/2.d0/rln1(1)/rln2(1)))/RAD
      end
!
!******************************************************************************
!
      subroutine LegPn01(mpn,mdp,minN,maxN,NF,dr)
      !计算minN~maxN阶NF+1组勒让德函数[minN~maxN]×[0,dr]
	implicit none
	integer::minN,maxN,NF,nn,i,k,n
	real*8::dr,t,dt,pi,RAD
	real*8::mpn(maxN-minN+1,NF+1),mdp(maxN-minN+1,NF+1)
	real*8,allocatable::pn(:),dp1(:)
!---------------------------------------------------------------------------
      pi=datan(1.d0)*4.d0;RAD=pi/180.d0
      dt=dr/dble(NF)*RAD
      allocate(pn(maxN+2),dp1(maxN+2))
      do k=1,NF+1
        t=dcos(dble(k-1)*dt)
        call LegPn_dt1(pn,dp1,maxN+1,t)
        do n=minN,maxN
          i=n-minN+1
          mpn(i,k)=pn(n+1)
          mdp(i,k)=dp1(n+1)
        enddo
      enddo
      deallocate(pn,dp1)
      end
!
!******************************************************************************
!
      subroutine Stat1d(dt,nn,rst)
      implicit none
	integer::nn,i,kk
      real*8::dt(nn),rst(4),pv,std,maxt,mint
!---------------------------------------------------------------------
	pv=0.d0;std=0.d0;maxt=-9.d28;mint=9.d28
      kk=0
      do i=1,nn
        if(dt(i)>9900.d0)goto 1001
        kk=kk+1;pv=pv+dt(i)
	  if(maxt<dt(i))maxt=dt(i)
        if(mint>dt(i))mint=dt(i)
1001    continue
      enddo
      pv=pv/dble(kk)
      do i=1,nn
        if(dt(i)>9900.d0)goto 1002
        std=std+(dt(i)-pv)**2
1002    continue
      enddo
      std=dsqrt(std/dble(kk))
      rst(1)=pv;rst(2)=std;rst(3)=mint;rst(4)=maxt
      end
