!  GravityfieldapproachSRBFtest.f90 
!
!  FUNCTIONS:
!  GravityfieldapproachSRBFtest - Entry point of console application.
!
!****************************************************************************
      program GravityfieldapproachSRBFtest
      implicit none
	character*800::observationfl,surfhgtgrdfl,checkpointfl
      real*8 para(20)
!---------------------------------------------------------------------
	write(observationfl,*) "resdbm541_1800.txt" !The discrete residual anomalous field element observation file
	write(surfhgtgrdfl,*) "surfhgt.dat"!The ellipsoidal height grid file of the calculation surface
	write(checkpointfl,*) "dbmchkpnt.txt"!The checkpoint file 
      para(1)=360;para(2)=1800!SRBF最小最大阶数 Minimum and maximum degree
      para(3)=1!多级次数 the order number m
      para(4)=0!0-径向多级子核函数，1-Poisson小波核函数
      para(5)=1800!SRBF的Reuter等级K the Reuter network level K for the SRBFs
      para(6)=100.d0!SRBF中心作用距离(km) the action distance of SRBF center
      para(7)=10.d0!Bjerhammar球面相对地面的平均深度 the Bjerhammar sphere burial depth (km)
      !观测场元和待估目标场元类型 the type of the observations and type of unknown target field element
      !=0扰动重力，=1高程异常，=2空间异常，=3扰动重力梯度，=4垂线偏差
      !=0 gravity disturbance (mGal), =1 height anomaly (m), =2 gravity anomaly (mGal), =3 disturbing gravity gradient
        !(E, radial) or =4 vertical deflection (″)
      para(8)=0;para(9)=1
      !观测文件记录中观测场元列序号 column ordinal number of the observations in the file record
      para(10)=7
      !当观测场元为垂线偏差时，Para(10:11)为垂线偏差南、西方向的列序号
      !When the observation is vertical deflection, para(10:11) are the column ordinal numbers of southward and westward
        !vertical deflection.
      para(11)=9
      !观测权列序号 column ordinal number of the observation weight in the file record
      !当权值属性列序号小于1，或超出记录列序号，或文件记录中权值属性小于零时，程序默认等权。
      !When the column ordinal number of the weight attribute is less than 1, exceeds the column number of the record, or
        ! the weight is less than zero, the program makes the weight equal to 1.
      !当文件记录中权值等于零时，该观测量不参与SRBF系数估计，程序可测定该观测量的外部精度指标。
      !When the weight in the file record is equal to zero, the observation will not participate in the estimation of the
        !SRBF coefficient, and the program can be employed to measure the external accuracy index of the observations.
      para(12)=0
      !选择法方程解算方法 Select the method of the solution of normal equation
      para(13)=1!1-LU分解,2-Cholesky分解,3-最小二乘QR分解,4-最小范数奇异值分解,5-岭估计
      write(*, *)"    Begin compulation, please wait......"
      call SRBFestimation(observationfl,surfhgtgrdfl,checkpointfl,para)
      write (*,*)'    Complete the computation! The approached target field elements are saved'
      write (*,*)'      in the grid file SRBFestimate.dat.'
      write (*,*)'    The program outputs also the residual observation file residuals.txt and '
      write (*,*)'      approached target field element file checkreslt.txt at checkpoints into'
      write (*,*)'      the current directory.'
      pause
      end
