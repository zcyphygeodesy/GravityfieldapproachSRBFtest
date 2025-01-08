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
      para(1)=360;para(2)=1800!SRBF��С������ Minimum and maximum degree
      para(3)=1!�༶���� the order number m
      para(4)=0!0-����༶�Ӻ˺�����1-PoissonС���˺���
      para(5)=1800!SRBF��Reuter�ȼ�K the Reuter network level K for the SRBFs
      para(6)=100.d0!SRBF�������þ���(km) the action distance of SRBF center
      para(7)=10.d0!Bjerhammar������Ե����ƽ����� the Bjerhammar sphere burial depth (km)
      !�۲ⳡԪ�ʹ���Ŀ�곡Ԫ���� the type of the observations and type of unknown target field element
      !=0�Ŷ�������=1�߳��쳣��=2�ռ��쳣��=3�Ŷ������ݶȣ�=4����ƫ��
      !=0 gravity disturbance (mGal), =1 height anomaly (m), =2 gravity anomaly (mGal), =3 disturbing gravity gradient
        !(E, radial) or =4 vertical deflection (��)
      para(8)=0;para(9)=1
      !�۲��ļ���¼�й۲ⳡԪ����� column ordinal number of the observations in the file record
      para(10)=7
      !���۲ⳡԪΪ����ƫ��ʱ��Para(10:11)Ϊ����ƫ���ϡ�������������
      !When the observation is vertical deflection, para(10:11) are the column ordinal numbers of southward and westward
        !vertical deflection.
      para(11)=9
      !�۲�Ȩ����� column ordinal number of the observation weight in the file record
      !��Ȩֵ���������С��1���򳬳���¼����ţ����ļ���¼��Ȩֵ����С����ʱ������Ĭ�ϵ�Ȩ��
      !When the column ordinal number of the weight attribute is less than 1, exceeds the column number of the record, or
        ! the weight is less than zero, the program makes the weight equal to 1.
      !���ļ���¼��Ȩֵ������ʱ���ù۲���������SRBFϵ�����ƣ�����ɲⶨ�ù۲������ⲿ����ָ�ꡣ
      !When the weight in the file record is equal to zero, the observation will not participate in the estimation of the
        !SRBF coefficient, and the program can be employed to measure the external accuracy index of the observations.
      para(12)=0
      !ѡ�񷨷��̽��㷽�� Select the method of the solution of normal equation
      para(13)=1!1-LU�ֽ�,2-Cholesky�ֽ�,3-��С����QR�ֽ�,4-��С��������ֵ�ֽ�,5-�����
      write(*, *)"    Begin compulation, please wait......"
      call SRBFestimation(observationfl,surfhgtgrdfl,checkpointfl,para)
      write (*,*)'    Complete the computation! The approached target field elements are saved'
      write (*,*)'      in the grid file SRBFestimate.dat.'
      write (*,*)'    The program outputs also the residual observation file residuals.txt and '
      write (*,*)'      approached target field element file checkreslt.txt at checkpoints into'
      write (*,*)'      the current directory.'
      pause
      end
