program wham
  implicit none
  integer, parameter :: pr=kind(1.D0)
  integer, parameter :: M=200   !** M : Number of bins
  integer :: S,i,j    !** S: Number of windows 
!  real(pr), dimension(M) :: x
  real(pr), dimension(:), allocatable :: k, x0
  real(pr), dimension(:),allocatable :: x,y,z
  real(pr), dimension(:),allocatable :: x_1neig,y_1neig,z_1neig
  real(pr), dimension(:,:),allocatable :: Coord1neig
  integer , dimension(:),allocatable :: traj,AcMove,atomNo,neig
!  real(pr), dimension(10000000):: x,y,z
!  integer , dimension(10000000):: traj,AcMove,atomNo
  integer :: Nt,io,No_ions,CumHisto,lineno,lol,mol1,mol2,dump,counter
  real(pr) :: rij,rcut,vec(3),box
  integer,dimension(0:20) :: Histogram
  integer,dimension(0:900) :: HistClustSize
  integer,dimension(0:900) :: HistClustSize1
  integer,dimension(0:900) :: HistClustSize2
  integer,dimension(0:900) :: HistClustSize3
  integer,dimension(0:900) :: HistClustSize4
  integer,dimension(0:900) :: HistClustSize5
  integer,dimension(0:900) :: HistClustSize6
  integer,dimension(0:900) :: HistClustSize7
  integer,dimension(0:900) :: HistClustSize8
  logical, dimension(:,:), allocatable :: ConnectionMatrix
  logical, dimension(:,:), allocatable :: ConnectionMatrix1, ConnectionMatrix2
  logical, dimension(:,:), allocatable :: ConnectionMatrix3, ConnectionMatrix4
  logical, dimension(:,:), allocatable :: ConnectionMatrix5, ConnectionMatrix6
  logical, dimension(:,:), allocatable :: ConnectionMatrix7, ConnectionMatrix8
  logical, dimension(:,:), allocatable :: ConnectionMatrix_LQ4crystal, ConnectionMatrix_q8crystal
  integer :: neigClustSize1,ClustSizeHead1
  integer :: neigClustSize2,ClustSizeHead2
  integer :: neigClustSize3,ClustSizeHead3
  integer :: neigClustSize4,ClustSizeHead4
  integer :: neigClustSize5,ClustSizeHead5
  integer :: neigClustSize6,ClustSizeHead6
  integer :: neigClustSize7,ClustSizeHead7
  integer :: neigClustSize8,ClustSizeHead8
  integer Largest_Cluster_Size
  integer Largest_Cluster_Head 
  integer :: count_1neig


            Nt = 2000000
            rcut=4.0
            box=52.1
            box=24.864736439521266*2
            Histogram=0
            HistClustSize1= 0
            HistClustSize2= 0
            HistClustSize3= 0
            HistClustSize4= 0
            HistClustSize5= 0
            HistClustSize6= 0
            HistClustSize7= 0
            HistClustSize8= 0

            allocate(x(Nt))
            allocate(y(Nt))
            allocate(z(Nt))
         allocate(traj(Nt))
       allocate(AcMove(Nt))
       allocate(atomNo(Nt))
               
            Nt=0
            i=0

            OPEN (1, file = 'lol')
            DO
              READ(1,*,iostat=io)
                IF (io/=0) EXIT
                Nt=Nt+1
            END DO
             CLOSE (1)

            OPEN (1, file = 'lol')
            DO
              READ(1,*,iostat=io)
                IF (io/=0) EXIT
                i=i+1
              if(i<2000)then
               read(1,*)traj(i),AcMove(i),atomNo(i),x(i),y(i),z(i)
               if(atomNo(i)<atomNo(i-1))then
                       No_ions=atomNo(i-1)
                       exit
               endif
              endif
            END DO
             CLOSE (1)
             print *,Nt,No_ions

           deallocate(x)
           deallocate(y)
           deallocate(z)
           deallocate(traj)
           deallocate(AcMove)
           deallocate(atomNo)


            allocate(x(No_ions))
            allocate(y(No_ions))
            allocate(z(No_ions))
         allocate(traj(No_ions))
       allocate(AcMove(No_ions))
       allocate(atomNo(No_ions))
       allocate(neig(No_ions))


            allocate(x_1neig(No_ions))
            allocate(y_1neig(No_ions))
            allocate(z_1neig(No_ions))
            allocate(Coord1neig(No_ions,3))


        allocate(ConnectionMatrix(No_ions,No_ions))
        allocate(ConnectionMatrix1(No_ions,No_ions))
        allocate(ConnectionMatrix2(No_ions,No_ions))
        allocate(ConnectionMatrix3(No_ions,No_ions))
        allocate(ConnectionMatrix4(No_ions,No_ions))
        allocate(ConnectionMatrix5(No_ions,No_ions))
        allocate(ConnectionMatrix6(No_ions,No_ions))
        allocate(ConnectionMatrix7(No_ions,No_ions))
        allocate(ConnectionMatrix8(No_ions,No_ions))
        allocate(ConnectionMatrix_LQ4crystal(No_ions,No_ions))
        allocate(ConnectionMatrix_q8crystal(No_ions,No_ions))

        counter=0
        lineno=0

        !Nt=No_ions*21

        open(unit=456,file='ToPlot')


        open(unit=100,file='blah')
        open(unit=101,file='ClusterCoordSize1.xyz')
        open(unit=102,file='ClusterCoordSize2.xyz')
        open(unit=103,file='ClusterCoordSize3.xyz')
        open(unit=104,file='ClusterCoordSize4.xyz')
        open(unit=105,file='ClusterCoordSize5.xyz')
        open(unit=106,file='ClusterCoordSize6.xyz')
        open(unit=107,file='ClusterCoordSize7.xyz')
        open(unit=108,file='ClusterCoordSize8.xyz')

        open(unit=201,file='Size1')
        open(unit=202,file='Size2')
        open(unit=203,file='Size3')
        open(unit=204,file='Size4')
        open(unit=205,file='Size5')
        open(unit=206,file='Size6')
        open(unit=207,file='Size7')
        open(unit=208,file='Size8')

     OPEN (1, file = 'lol')
     DO lol=0,Nt
        if(lineno>=Nt)exit

        if(mod(lineno,No_ions).eq. 0)then      
                counter=counter+1
         !  open(unit=10,file='lol')
                 do i=1,No_ions
                         read(1,*)traj(i),AcMove(i),atomNo(i),x(i),y(i),z(i)
                 enddo                
         ! close(unit=10)
         neig=0
                do i=1,No_ions
                   do j=1,No_ions
                      vec(1)=x(i)-x(j)
                      vec(2)=y(i)-y(j)
                      vec(3)=z(i)-z(j)
                      vec=vec-box*anint(vec/box)
                      rij=(vec(1)**2+vec(2)**2+vec(3)**2)**0.5
                        if(rij<rcut .and. i<=No_ions/2 .and. j>No_ions/2)then
                           neig(i)=neig(i)+1
                        elseif(rij<rcut .and. i>No_ions/2 .and. j<=No_ions/2)then
                           neig(i)=neig(i)+1
                        endif
                   enddo
                enddo                

            !** Determine connectivity among neighboring particles**
            ConnectionMatrix=.false.
            do mol1=1,No_ions-1
              ConnectionMatrix(mol1,mol1)=.true.
              do mol2=mol1+1,No_ions
                vec(1)=x(mol1)-x(mol2)
                vec(2)=y(mol1)-y(mol2)
                vec(3)=z(mol1)-z(mol2)
                vec=vec-box*anint(vec/box)
                rij=(vec(1)**2+vec(2)**2+vec(3)**2)**0.5
                 if(rij<rcut)then
                  ConnectionMatrix(mol1,mol2)=.true.
                  ConnectionMatrix(mol2,mol1)=.true.
                 end if
              end do
            end do
            ConnectionMatrix(No_ions,No_ions)=.true.
            ConnectionMatrix1=ConnectionMatrix
            ConnectionMatrix2=ConnectionMatrix
            ConnectionMatrix3=ConnectionMatrix
            ConnectionMatrix4=ConnectionMatrix
            ConnectionMatrix5=ConnectionMatrix
            ConnectionMatrix6=ConnectionMatrix
            ConnectionMatrix7=ConnectionMatrix
            ConnectionMatrix8=ConnectionMatrix
            ConnectionMatrix_LQ4crystal=ConnectionMatrix
            ConnectionMatrix_q8crystal=ConnectionMatrix

            do mol1=1,No_ions
                if(neig(mol1) <= 1)then
                ConnectionMatrix1(mol1,:)=.false.
                ConnectionMatrix1(:,mol1)=.false.
                end if
                if(neig(mol1) <= 2)then
                ConnectionMatrix2(mol1,:)=.false.
                ConnectionMatrix2(:,mol1)=.false.
                end if
                if(neig(mol1) <= 3)then
                ConnectionMatrix3(mol1,:)=.false.
                ConnectionMatrix3(:,mol1)=.false.
                end if
                if(neig(mol1) <= 4)then
                ConnectionMatrix4(mol1,:)=.false.
                ConnectionMatrix4(:,mol1)=.false.
                end if
                if(neig(mol1) <= 5)then
                ConnectionMatrix5(mol1,:)=.false.
                ConnectionMatrix5(:,mol1)=.false.
                end if
                if(neig(mol1) <= 6)then
                ConnectionMatrix6(mol1,:)=.false.
                ConnectionMatrix6(:,mol1)=.false.
                end if
                if(neig(mol1) <= 7)then
                ConnectionMatrix7(mol1,:)=.false.
                ConnectionMatrix7(:,mol1)=.false.
                end if
                if(neig(mol1) <= 8)then
                ConnectionMatrix8(mol1,:)=.false.
                ConnectionMatrix8(:,mol1)=.false.
                end if
            end do

              neigClustSize1 = Largest_Cluster_Size(ConnectionMatrix1,No_ions)
              neigClustSize2 = Largest_Cluster_Size(ConnectionMatrix2,No_ions)
              neigClustSize3 = Largest_Cluster_Size(ConnectionMatrix3,No_ions)
              neigClustSize4 = Largest_Cluster_Size(ConnectionMatrix4,No_ions)
              neigClustSize5 = Largest_Cluster_Size(ConnectionMatrix5,No_ions)
              neigClustSize6 = Largest_Cluster_Size(ConnectionMatrix6,No_ions)
              neigClustSize7 = Largest_Cluster_Size(ConnectionMatrix7,No_ions)
              neigClustSize8 = Largest_Cluster_Size(ConnectionMatrix8,No_ions)
!print *,real(lineno/No_ions,pr),neigClustSize1,neigClustSize4,neigClustSize5,neigClustSize6,neigClustSize7,neigClustSize8

              ClustSizeHead1 = Largest_Cluster_Head(ConnectionMatrix1,No_ions)
              ClustSizeHead2 = Largest_Cluster_Head(ConnectionMatrix2,No_ions)
              ClustSizeHead3 = Largest_Cluster_Head(ConnectionMatrix3,No_ions)
              ClustSizeHead4 = Largest_Cluster_Head(ConnectionMatrix4,No_ions)
              ClustSizeHead5 = Largest_Cluster_Head(ConnectionMatrix5,No_ions)
              ClustSizeHead6 = Largest_Cluster_Head(ConnectionMatrix6,No_ions)
              ClustSizeHead7 = Largest_Cluster_Head(ConnectionMatrix7,No_ions)
              ClustSizeHead8 = Largest_Cluster_Head(ConnectionMatrix8,No_ions)

              write(101,*)'Size',NeigClustSize1
              write(102,*)'Size',NeigClustSize2
              write(103,*)'Size',NeigClustSize3
              write(104,*)'Size',NeigClustSize4
              write(105,*)'Size',NeigClustSize5
              write(106,*)'Size',NeigClustSize6
              write(107,*)'Size',NeigClustSize7
              write(108,*)'Size',NeigClustSize8

              write(101,*)'timestep',counter,box
              write(102,*)'timestep',counter,box
              write(103,*)'timestep',counter,box
              write(104,*)'timestep',counter,box
              write(105,*)'timestep',counter,box
              write(106,*)'timestep',counter,box
              write(107,*)'timestep',counter,box
              write(108,*)'timestep',counter,box


            deallocate(x_1neig)
            deallocate(y_1neig)
            deallocate(z_1neig)
            deallocate(Coord1neig)

            allocate(x_1neig(No_ions))
            allocate(y_1neig(No_ions))
            allocate(z_1neig(No_ions))
            allocate(Coord1neig(No_ions,3))
            count_1neig=1
              do mol1=1,No_ions
                      if(connectionmatrix1(ClustSizeHead1,mol1))then
                               if(mol1 <= real(No_ions/2,pr))  write(101,*)'Na',x(mol1),y(mol1),z(mol1)
                               if(mol1 >  real(No_ions/2,pr))  write(101,*)'Cl',x(mol1),y(mol1),z(mol1)
                                x_1neig(count_1neig)=x(mol1)
                                y_1neig(count_1neig)=y(mol1)
                                z_1neig(count_1neig)=z(mol1)
                                Coord1neig(count_1neig,1)=x(mol1)
                                Coord1neig(count_1neig,2)=y(mol1)
                                Coord1neig(count_1neig,3)=z(mol1)
                                count_1neig=count_1neig+1
                      endif
                      if(connectionmatrix2(ClustSizeHead2,mol1))then
                               if(mol1 <= real(No_ions/2,pr))  write(102,*)'Na',x(mol1),y(mol1),z(mol1)
                               if(mol1 >  real(No_ions/2,pr))  write(102,*)'Cl',x(mol1),y(mol1),z(mol1)
                      endif
                      if(connectionmatrix3(ClustSizeHead3,mol1))then
                               if(mol1 <= real(No_ions/2,pr))  write(103,*)'Na',x(mol1),y(mol1),z(mol1)
                               if(mol1 >  real(No_ions/2,pr))  write(103,*)'Cl',x(mol1),y(mol1),z(mol1)
                      endif
                      if(connectionmatrix4(ClustSizeHead4,mol1))then
                               if(mol1 <= real(No_ions/2,pr))  write(104,*)'Na',x(mol1),y(mol1),z(mol1)
                               if(mol1 >  real(No_ions/2,pr))  write(104,*)'Cl',x(mol1),y(mol1),z(mol1)
                      endif
                      if(connectionmatrix5(ClustSizeHead5,mol1))then
                               if(mol1 <= real(No_ions/2,pr))  write(105,*)'Na',x(mol1),y(mol1),z(mol1)
                               if(mol1 >  real(No_ions/2,pr))  write(105,*)'Cl',x(mol1),y(mol1),z(mol1)
                      endif
                      if(connectionmatrix6(ClustSizeHead6,mol1))then
                               if(mol1 <= real(No_ions/2,pr))  write(106,*)'Na',x(mol1),y(mol1),z(mol1)
                               if(mol1 >  real(No_ions/2,pr))  write(106,*)'Cl',x(mol1),y(mol1),z(mol1)
                      endif
                      if(connectionmatrix7(ClustSizeHead7,mol1))then
                               if(mol1 <= real(No_ions/2,pr))  write(107,*)'Na',x(mol1),y(mol1),z(mol1)
                               if(mol1 >  real(No_ions/2,pr))  write(107,*)'Cl',x(mol1),y(mol1),z(mol1)
                      endif
                      if(connectionmatrix8(ClustSizeHead7,mol1))then
                               if(mol1 <= real(No_ions/2,pr))  write(108,*)'Na',x(mol1),y(mol1),z(mol1)
                               if(mol1 >  real(No_ions/2,pr))  write(108,*)'Cl',x(mol1),y(mol1),z(mol1)
                      endif
              end do

!print *,real(lineno/No_ions,pr),neigClustSize1,neigClustSize2,neigClustSize3,neigClustSize4,neigClustSize5,neigClustSize6
write (456,*)real(lineno/No_ions,pr),neigClustSize1,neigClustSize2,neigClustSize3,neigClustSize4,neigClustSize5,neigClustSize6&
,neigClustSize7,neigClustSize8


               write(201,*)real(lineno/No_ions,pr),neigClustSize1,ClustSizeHead1
               write(202,*)real(lineno/No_ions,pr),neigClustSize2,ClustSizeHead2
               write(203,*)real(lineno/No_ions,pr),neigClustSize3,ClustSizeHead3
               write(204,*)real(lineno/No_ions,pr),neigClustSize4,ClustSizeHead4
               write(205,*)real(lineno/No_ions,pr),neigClustSize5,ClustSizeHead5
               write(206,*)real(lineno/No_ions,pr),neigClustSize6,ClustSizeHead6
               write(207,*)real(lineno/No_ions,pr),neigClustSize7,ClustSizeHead7
               write(208,*)real(lineno/No_ions,pr),neigClustSize8,ClustSizeHead8

              HistClustSize1(neigClustSize1)= HistClustSize1(neigClustSize1)+1
              HistClustSize2(neigClustSize2)= HistClustSize2(neigClustSize2)+1
              HistClustSize3(neigClustSize3)= HistClustSize3(neigClustSize3)+1
              HistClustSize4(neigClustSize4)= HistClustSize4(neigClustSize4)+1
              HistClustSize5(neigClustSize5)= HistClustSize5(neigClustSize5)+1
              HistClustSize6(neigClustSize6)= HistClustSize6(neigClustSize6)+1
              HistClustSize7(neigClustSize7)= HistClustSize7(neigClustSize7)+1
              HistClustSize8(neigClustSize8)= HistClustSize8(neigClustSize8)+1

                do i=1,No_ions
                        Histogram(neig(i))=Histogram(neig(i))+1
                enddo


                write(100,*)No_ions
                write(100,*)'timestep',counter
                do i=1,neigClustSize1
                      if(i <= real(No_ions/2,pr))  write(100,*)'Na',x(i),y(i),z(i)
                      if(i > real(No_ions/2,pr))  write(100,*)'Cl',x(i),y(i),z(i)
                      !write(100,*)x_1neig(i),Coord1neig(i,1)
                enddo                




        endif

     lineno=lineno+1
if(lineno>=Nt)exit
    enddo        

close(unit=1)
close(unit=100)
close(unit=456)

                open(unit=550,file='Histogram_rcut4p5')
                        do i=0,20
                                write(550,*)i,Histogram(i),sum(Histogram)-sum(Histogram(0:i))
                        enddo
                close(unit=550)

               open(unit=551,file='ClustSize1_Hist_rcut4p5')
               open(unit=552,file='ClustSize2_Hist_rcut4p5')
               open(unit=553,file='ClustSize3_Hist_rcut4p5')
               open(unit=554,file='ClustSize4_Hist_rcut4p5')
               open(unit=555,file='ClustSize5_Hist_rcut4p5')
               open(unit=556,file='ClustSize6_Hist_rcut4p5')
               open(unit=557,file='ClustSize7_Hist_rcut4p5')
               open(unit=558,file='ClustSize8_Hist_rcut4p5')
                       do i=0,200
                        !      write(551,*)i,real(HistClustSize1(i),pr)/real(1000,pr),sum(HistClustSize1)-sum(HistClustSize1(0:i))
                        !      write(552,*)i,real(HistClustSize2(i),pr)/real(1000,pr),sum(HistClustSize2)-sum(HistClustSize2(0:i))
                        !      write(553,*)i,real(HistClustSize3(i),pr)/real(1000,pr),sum(HistClustSize3)-sum(HistClustSize3(0:i))
                        !      write(554,*)i,real(HistClustSize4(i),pr)/real(1000,pr),sum(HistClustSize4)-sum(HistClustSize4(0:i))
                        !      write(555,*)i,real(HistClustSize5(i),pr)/real(1000,pr),sum(HistClustSize5)-sum(HistClustSize5(0:i))
                        !      write(556,*)i,real(HistClustSize6(i),pr)/real(1000,pr),sum(HistClustSize6)-sum(HistClustSize6(0:i))
    write(551,*)i,HistClustSize1(i),sum(HistClustSize1)-sum(HistClustSize1(0:i)),real(HistClustSize1(i),pr)/real(lineno/No_ions,pr)
    write(552,*)i,HistClustSize2(i),sum(HistClustSize2)-sum(HistClustSize2(0:i)),real(HistClustSize2(i),pr)/real(lineno/No_ions,pr)
    write(553,*)i,HistClustSize3(i),sum(HistClustSize3)-sum(HistClustSize3(0:i)),real(HistClustSize3(i),pr)/real(lineno/No_ions,pr)
    write(554,*)i,HistClustSize4(i),sum(HistClustSize4)-sum(HistClustSize4(0:i)),real(HistClustSize4(i),pr)/real(lineno/No_ions,pr)
    write(555,*)i,HistClustSize5(i),sum(HistClustSize5)-sum(HistClustSize5(0:i)),real(HistClustSize5(i),pr)/real(lineno/No_ions,pr)
    write(556,*)i,HistClustSize6(i),sum(HistClustSize6)-sum(HistClustSize6(0:i)),real(HistClustSize6(i),pr)/real(lineno/No_ions,pr)
    write(557,*)i,HistClustSize7(i),sum(HistClustSize7)-sum(HistClustSize7(0:i)),real(HistClustSize7(i),pr)/real(lineno/No_ions,pr)
    write(558,*)i,HistClustSize8(i),sum(HistClustSize8)-sum(HistClustSize8(0:i)),real(HistClustSize8(i),pr)/real(lineno/No_ions,pr)
                       enddo
               close(unit=551)
               close(unit=552)
               close(unit=553)
               close(unit=554)
               close(unit=555)
               close(unit=556)
               close(unit=557)
               close(unit=558)


end program wham





  integer function Largest_Cluster_Size(ConnMat,No_ions)

    integer,intent(in) :: No_ions
    logical,dimension(No_ions,No_ions),intent(inout) :: ConnMat
    logical, dimension(No_ions) :: Row
    integer, parameter :: pr=kind(1.D0)
          integer :: mol1, mol2, icount
          integer :: NumberOfClusters, ClusterSize
          integer, dimension(No_ions) :: ClusterHead
          real(PR) :: rij,rijsq,vec(3),TempVec(3)
          integer :: LargestClusterHeadID,counter
!       
          do mol1=1,No_ions
            Row(mol1)=ConnMat(mol1,mol1)
          end do

!       
!       
        Largest_Cluster_Size=0
        LargestClusterHeadID=1
!          !** identify the clusters
        NumberOfClusters=0
        ClusterHead=0
        do mol1=1,No_ions-1
          if (Row(mol1)) then
            NumberOfClusters=NumberOfClusters+1
            ClusterHead(NumberOfClusters)=mol1
            Row(mol1)=.false.
            !**keep combining rows until all a(mol1,mol2)=1 have solid(mol2)=false**
           do
             icount=0        
             do mol2=mol1+1,No_ions
               if (ConnMat(mol1,mol2) .and. Row(mol2)) then
                 icount=icount+1
                 Row(mol2)=.false.
                 !**determine union of the two rows**
                 ConnMat(mol1,:)=ConnMat(mol1,:) .or. ConnMat(mol2,:)
               endif
             enddo
             if(icount == 0)exit
          end do
     
            ClusterSize=0
            do mol2=1,No_ions
              if(ConnMat(mol1,mol2))ClusterSize=ClusterSize+1
            end do
           if(ClusterSize>Largest_Cluster_Size)LargestClusterHeadID=mol1
           Largest_Cluster_Size=max(ClusterSize,Largest_Cluster_Size)
          endif
        enddo
        if (Row(No_ions)) then
          NumberOfClusters=NumberOfClusters+1
           ClusterSize=1
           ClusterHead(NumberOfClusters)=No_ions
           if(ClusterSize>Largest_Cluster_Size)LargestClusterHeadID=No_ions
           Largest_Cluster_Size=max(ClusterSize,Largest_Cluster_Size)
         endif
  !print *,Largest_Cluster_Size,LargestClusterHeadID     

  end function Largest_Cluster_Size




  integer function Largest_Cluster_Head(ConnMat,No_ions)

    integer,intent(in) :: No_ions
    logical,dimension(No_ions,No_ions),intent(inout) :: ConnMat
    logical, dimension(No_ions) :: Row
    integer, parameter :: pr=kind(1.D0)
          integer :: mol1, mol2, icount
          integer :: NumberOfClusters, ClusterSize
          integer, dimension(No_ions) :: ClusterHead
          real(PR) :: rij,rijsq,vec(3),TempVec(3)
          integer :: LargestClusterHeadID,counter
!       
          do mol1=1,No_ions
            Row(mol1)=ConnMat(mol1,mol1)
          end do

!       
!       
        Largest_Cluster_Size=0
        LargestClusterHeadID=1
!          !** identify the clusters
        NumberOfClusters=0
        ClusterHead=0
        do mol1=1,No_ions-1
          if (Row(mol1)) then
            NumberOfClusters=NumberOfClusters+1
            ClusterHead(NumberOfClusters)=mol1
            Row(mol1)=.false.
            !**keep combining rows until all a(mol1,mol2)=1 have solid(mol2)=false**
           do
             icount=0        
             do mol2=mol1+1,No_ions
               if (ConnMat(mol1,mol2) .and. Row(mol2)) then
                 icount=icount+1
                 Row(mol2)=.false.
                 !**determine union of the two rows**
                 ConnMat(mol1,:)=ConnMat(mol1,:) .or. ConnMat(mol2,:)
               endif
             enddo
             if(icount == 0)exit
          end do
     
            ClusterSize=0
            do mol2=1,No_ions
              if(ConnMat(mol1,mol2))ClusterSize=ClusterSize+1
            end do
           if(ClusterSize>Largest_Cluster_Size)LargestClusterHeadID=mol1
           Largest_Cluster_Size=max(ClusterSize,Largest_Cluster_Size)
          endif
        enddo
        if (Row(No_ions)) then
          NumberOfClusters=NumberOfClusters+1
           ClusterSize=1
           ClusterHead(NumberOfClusters)=No_ions
           if(ClusterSize>Largest_Cluster_Size)LargestClusterHeadID=No_ions
           Largest_Cluster_Size=max(ClusterSize,Largest_Cluster_Size)
         endif
       
         Largest_Cluster_Head=LargestClusterHeadID
  end function Largest_Cluster_Head

        !====================================================================================================
!      
!              subroutine Compute_q8_CrystalSize(IonCenterOfMass,NumberOfIons)
!      
!                integer, parameter :: pr=kind(1.D0)
!                complex(PR), dimension(:,:), allocatable :: q8m
!                real(PR), dimension(:), allocatable :: q8
!               ! logical, dimension(:), allocatable :: solid
!                real(pr), dimension(:,:),allocatable :: IonCenterOfMass
!                integer :: mol1, mol2
!                integer :: l,m,i,j
!                real(PR), dimension(3) :: vec
!                complex(PR) :: ylm
!                integer :: NumberOfNeighbors,NumberOfIons
!                real(PR), dimension(:), allocatable :: rijsq,theta,phi
!                allocate(IonCenterOfMass(NumberOfIons,3))
!                allocate(q8m(-8:8,NumberOfIons))
!                allocate(q8(NumberOfIons))
!             
!                NumberOfNeighbors=12
!                allocate(rijsq(NumberOfIons-1),theta(NumberOfIons-1),phi(NumberOfIons-1))
!                q8m=(0._PR, 0._PR)
!                do mol1=1,NumberOfIons
!                  !** Calculate distances
!                  i=1
!                  do mol2=1,NumberOfIons
!                    if(mol1 /= mol2)then
!                      vec=IonCenterOfMass(:,mol1)-IonCenterOfMass(:,mol2)
!                      vec=vec-box*anint(vec/box)
!                      rijsq(i)=vec(1)**2+vec(2)**2+vec(3)**2
!                      theta(i)=atan2(sqrt(vec(1)**2+vec(2)**2),vec(3))
!                      phi(i)=atan2(vec(2),vec(1))
!                      i=i+1
!                    end if
!                  end do
!             
!                  !** Sort arrays
!                  i=2
!                  do while(i<NumberOfIons)
!                    j=i
!                    do !while(j > 1 .and. (rijsq(j-1) > rijsq(j)))
!                      if(rijsq(j-1) > rijsq(j))then
!                        call swap(rijsq(j-1),rijsq(j))
!                        call swap(theta(j-1),theta(j))
!                        call swap(phi(j-1),phi(j))
!                      end if
!                      j=j-1
!                      if(j == 1)exit
!                    end do
!                    i=i+1
!                  end do
!             
!                  !** Compute Spherical harmonics for 12 nearest neighbors
!                  do i=1,NumberOfNeighbors
!                    !** q8
!                    l=8
!                    do m=-l,l
!                      ylm=sp_hrmcs(l,m,theta(i),phi(i))
!                      q8m(m,mol1)=q8m(m,mol1)+ylm
!                    enddo
!                  end do
!                  q8m(:,mol1)=q8m(:,mol1)/NumberOfNeighbors
!                  q8(mol1)=sqrt(4._PR*PI/real(2*l+1,PR)*real(dot_product(q8m(:,mol1),q8m(:,mol1))))
!                end do 
!             
!                deallocate(rijsq,theta,phi)
!             
!             
!              end subroutine Compute_q8_CrystalSize
!      
