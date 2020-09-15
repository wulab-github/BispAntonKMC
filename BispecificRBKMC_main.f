      program RBKMC_TwoSpeciesPlm3D3D
      implicit none

      integer simu_step
      parameter (simu_step=100000000)
      real*8 time_step
      parameter (time_step=0.000000001)! Unit: Second
      real*8 distance_step,distance_amp
      parameter (distance_step=1.0)! Unit: nm
      real*8 cell_range_x
      parameter (cell_range_x=1000.0)! Unit: nm
      real*8 cell_range_y
      parameter (cell_range_y=1000.0)! Unit: nm
      real*8 cell_range_z
      parameter (cell_range_z=250.0)! Unit: nm
c>>   the z-coordinate of 2D membrane surface is 0, the extracellular space is from 0 to cell_range_z
c>>   Periodic boundary is applied in x and y axis, but not z, the z coordinates of 3D molecules are not allowed to be smaller than 0 and larger than cell_range_z 
      real*8 TCR_radius
      parameter (TCR_radius=4.0)! Unit: nm
      real*8 TCR_D
      parameter (TCR_D=10000000.0) ! Unit: nm square per second
      real*8 TCR_rot_D
      parameter (TCR_rot_D=0.25)! unit: degree
      real*8 CD_radius
      parameter (CD_radius=4.0)! Unit: nm
      real*8 CD_D
      parameter (CD_D=10000000.0)! Unit: nm square per second
      real*8 CD_rot_D
      parameter (CD_rot_D=0.25)! unit: degree
      real*8 MHC_radius
      parameter (MHC_radius=4.0)! Unit: nm
      real*8 B7_radius
      parameter (B7_radius=4.0)! Unit: nm
      real*8 Plm_radius
      parameter (Plm_radius=6.6)! Unit: nm  ! the average interdomain distance of the ligand
CC>>>>>>> Linker specific
      real*8 Plm_D
      parameter (Plm_D=67000000.0)! Unit: nm square per second
      real*8 Plm_rot_D
      parameter (Plm_rot_D=0.25)! unit: degree
      real*8 Plm_Flx_TlRg
c      parameter (Plm_Flx_TlRg=2.0)! Unit: nm
      real*8 Plm_Flx_RtRg
c      parameter (Plm_Flx_RtRg=10.0)! unit: degree
      real*8 complex_D
      parameter (complex_D=5000000.0)! Unit: nm square per second
      real*8 Complex_rot_D
      parameter (Complex_rot_D=0.1)! unit: degree
      real*8 EffectiveVolume
      parameter (EffectiveVolume=1.0) ! Unit is molecule per nm cubic, the DissCalibration will be estimated by molecular simulation by evaluating the entropy effect during dimerization.
      real*8 DissCalibration
      parameter (DissCalibration=EffectiveVolume/4.0) ! unit is per Mole, the unit of 1.6 is Mole nm cubic per molecule
      real*8 AssCalibration
      parameter (AssCalibration=1.0)
      real*8 Ass_Rate_TM
      parameter (Ass_Rate_TM=35000000.0) ! unit is molecule per second ! WT = 3.5E7; MT1 = 8E6;
      real*8 BindingAffinity_TM
      parameter (BindingAffinity_TM=-13.0) ! WT = -13RT; MT1 = -8RT;
      real*8 Diss_Rate_TM
      parameter (Diss_Rate_TM=Ass_Rate_TM*DissCalibration*1
     &     *DEXP(BindingAffinity_TM))
      real*8 Ass_Rate_CB
      parameter (Ass_Rate_CB=30000000.0)  ! WT = 3E7; MT2 = 1.5E7; MT3 = 4.5E7
      real*8 BindingAffinity_CB
      parameter (BindingAffinity_CB=-12.0) ! WT = -12RT; MT2 = -10RT; MT3 = -14RT; MT4 is double mutant
      real*8 Diss_Rate_CB
      parameter (Diss_Rate_CB=Ass_Rate_CB*DissCalibration*1
     &     *DEXP(BindingAffinity_CB))
      integer TCR_tot_num
      parameter (TCR_tot_num=300) ! TgCel = 300; NtCel = 30
      integer TCR_res_num
      parameter (TCR_res_num=4)
      integer CD_tot_num
      parameter (CD_tot_num=300)
      integer CD_res_num
      parameter (CD_res_num=4)
      integer Plm_tot_num
      parameter (Plm_tot_num=300)
      integer Plm_res_num
      parameter (Plm_res_num=4)
      integer Plm_MHC_num
      parameter (Plm_MHC_num=1)
      integer Plm_B7_num
      parameter (Plm_B7_num=1)
      integer num_trajec
      parameter (num_trajec=1)
      real*8 pai
      parameter (pai=3.1415926) 
      real*8 bond_dist_cutoff_TM
      parameter (bond_dist_cutoff_TM=2.0)
      real*8 bond_thetapd_TM,bond_thetapd_cutoff_TM
      parameter (bond_thetapd_TM=180.0,bond_thetapd_cutoff_TM=30.0)
      real*8 bond_thetaot_TM,bond_thetaot_cutoff_TM
      parameter (bond_thetaot_TM=90.0,bond_thetaot_cutoff_TM=30.0)
      real*8 bond_dist_cutoff_CB
      parameter (bond_dist_cutoff_CB=1.5)
      real*8 bond_thetapd_CB,bond_thetapd_cutoff_CB
      parameter (bond_thetapd_CB=180.0,bond_thetapd_cutoff_CB=30.0)
      real*8 bond_thetaot_CB,bond_thetaot_cutoff_CB
      parameter (bond_thetaot_CB=90.0,bond_thetaot_cutoff_CB=30.0)
      integer ene_output_flag,trj_output_flag
      parameter (ene_output_flag=1,trj_output_flag=1)
      integer nbin
      parameter (nbin=10)

cccccccccccccccccccccccccccccccccccccccc

      real*8 TCR_x(TCR_tot_num,TCR_res_num)
      real*8 TCR_y(TCR_tot_num,TCR_res_num)
      real*8 TCR_z(TCR_tot_num,TCR_res_num)
      real*8 TCR_x_0(TCR_tot_num,TCR_res_num)
      real*8 TCR_y_0(TCR_tot_num,TCR_res_num)
      real*8 TCR_z_0(TCR_tot_num,TCR_res_num)
      integer TCR_status(TCR_tot_num)
      real*8 TCR_x_new(TCR_tot_num,TCR_res_num)
      real*8 TCR_y_new(TCR_tot_num,TCR_res_num)
      real*8 TCR_z_new(TCR_tot_num,TCR_res_num)
      real*8 TCR_x_new0(TCR_tot_num,TCR_res_num)
      real*8 TCR_y_new0(TCR_tot_num,TCR_res_num)
      real*8 TCR_z_new0(TCR_tot_num,TCR_res_num)
      integer TCR_status_new(TCR_tot_num)
      real*8 CD_x(CD_tot_num,CD_res_num)
      real*8 CD_y(CD_tot_num,CD_res_num)
      real*8 CD_z(CD_tot_num,CD_res_num)
      real*8 CD_x_0(CD_tot_num,CD_res_num)
      real*8 CD_y_0(CD_tot_num,CD_res_num)
      real*8 CD_z_0(CD_tot_num,CD_res_num)
      integer CD_status(CD_tot_num)
      real*8 CD_x_new(CD_tot_num,CD_res_num)
      real*8 CD_y_new(CD_tot_num,CD_res_num)
      real*8 CD_z_new(CD_tot_num,CD_res_num)
      real*8 CD_x_new0(CD_tot_num,CD_res_num)
      real*8 CD_y_new0(CD_tot_num,CD_res_num)
      real*8 CD_z_new0(CD_tot_num,CD_res_num)
      integer CD_status_new(CD_tot_num)
      real*8 Plm_MHC_x(Plm_tot_num,Plm_MHC_num,Plm_res_num)
      real*8 Plm_MHC_y(Plm_tot_num,Plm_MHC_num,Plm_res_num)
      real*8 Plm_MHC_z(Plm_tot_num,Plm_MHC_num,Plm_res_num)
      real*8 Plm_MHC_xo(Plm_tot_num,Plm_MHC_num,Plm_res_num)
      real*8 Plm_MHC_yo(Plm_tot_num,Plm_MHC_num,Plm_res_num)
      real*8 Plm_MHC_zo(Plm_tot_num,Plm_MHC_num,Plm_res_num)
      real*8 Plm_MHC_x_0(Plm_tot_num,Plm_MHC_num,Plm_res_num)
      real*8 Plm_MHC_y_0(Plm_tot_num,Plm_MHC_num,Plm_res_num)
      real*8 Plm_MHC_z_0(Plm_tot_num,Plm_MHC_num,Plm_res_num)
      real*8 Plm_MHC_xo_0(Plm_tot_num,Plm_MHC_num,Plm_res_num)
      real*8 Plm_MHC_yo_0(Plm_tot_num,Plm_MHC_num,Plm_res_num)
      real*8 Plm_MHC_zo_0(Plm_tot_num,Plm_MHC_num,Plm_res_num)
      real*8 Plm_MHC_x_new(Plm_tot_num,Plm_MHC_num,Plm_res_num)
      real*8 Plm_MHC_y_new(Plm_tot_num,Plm_MHC_num,Plm_res_num)
      real*8 Plm_MHC_z_new(Plm_tot_num,Plm_MHC_num,Plm_res_num)
      real*8 Plm_MHC_xo_new(Plm_tot_num,Plm_MHC_num,Plm_res_num)
      real*8 Plm_MHC_yo_new(Plm_tot_num,Plm_MHC_num,Plm_res_num)
      real*8 Plm_MHC_zo_new(Plm_tot_num,Plm_MHC_num,Plm_res_num)
      real*8 Plm_MHC_x_new0(Plm_tot_num,Plm_MHC_num,Plm_res_num)
      real*8 Plm_MHC_y_new0(Plm_tot_num,Plm_MHC_num,Plm_res_num)
      real*8 Plm_MHC_z_new0(Plm_tot_num,Plm_MHC_num,Plm_res_num)
      real*8 Plm_MHC_xo_new0(Plm_tot_num,Plm_MHC_num,Plm_res_num)
      real*8 Plm_MHC_yo_new0(Plm_tot_num,Plm_MHC_num,Plm_res_num)
      real*8 Plm_MHC_zo_new0(Plm_tot_num,Plm_MHC_num,Plm_res_num)
      integer Plm_MHC_status(Plm_tot_num,Plm_MHC_num)
      integer Plm_MHC_status_new(Plm_tot_num,Plm_MHC_num)
      integer Plm_MHC_idx(Plm_tot_num,Plm_MHC_num)
      integer Plm_MHC_idx_new(Plm_tot_num,Plm_MHC_num)
      real*8 Plm_B7_x(Plm_tot_num,Plm_B7_num,Plm_res_num)
      real*8 Plm_B7_y(Plm_tot_num,Plm_B7_num,Plm_res_num)
      real*8 Plm_B7_z(Plm_tot_num,Plm_B7_num,Plm_res_num)
      real*8 Plm_B7_xo(Plm_tot_num,Plm_B7_num,Plm_res_num)
      real*8 Plm_B7_yo(Plm_tot_num,Plm_B7_num,Plm_res_num)
      real*8 Plm_B7_zo(Plm_tot_num,Plm_B7_num,Plm_res_num)
      real*8 Plm_B7_x_0(Plm_tot_num,Plm_B7_num,Plm_res_num)
      real*8 Plm_B7_y_0(Plm_tot_num,Plm_B7_num,Plm_res_num)
      real*8 Plm_B7_z_0(Plm_tot_num,Plm_B7_num,Plm_res_num)
      real*8 Plm_B7_xo_0(Plm_tot_num,Plm_B7_num,Plm_res_num)
      real*8 Plm_B7_yo_0(Plm_tot_num,Plm_B7_num,Plm_res_num)
      real*8 Plm_B7_zo_0(Plm_tot_num,Plm_B7_num,Plm_res_num)
      real*8 Plm_B7_x_new(Plm_tot_num,Plm_B7_num,Plm_res_num)
      real*8 Plm_B7_y_new(Plm_tot_num,Plm_B7_num,Plm_res_num)
      real*8 Plm_B7_z_new(Plm_tot_num,Plm_B7_num,Plm_res_num)
      real*8 Plm_B7_xo_new(Plm_tot_num,Plm_B7_num,Plm_res_num)
      real*8 Plm_B7_yo_new(Plm_tot_num,Plm_B7_num,Plm_res_num)
      real*8 Plm_B7_zo_new(Plm_tot_num,Plm_B7_num,Plm_res_num)
      real*8 Plm_B7_x_new0(Plm_tot_num,Plm_B7_num,Plm_res_num)
      real*8 Plm_B7_y_new0(Plm_tot_num,Plm_B7_num,Plm_res_num)
      real*8 Plm_B7_z_new0(Plm_tot_num,Plm_B7_num,Plm_res_num)
      real*8 Plm_B7_xo_new0(Plm_tot_num,Plm_B7_num,Plm_res_num)
      real*8 Plm_B7_yo_new0(Plm_tot_num,Plm_B7_num,Plm_res_num)
      real*8 Plm_B7_zo_new0(Plm_tot_num,Plm_B7_num,Plm_res_num)
      integer Plm_B7_status(Plm_tot_num,Plm_B7_num)
      integer Plm_B7_status_new(Plm_tot_num,Plm_B7_num)
      integer Plm_B7_idx(Plm_tot_num,Plm_B7_num)
      integer Plm_B7_idx_new(Plm_tot_num,Plm_B7_num)

      integer i,j,k,n_t,j2
      real*8 temp_i,temp_j,temp_k
      real*8 theta,phi,psai,phai
      real*8 t(3,3)
      real*8 dist
      real*8 cm1_a_x
      real*8 cm1_a_y
      real*8 cm1_a_z
      real*8 cm1_b_x
      real*8 cm1_b_y
      real*8 cm1_b_z
      real*8 PB_x,PB_y,PB_z
      real*8 initial_simu_time,current_simu_time
      integer mc_time_step
      integer iteration_mole_step
      integer selecting_mole_index
      integer TCR_index,CD_index
      integer Plm_index
      real*8 Prob_Diff
      integer collision_flag
      real*8 point_x(3),point_y(3),point_z(3)
      real*8 theta_pd,theta_ot
      real*8 Prob_Ass,Prob_Diss
      integer selected_Plm,selected_TCR,selected_CD
      real*8 rot_angle
      real*8 x_axis_bg,y_axis_bg,z_axis_bg
      real*8 x_axis_ed,y_axis_ed,z_axis_ed
      real*8 x_bf_rot,y_bf_rot,z_bf_rot
      real*8 x_af_rot,y_af_rot,z_af_rot
      integer upper_bound_flag,lower_bound_flag
      real*8 Flx_Rg_x,Flx_Rg_y,Flx_Rg_z,Flx_Rg
      integer Plm_bound_flag
      real*8 temp_PB_x,temp_PB_y
      integer TCR_complex_num,CD_complex_num
      real*8 temp,prob
      integer part1,part2,part3
      integer part4,part5,part6
      integer part7,part8,part9
      real*8 Tlprob_dist_dom1(nbin),Tlprob_dist_dom2(nbin)
      real*8 Tlcutoff_dist_dom1(nbin),Tlcutoff_dist_dom2(nbin)
      real*8 Tlculm_dist_dom1(nbin),Tlculm_dist_dom2(nbin)
      real*8 Rtprob_dist_dom1(nbin),Rtprob_dist_dom2(nbin)
      real*8 Rtcutoff_dist_dom1(nbin),Rtcutoff_dist_dom2(nbin)
      real*8 Rtculm_dist_dom1(nbin),Rtculm_dist_dom2(nbin)
      real*8 amplitude_binselect1,amplitude_binselect2
      integer amplitude_binindex1,amplitude_binindex2
      integer DoubleBind_num,MHC_flag,B7_flag

      character*25 serialnum

      real rand3
      double precision r3
      real rand4
      double precision r4
      real rand5
      double precision r5

      r4=5.0
      r5=5.0    
      r3=5.0      

cccccccccccccccccccccccccccccccccccccccc
c   record simulation parameter
cccccccccccccccccccccccccccccccccccccccc      

      serialnum='SamplTry_TgCel_GS15_PDWT0'

ccccccccccccccccccccccccccccccccccccccccccccccc

      open(unit=10,file=
     &     'BispRBKMCresult_simupara_'
     &     //serialnum//'.txt',
     &     status='unknown')
      write(10,*) 'cell_range_x,cell_range_y,cell_range_z',
     &     cell_range_x,cell_range_y,cell_range_z
      write(10,*) 'time_step',time_step
      write(10,*) 'distance_step',distance_step
      write(10,*) 'TCR_radius',TCR_radius,'TCR_D',TCR_D
      write(10,*) 'TCR_rot_D',TCR_rot_D
      write(10,*) 'CD_radius',CD_radius,'CD_D',CD_D
      write(10,*) 'CD_rot_D',CD_rot_D
      write(10,*) 'MHC_radius',MHC_radius
      write(10,*) 'B7_radius',B7_radius
      write(10,*) 'Plm_radius,Plm_D',Plm_radius,Plm_D
      write(10,*) 'Plm_rot_D',Plm_rot_D
      write(10,*) 'complex_D',complex_D
      write(10,*) 'Complex_rot_D',Complex_rot_D
      write(10,*) 'Plm_Flx_TlRg',Plm_Flx_TlRg
      write(10,*) 'Plm_Flx_RtRg',Plm_Flx_RtRg
      write(10,*) 'Ass_Rate_TM',Ass_Rate_TM
      write(10,*) 'Ass_Rate_CB',Ass_Rate_CB
      write(10,*) 'BindingAffinity_TM',BindingAffinity_TM
      write(10,*) 'BindingAffinity_CB',BindingAffinity_CB
      write(10,*) 'TCR_tot_num',TCR_tot_num
      write(10,*) 'CD_tot_num',CD_tot_num
      write(10,*) 'Plm_tot_num',Plm_tot_num
      write(10,*) 'Plm_MHC_num',Plm_MHC_num
      write(10,*) 'Plm_B7_num',Plm_B7_num
      write(10,*) 'bond_dist_cutoff_TM',bond_dist_cutoff_TM
      write(10,*) 'bond_thetapd_TM,bond_thetapd_cutoff_TM',
     &     bond_thetapd_TM,bond_thetapd_cutoff_TM
      write(10,*) 'bond_thetaot_TM,bond_thetaot_cutoff_TM',
     &     bond_thetaot_TM,bond_thetaot_cutoff_TM
      write(10,*) 'bond_dist_cutoff_CB',bond_dist_cutoff_CB
      write(10,*) 'bond_thetapd_CB,bond_thetapd_cutoff_CB',
     &     bond_thetapd_CB,bond_thetapd_cutoff_CB
      write(10,*) 'bond_thetaot_CB,bond_thetaot_cutoff_CB',
     &     bond_thetaot_CB,bond_thetaot_cutoff_CB

      close(10)
                       
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c>>   read fluctuation probability
ccccccccccccccccccccccccccccccccccccccccccccccccccccc

      open(unit=10,file='GS15_TlTotFl_ProbDistribution.dat',
     &     status='old')
      
      do i=1,nbin
         read(10,1000) Tlprob_dist_dom1(i),Tlculm_dist_dom1(i),
     &        Tlcutoff_dist_dom1(i),
     &        Tlprob_dist_dom2(i),Tlculm_dist_dom2(i),
     &        Tlcutoff_dist_dom2(i)
      enddo
      close(10)

      open(unit=10,file='GS15_RtTotFl_ProbDistribution.dat',
     &     status='old')
      
      do i=1,nbin
         read(10,1000) Rtprob_dist_dom1(i),Rtculm_dist_dom1(i),
     &        Rtcutoff_dist_dom1(i),
     &        Rtprob_dist_dom2(i),Rtculm_dist_dom2(i),
     &        Rtcutoff_dist_dom2(i)
      enddo
      close(10)

 1000 format(5x,2F8.5,1x,F10.3,1x,2F8.5,1x,F10.3)


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc   generate multiple trajectories
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do n_t=1,num_trajec


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c   construct the initialized position and conformation of molecules in 3D and 2D
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   construct the position and conformation of TCR on 2D

c>>>>>   random position

         do i=1,TCR_tot_num
 100        continue
            temp_i=rand3(r3)*cell_range_x-cell_range_x/2  
            temp_j=rand3(r3)*cell_range_y-cell_range_y/2             
            do j=1,i-1
               dist=sqrt((temp_i-TCR_x(j,1))**2+
     &              (temp_j-TCR_y(j,1))**2)
               if(dist.le.
     &              (TCR_radius+TCR_radius))then
                  goto 100
               endif
            enddo
            TCR_x(i,1)=temp_i
            TCR_y(i,1)=temp_j
            TCR_z(i,1)=0
            TCR_x_0(i,2)=temp_i+TCR_radius
            TCR_y_0(i,2)=temp_j
            TCR_z_0(i,2)=0
            TCR_x_0(i,3)=temp_i
            TCR_y_0(i,3)=temp_j+TCR_radius
            TCR_z_0(i,3)=0
            TCR_x(i,4)=temp_i
            TCR_y(i,4)=temp_j
            TCR_z(i,4)=TCR_radius
            TCR_status(i)=0

c>>>>>   random orientation along membrane normal

            rot_angle=2*rand3(r3)*180-180
            x_axis_bg=TCR_x(i,1)
            y_axis_bg=TCR_y(i,1)
            z_axis_bg=-100.0
            x_axis_ed=TCR_x(i,TCR_res_num)
            y_axis_ed=TCR_y(i,TCR_res_num)
            z_axis_ed=100.0
            do j=2,TCR_res_num-1
               x_bf_rot=TCR_x_0(i,j)
               y_bf_rot=TCR_y_0(i,j)
               z_bf_rot=TCR_z_0(i,j)
               x_af_rot=0
               y_af_rot=0
               z_af_rot=0
               call rot_along_axis(rot_angle,
     &              x_axis_bg,y_axis_bg,z_axis_bg,
     &              x_axis_ed,y_axis_ed,z_axis_ed,
     &              x_bf_rot,y_bf_rot,z_bf_rot,
     &              x_af_rot,y_af_rot,z_af_rot)
               TCR_x(i,j)=x_af_rot
               TCR_y(i,j)=y_af_rot
               TCR_z(i,j)=z_af_rot
            enddo

         enddo

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   construct the position and conformation of CD28 on 2D

c>>>>>   random position

         do i=1,CD_tot_num
 101        continue
            temp_i=rand3(r3)*cell_range_x-cell_range_x/2  
            temp_j=rand3(r3)*cell_range_y-cell_range_y/2             
            do j=1,i-1
               dist=sqrt((temp_i-CD_x(j,1))**2+
     &              (temp_j-CD_y(j,1))**2)
               if(dist.le.
     &              (CD_radius+CD_radius))then
                  goto 101
               endif
            enddo
            do j=1,TCR_tot_num
               dist=sqrt((temp_i-TCR_x(j,1))**2+
     &              (temp_j-TCR_y(j,1))**2)
               if(dist.le.
     &              (TCR_radius+CD_radius))then
                  goto 101
               endif
            enddo
            CD_x(i,1)=temp_i
            CD_y(i,1)=temp_j
            CD_z(i,1)=0
            CD_x_0(i,2)=temp_i+CD_radius
            CD_y_0(i,2)=temp_j
            CD_z_0(i,2)=0
            CD_x_0(i,3)=temp_i
            CD_y_0(i,3)=temp_j+CD_radius
            CD_z_0(i,3)=0
            CD_x(i,4)=temp_i
            CD_y(i,4)=temp_j
            CD_z(i,4)=CD_radius
            CD_status(i)=0

c>>>>>   random orientation along membrane normal

            rot_angle=2*rand3(r3)*180-180
            x_axis_bg=CD_x(i,1)
            y_axis_bg=CD_y(i,1)
            z_axis_bg=-100.0
            x_axis_ed=CD_x(i,TCR_res_num)
            y_axis_ed=CD_y(i,TCR_res_num)
            z_axis_ed=100.0
            do j=2,CD_res_num-1
               x_bf_rot=CD_x_0(i,j)
               y_bf_rot=CD_y_0(i,j)
               z_bf_rot=CD_z_0(i,j)
               x_af_rot=0
               y_af_rot=0
               z_af_rot=0
               call rot_along_axis(rot_angle,
     &              x_axis_bg,y_axis_bg,z_axis_bg,
     &              x_axis_ed,y_axis_ed,z_axis_ed,
     &              x_bf_rot,y_bf_rot,z_bf_rot,
     &              x_af_rot,y_af_rot,z_af_rot)
               CD_x(i,j)=x_af_rot
               CD_y(i,j)=y_af_rot
               CD_z(i,j)=z_af_rot
            enddo

         enddo

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   construct the position and conformation of MHC-B7 Polymer

c>>>>>   random position

         do i=1,Plm_tot_num
 200        continue
            temp_i=rand3(r3)*cell_range_x-cell_range_x/2  
            temp_j=rand3(r3)*cell_range_y-cell_range_y/2  
            temp_k=rand3(r3)*cell_range_z           
            do j=1,i-1
               cm1_a_x=0
               cm1_a_y=0
               cm1_a_z=0
               do k=1,Plm_MHC_num
                  cm1_a_x=cm1_a_x+Plm_MHC_xo(j,k,1)
                  cm1_a_y=cm1_a_y+Plm_MHC_yo(j,k,1)
                  cm1_a_z=cm1_a_z+Plm_MHC_zo(j,k,1)
               enddo
               cm1_a_x=cm1_a_x/real(Plm_MHC_num)
               cm1_a_y=cm1_a_y/real(Plm_MHC_num)
               cm1_a_z=cm1_a_z/real(Plm_MHC_num)
               dist=sqrt((temp_i-cm1_a_x)**2+
     &              (temp_j-cm1_a_y)**2+
     &              (temp_k-cm1_a_z)**2)
               if(dist.le.
     &              (2*(Plm_radius+1.0*(MHC_radius+B7_radius))))then
                  goto 200
               endif
            enddo
            do j=1,TCR_tot_num
               dist=sqrt((temp_i-TCR_x(j,1))**2+
     &              (temp_j-TCR_y(j,1))**2+(temp_k-TCR_z(j,1))**2)
               if(dist.le.
     &              (TCR_radius+(Plm_radius+0.5*(MHC_radius+B7_radius)))
     &              )then
                  goto 200
               endif
            enddo
            do j=1,CD_tot_num
               dist=sqrt((temp_i-CD_x(j,1))**2+
     &              (temp_j-CD_y(j,1))**2+(temp_k-CD_z(j,1))**2)
               if(dist.le.
     &              (CD_radius+(Plm_radius+0.5*(MHC_radius+B7_radius)))
     &              )then
                  goto 200
               endif
            enddo
            Plm_MHC_xo_0(i,1,1)=temp_i+MHC_radius+0.5*Plm_radius
            Plm_MHC_yo_0(i,1,1)=temp_j
            Plm_MHC_zo_0(i,1,1)=temp_k
            Plm_MHC_xo_0(i,1,2)=temp_i+MHC_radius+0.5*Plm_radius
     &           +MHC_radius
            Plm_MHC_yo_0(i,1,2)=temp_j
            Plm_MHC_zo_0(i,1,2)=temp_k
            Plm_MHC_xo_0(i,1,3)=temp_i+MHC_radius+0.5*Plm_radius
            Plm_MHC_yo_0(i,1,3)=temp_j
     &           +MHC_radius
            Plm_MHC_zo_0(i,1,3)=temp_k
            Plm_MHC_xo_0(i,1,4)=temp_i+MHC_radius+0.5*Plm_radius
            Plm_MHC_yo_0(i,1,4)=temp_j
            Plm_MHC_zo_0(i,1,4)=temp_k+MHC_radius

            Plm_B7_xo_0(i,1,1)=temp_i-B7_radius-0.5*Plm_radius
            Plm_B7_yo_0(i,1,1)=temp_j
            Plm_B7_zo_0(i,1,1)=temp_k
            Plm_B7_xo_0(i,1,2)=temp_i-B7_radius-0.5*Plm_radius
     &           -B7_radius
            Plm_B7_yo_0(i,1,2)=temp_j-B7_radius
            Plm_B7_zo_0(i,1,2)=temp_k
            Plm_B7_xo_0(i,1,3)=temp_i-B7_radius-0.5*Plm_radius
            Plm_B7_yo_0(i,1,3)=temp_j
     &           -B7_radius
            Plm_B7_zo_0(i,1,3)=temp_k
            Plm_B7_xo_0(i,1,4)=temp_i-B7_radius-0.5*Plm_radius
            Plm_B7_yo_0(i,1,4)=temp_j
            Plm_B7_zo_0(i,1,4)=temp_k+B7_radius

            
c>>>>>   random orientation

            cm1_a_x=0
            cm1_a_y=0
            cm1_a_z=0
            do j=1,Plm_MHC_num
               do k=1,Plm_res_num
                  cm1_a_x=cm1_a_x+Plm_MHC_xo_0(i,j,k)
                  cm1_a_y=cm1_a_y+Plm_MHC_yo_0(i,j,k)
                  cm1_a_z=cm1_a_z+Plm_MHC_zo_0(i,j,k)
               enddo
            enddo
            do j=1,Plm_B7_num
               do k=1,Plm_res_num
                  cm1_a_x=cm1_a_x+Plm_B7_xo_0(i,j,k)
                  cm1_a_y=cm1_a_y+Plm_B7_yo_0(i,j,k)
                  cm1_a_z=cm1_a_z+Plm_B7_zo_0(i,j,k)
               enddo
            enddo
            cm1_a_x=cm1_a_x/real(Plm_MHC_num*Plm_res_num
     &           +Plm_B7_num*Plm_res_num)
            cm1_a_y=cm1_a_y/real(Plm_MHC_num*Plm_res_num
     &           +Plm_B7_num*Plm_res_num)
            cm1_a_z=cm1_a_z/real(Plm_MHC_num*Plm_res_num
     &           +Plm_B7_num*Plm_res_num)

            theta=(2*rand3(r3)-1)*pai
            phi=(2*rand3(r3)-1)*pai
            psai=(2*rand3(r3)-1)*pai
            
            t(1,1)=cos(psai)*cos(phi)-cos(theta)*sin(phi)*sin(psai)
            t(1,2)=-sin(psai)*cos(phi)-cos(theta)*sin(phi)*cos(psai)
            t(1,3)=sin(theta)*sin(phi)
            
            t(2,1)=cos(psai)*sin(phi)+cos(theta)*cos(phi)*sin(psai)
            t(2,2)=-sin(psai)*sin(phi)+cos(theta)*cos(phi)*cos(psai)
            t(2,3)=-sin(theta)*cos(phi)
            
            t(3,1)=sin(psai)*sin(theta)
            t(3,2)=cos(psai)*sin(theta)
            t(3,3)=cos(theta)

            do j=1,Plm_MHC_num
               do k=1,Plm_res_num
                  Plm_MHC_xo(i,j,k)=t(1,1)*
     &                 (Plm_MHC_xo_0(i,j,k)-cm1_a_x)+
     &                 t(1,2)*(Plm_MHC_yo_0(i,j,k)-cm1_a_y)
     &                 +t(1,3)*(Plm_MHC_zo_0(i,j,k)-cm1_a_z)
     &                 +cm1_a_x
                  Plm_MHC_yo(i,j,k)=t(2,1)*
     &                 (Plm_MHC_xo_0(i,j,k)-cm1_a_x)+
     &                 t(2,2)*(Plm_MHC_yo_0(i,j,k)-cm1_a_y)
     &                 +t(2,3)*(Plm_MHC_zo_0(i,j,k)-cm1_a_z)
     &                 +cm1_a_y
                  Plm_MHC_zo(i,j,k)=t(3,1)*
     &                 (Plm_MHC_xo_0(i,j,k)-cm1_a_x)+
     &                 t(3,2)*(Plm_MHC_yo_0(i,j,k)-cm1_a_y)
     &                 +t(3,3)*(Plm_MHC_zo_0(i,j,k)-cm1_a_z)
     &                 +cm1_a_z
               enddo
            enddo

            do j=1,Plm_B7_num
               do k=1,Plm_res_num
                  Plm_B7_xo(i,j,k)=t(1,1)*(Plm_B7_xo_0(i,j,k)-cm1_a_x)+
     &                 t(1,2)*(Plm_B7_yo_0(i,j,k)-cm1_a_y)
     &                 +t(1,3)*(Plm_B7_zo_0(i,j,k)-cm1_a_z)
     &                 +cm1_a_x
                  Plm_B7_yo(i,j,k)=t(2,1)*(Plm_B7_xo_0(i,j,k)-cm1_a_x)+
     &                 t(2,2)*(Plm_B7_yo_0(i,j,k)-cm1_a_y)
     &                 +t(2,3)*(Plm_B7_zo_0(i,j,k)-cm1_a_z)
     &                 +cm1_a_y
                  Plm_B7_zo(i,j,k)=t(3,1)*(Plm_B7_xo_0(i,j,k)-cm1_a_x)+
     &                 t(3,2)*(Plm_B7_yo_0(i,j,k)-cm1_a_y)
     &                 +t(3,3)*(Plm_B7_zo_0(i,j,k)-cm1_a_z)
     &                 +cm1_a_z
               enddo
            enddo

c>>>>>   random translational flexibility



            do j=1,Plm_MHC_num

               amplitude_binselect1=rand5(r5)
               amplitude_binindex1=0
               do k=1,nbin
                  if(k.eq.1)then
                     if(amplitude_binselect1.lt.Tlculm_dist_dom1(k))then
                        amplitude_binindex1=k
                     endif
                  elseif(k.gt.1)then
                     if((amplitude_binselect1.lt.Tlculm_dist_dom1(k))
     &                    .and.
     &                    (amplitude_binselect1.ge.
     &                    Tlculm_dist_dom1(k-1))
     &                    )then
                        amplitude_binindex1=k
                     endif           
                  endif
               enddo       

               if(amplitude_binindex1.gt.1)then
                  Flx_Rg=(rand3(r3)*(Tlcutoff_dist_dom1
     &                 (amplitude_binindex1)-
     &                 Tlcutoff_dist_dom1
     &                 (amplitude_binindex1-1))
     &                 +Tlcutoff_dist_dom1(amplitude_binindex1-1))*0.1
               elseif(amplitude_binindex1.eq.1)then
                  Flx_Rg=rand3(r3)*Tlcutoff_dist_dom1
     &                 (amplitude_binindex1)*0.1
               endif
               theta=rand3(r3)*pai
               phai=rand3(r3)*2*pai
               do k=1,Plm_res_num
                  Plm_MHC_x_0(i,j,k)=Plm_MHC_xo(i,j,k)+
     &                 Flx_Rg*sin(theta)*cos(phai)
                  Plm_MHC_y_0(i,j,k)=Plm_MHC_yo(i,j,k)+
     &                 Flx_Rg*sin(theta)*sin(phai)
                  Plm_MHC_z_0(i,j,k)=Plm_MHC_zo(i,j,k)+
     &                 Flx_Rg*cos(theta)
               enddo
            enddo


            do j=1,Plm_B7_num

               amplitude_binselect2=rand5(r5)
               amplitude_binindex2=0
               do k=1,nbin
                  if(k.eq.1)then
                     if(amplitude_binselect2.lt.Tlculm_dist_dom2(k))then
                        amplitude_binindex2=k
                     endif
                  elseif(k.gt.1)then
                     if((amplitude_binselect2.lt.Tlculm_dist_dom2(k))
     &                    .and.
     &                    (amplitude_binselect2.ge.
     &                    Tlculm_dist_dom2(k-1))
     &                    )then
                        amplitude_binindex2=k
                     endif           
                  endif
               enddo      
 
               if(amplitude_binindex2.gt.1)then
                  Flx_Rg=(rand3(r3)*(Tlcutoff_dist_dom2
     &                 (amplitude_binindex2)-
     &                 Tlcutoff_dist_dom2
     &                 (amplitude_binindex2-1))
     &                 +Tlcutoff_dist_dom2(amplitude_binindex2-1))*0.1
               elseif(amplitude_binindex2.eq.1)then
                  Flx_Rg=rand3(r3)*Tlcutoff_dist_dom2
     &                 (amplitude_binindex2)*0.1
               endif
               theta=rand3(r3)*pai
               phai=rand3(r3)*2*pai
               do k=1,Plm_res_num
                  Plm_B7_x_0(i,j,k)=Plm_B7_xo(i,j,k)+
     &                 Flx_Rg*sin(theta)*cos(phai)
                  Plm_B7_y_0(i,j,k)=Plm_B7_yo(i,j,k)+
     &                 Flx_Rg*sin(theta)*sin(phai)
                  Plm_B7_z_0(i,j,k)=Plm_B7_zo(i,j,k)+
     &                 Flx_Rg*cos(theta)
               enddo
            enddo

c>>>>>   random rotational flexibility

            do j=1,Plm_MHC_num

               Plm_MHC_x(i,j,1)=Plm_MHC_x_0(i,j,1)
               Plm_MHC_y(i,j,1)=Plm_MHC_y_0(i,j,1)
               Plm_MHC_z(i,j,1)=Plm_MHC_z_0(i,j,1)

               amplitude_binselect1=rand5(r5)
               amplitude_binindex1=0
               do k=1,nbin
                  if(k.eq.1)then
                     if(amplitude_binselect1.lt.Rtculm_dist_dom1(k))then
                        amplitude_binindex1=k
                     endif
                  elseif(k.gt.1)then
                     if((amplitude_binselect1.lt.Rtculm_dist_dom1(k))
     &                    .and.
     &                    (amplitude_binselect1.ge.
     &                    Rtculm_dist_dom1(k-1))
     &                    )then
                        amplitude_binindex1=k
                     endif           
                  endif
               enddo       

               if(amplitude_binindex1.gt.1)then
                  Plm_Flx_RtRg=rand3(r3)*(Rtcutoff_dist_dom1
     &                 (amplitude_binindex1)-
     &                 Rtcutoff_dist_dom1
     &                 (amplitude_binindex1-1))
     &                 +Rtcutoff_dist_dom1(amplitude_binindex1-1)
               elseif(amplitude_binindex1.eq.1)then
                  Plm_Flx_RtRg=rand3(r3)*Rtcutoff_dist_dom1
     &                 (amplitude_binindex1)
               endif

               if(rand5(r5).lt.0.5)then
                  Plm_Flx_RtRg=-Plm_Flx_RtRg
               endif

               theta=Plm_Flx_RtRg*pai/180.0
               phi=Plm_Flx_RtRg*pai/180.0
               psai=Plm_Flx_RtRg*pai/180.0
               
               t(1,1)=cos(psai)*cos(phi)-cos(theta)*sin(phi)*sin(psai)
               t(1,2)=-sin(psai)*cos(phi)-cos(theta)*sin(phi)*cos(psai)
               t(1,3)=sin(theta)*sin(phi)
               
               t(2,1)=cos(psai)*sin(phi)+cos(theta)*cos(phi)*sin(psai)
               t(2,2)=-sin(psai)*sin(phi)+cos(theta)*cos(phi)*cos(psai)
               t(2,3)=-sin(theta)*cos(phi)
               
               t(3,1)=sin(psai)*sin(theta)
               t(3,2)=cos(psai)*sin(theta)
               t(3,3)=cos(theta)
               
               do k=2,Plm_res_num
                  Plm_MHC_x(i,j,k)=t(1,1)*(Plm_MHC_x_0(i,j,k)
     &                 -Plm_MHC_x(i,j,1))+
     &                 t(1,2)*(Plm_MHC_y_0(i,j,k)-Plm_MHC_y(i,j,1))
     &                 +t(1,3)*(Plm_MHC_z_0(i,j,k)-Plm_MHC_z(i,j,1))
     &                 +Plm_MHC_x(i,j,1)
                  Plm_MHC_y(i,j,k)=t(2,1)*(Plm_MHC_x_0(i,j,k)
     &                 -Plm_MHC_x(i,j,1))+
     &                 t(2,2)*(Plm_MHC_y_0(i,j,k)-Plm_MHC_y(i,j,1))
     &                 +t(2,3)*(Plm_MHC_z_0(i,j,k)-Plm_MHC_z(i,j,1))
     &                 +Plm_MHC_y(i,j,1)
                  Plm_MHC_z(i,j,k)=t(3,1)*(Plm_MHC_x_0(i,j,k)
     &                 -Plm_MHC_x(i,j,1))+
     &                 t(3,2)*(Plm_MHC_y_0(i,j,k)-Plm_MHC_y(i,j,1))
     &                 +t(3,3)*(Plm_MHC_z_0(i,j,k)-Plm_MHC_z(i,j,1))
     &                 +Plm_MHC_z(i,j,1)
               enddo

            enddo

            do j=1,Plm_B7_num

               Plm_B7_x(i,j,1)=Plm_B7_x_0(i,j,1)
               Plm_B7_y(i,j,1)=Plm_B7_y_0(i,j,1)
               Plm_B7_z(i,j,1)=Plm_B7_z_0(i,j,1)

               amplitude_binselect2=rand5(r5)
               amplitude_binindex2=0
               do k=1,nbin
                  if(k.eq.1)then
                     if(amplitude_binselect2.lt.Tlculm_dist_dom2(k))then
                        amplitude_binindex2=k
                     endif
                  elseif(k.gt.1)then
                     if((amplitude_binselect2.lt.Tlculm_dist_dom2(k))
     &                    .and.
     &                    (amplitude_binselect2.ge.
     &                    Tlculm_dist_dom2(k-1))
     &                    )then
                        amplitude_binindex2=k
                     endif           
                  endif
               enddo      

               if(amplitude_binindex2.gt.1)then
                  Plm_Flx_RtRg=rand3(r3)*(Rtcutoff_dist_dom2
     &                 (amplitude_binindex2)-
     &                 Rtcutoff_dist_dom2
     &                 (amplitude_binindex2-1))
     &                 +Rtcutoff_dist_dom2(amplitude_binindex2-1)
               elseif(amplitude_binindex2.eq.1)then
                  Plm_Flx_RtRg=rand3(r3)*Rtcutoff_dist_dom2
     &                 (amplitude_binindex2)
               endif

               if(rand5(r5).lt.0.5)then
                  Plm_Flx_RtRg=-Plm_Flx_RtRg
               endif

               theta=Plm_Flx_RtRg*pai/180.0
               phi=Plm_Flx_RtRg*pai/180.0
               psai=Plm_Flx_RtRg*pai/180.0
               
               t(1,1)=cos(psai)*cos(phi)-cos(theta)*sin(phi)*sin(psai)
               t(1,2)=-sin(psai)*cos(phi)-cos(theta)*sin(phi)*cos(psai)
               t(1,3)=sin(theta)*sin(phi)
               
               t(2,1)=cos(psai)*sin(phi)+cos(theta)*cos(phi)*sin(psai)
               t(2,2)=-sin(psai)*sin(phi)+cos(theta)*cos(phi)*cos(psai)
               t(2,3)=-sin(theta)*cos(phi)
               
               t(3,1)=sin(psai)*sin(theta)
               t(3,2)=cos(psai)*sin(theta)
               t(3,3)=cos(theta)
               
               do k=2,Plm_res_num
                  Plm_B7_x(i,j,k)=t(1,1)*(Plm_B7_x_0(i,j,k)
     &                 -Plm_B7_x(i,j,1))+
     &                 t(1,2)*(Plm_B7_y_0(i,j,k)-Plm_B7_y(i,j,1))
     &                 +t(1,3)*(Plm_B7_z_0(i,j,k)-Plm_B7_z(i,j,1))
     &                 +Plm_B7_x(i,j,1)
                  Plm_B7_y(i,j,k)=t(2,1)*(Plm_B7_x_0(i,j,k)
     &                 -Plm_B7_x(i,j,1))+
     &                 t(2,2)*(Plm_B7_y_0(i,j,k)-Plm_B7_y(i,j,1))
     &                 +t(2,3)*(Plm_B7_z_0(i,j,k)-Plm_B7_z(i,j,1))
     &                 +Plm_B7_y(i,j,1)
                  Plm_B7_z(i,j,k)=t(3,1)*(Plm_B7_x_0(i,j,k)
     &                 -Plm_B7_x(i,j,1))+
     &                 t(3,2)*(Plm_B7_y_0(i,j,k)-Plm_B7_y(i,j,1))
     &                 +t(3,3)*(Plm_B7_z_0(i,j,k)-Plm_B7_z(i,j,1))
     &                 +Plm_B7_z(i,j,1)
               enddo

            enddo


         enddo


         do i=1,Plm_tot_num
            do j=1,Plm_MHC_num
               Plm_MHC_status(i,j)=0
               Plm_MHC_idx(i,j)=0
            enddo
            do j=1,Plm_B7_num
               Plm_B7_status(i,j)=0
               Plm_B7_idx(i,j)=0
            enddo
         enddo

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      begin  main  loop of Diffusion-Reaction simulation
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         initial_simu_time=0.0

         do mc_time_step=1,simu_step

            current_simu_time=current_simu_time+time_step

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   perform the diffusion for molecules
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   

            do i=1,TCR_tot_num
               do j=1,TCR_res_num
                  TCR_x_new(i,j)=TCR_x(i,j)
                  TCR_y_new(i,j)=TCR_y(i,j)
                  TCR_z_new(i,j)=TCR_z(i,j)
               enddo
               TCR_status_new(i)=TCR_status(i)
            enddo
            do i=1,CD_tot_num
               do j=1,CD_res_num
                  CD_x_new(i,j)=CD_x(i,j)
                  CD_y_new(i,j)=CD_y(i,j)
                  CD_z_new(i,j)=CD_z(i,j)
               enddo
               CD_status_new(i)=CD_status(i)
            enddo

            do i=1,Plm_tot_num
               do j=1,Plm_MHC_num
                  do k=1,Plm_res_num
                     Plm_MHC_x_new(i,j,k)=Plm_MHC_x(i,j,k)
                     Plm_MHC_y_new(i,j,k)=Plm_MHC_y(i,j,k)
                     Plm_MHC_z_new(i,j,k)=Plm_MHC_z(i,j,k)
                     Plm_MHC_xo_new(i,j,k)=Plm_MHC_xo(i,j,k)
                     Plm_MHC_yo_new(i,j,k)=Plm_MHC_yo(i,j,k)
                     Plm_MHC_zo_new(i,j,k)=Plm_MHC_zo(i,j,k)
                  enddo
                  Plm_MHC_status_new(i,j)=Plm_MHC_status(i,j)
                  Plm_MHC_idx_new(i,j)=Plm_MHC_idx(i,j)
               enddo
               do j=1,Plm_B7_num
                  do k=1,Plm_res_num
                     Plm_B7_x_new(i,j,k)=Plm_B7_x(i,j,k)
                     Plm_B7_y_new(i,j,k)=Plm_B7_y(i,j,k)
                     Plm_B7_z_new(i,j,k)=Plm_B7_z(i,j,k)
                     Plm_B7_xo_new(i,j,k)=Plm_B7_xo(i,j,k)
                     Plm_B7_yo_new(i,j,k)=Plm_B7_yo(i,j,k)
                     Plm_B7_zo_new(i,j,k)=Plm_B7_zo(i,j,k)
                  enddo
                  Plm_B7_status_new(i,j)=Plm_B7_status(i,j)
                  Plm_B7_idx_new(i,j)=Plm_B7_idx(i,j)
               enddo
            enddo

            do iteration_mole_step=1,
     &           TCR_tot_num+CD_tot_num
     &           +Plm_tot_num

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   randomly select one molecule
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c               selecting_mole_index=int(rand3(r3)*(TCR_tot_num+
c     &              CD_tot_num+Plm_tot_num))+1

               selecting_mole_index=iteration_mole_step

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   if the selected molecule is an MHC_B7 Polymer

               if(selecting_mole_index.le.Plm_tot_num)then
c                  Plm_index=int(rand3(r3)*Plm_tot_num)+1
                  Plm_index=selecting_mole_index
cc>>>>   make random diffusion for selected molecule
cc>>>>   if it is unbound to any substrate
                  Plm_bound_flag=0
                  do i=1,Plm_MHC_num
                     if(Plm_MHC_status_new(Plm_index,i).eq.1)then
                        Plm_bound_flag=1
                     endif
                  enddo
                  do i=1,Plm_B7_num
                     if(Plm_B7_status_new(Plm_index,i).eq.1)then
                        Plm_bound_flag=1
                     endif
                  enddo
c                  if(Plm_index.eq.18)then
c                     print*,'18 bopund',Plm_bound_flag
c                  endif
                  if(Plm_bound_flag.eq.0)then
                     Prob_Diff=(6*Plm_D*time_step)/(distance_step**2)
                     if(rand3(r3).lt.Prob_Diff)then
cc>>     if smaller than the diffusion probability, move alone
                        theta=rand3(r3)*pai
                        phai=rand3(r3)*2*pai
                        distance_amp=2*rand4(r4)*distance_step
                        do i=1,Plm_MHC_num
                           do j=1,Plm_res_num
                              Plm_MHC_xo_new0(Plm_index,i,j)=
     &                             Plm_MHC_xo(Plm_index,i,j)+
     &                             distance_amp*sin(theta)*cos(phai)
                              Plm_MHC_yo_new0(Plm_index,i,j)=
     &                             Plm_MHC_yo(Plm_index,i,j)+
     &                             distance_amp*sin(theta)*sin(phai)
                              Plm_MHC_zo_new0(Plm_index,i,j)=
     &                             Plm_MHC_zo(Plm_index,i,j)+
     &                             distance_amp*cos(theta)
                           enddo
                        enddo
                        do i=1,Plm_B7_num
                           do j=1,Plm_res_num
                              Plm_B7_xo_new0(Plm_index,i,j)=
     &                             Plm_B7_xo(Plm_index,i,j)+
     &                             distance_amp*sin(theta)*cos(phai)
                              Plm_B7_yo_new0(Plm_index,i,j)=
     &                             Plm_B7_yo(Plm_index,i,j)+
     &                             distance_amp*sin(theta)*sin(phai)
                              Plm_B7_zo_new0(Plm_index,i,j)=
     &                             Plm_B7_zo(Plm_index,i,j)+
     &                             distance_amp*cos(theta)
                           enddo
                        enddo
                        lower_bound_flag=0
                        do i=1,Plm_MHC_num
                           do j=1,Plm_res_num
                              if(Plm_MHC_zo_new0(Plm_index,i,j)
     &                             .le.0)then
                                 lower_bound_flag=1
                              endif
                           enddo
                        enddo
                        do i=1,Plm_B7_num
                           do j=1,Plm_res_num
                              if(Plm_B7_zo_new0(Plm_index,i,j)
     &                             .le.0)then
                                 lower_bound_flag=1
                              endif
                           enddo
                        enddo
                        upper_bound_flag=0
                        do i=1,Plm_MHC_num
                           do j=1,Plm_res_num
                              if(Plm_MHC_zo_new0(Plm_index,i,j)
     &                             .ge.cell_range_z)then
                                 upper_bound_flag=1
                              endif
                           enddo
                        enddo
                        do i=1,Plm_B7_num
                           do j=1,Plm_res_num
                              if(Plm_B7_zo_new0(Plm_index,i,j)
     &                             .ge.cell_range_z)then
                                 upper_bound_flag=1
                              endif
                           enddo
                        enddo
                        if(lower_bound_flag.eq.1)then
                           do i=1,Plm_MHC_num
                              do j=1,Plm_res_num
                                 Plm_MHC_zo_new0(Plm_index,i,j)=
     &                                Plm_MHC_zo(Plm_index,i,j)+
     &                                abs(2*distance_step*cos(theta))
                              enddo
                           enddo
                           do i=1,Plm_B7_num
                              do j=1,Plm_res_num
                                 Plm_B7_zo_new0(Plm_index,i,j)=
     &                                Plm_B7_zo(Plm_index,i,j)+
     &                                abs(2*distance_step*cos(theta))
                              enddo
                           enddo
                        elseif(upper_bound_flag.eq.1)then
                           do i=1,Plm_MHC_num
                              do j=1,Plm_res_num
                                 Plm_MHC_zo_new0(Plm_index,i,j)=
     &                                Plm_MHC_zo(Plm_index,i,j)-
     &                                abs(2*distance_step*cos(theta))
                              enddo
                           enddo
                           do i=1,Plm_B7_num
                              do j=1,Plm_res_num
                                 Plm_B7_zo_new0(Plm_index,i,j)=
     &                                Plm_B7_zo(Plm_index,i,j)-
     &                                abs(2*distance_step*cos(theta))
                              enddo
                           enddo
                        endif
                        temp_PB_x=0
                        temp_PB_y=0
                        do i=1,Plm_MHC_num
                           temp_PB_x=temp_PB_x+
     &                          Plm_MHC_xo_new0(Plm_index,i,1)
                           temp_PB_y=temp_PB_y+
     &                          Plm_MHC_yo_new0(Plm_index,i,1)
                        enddo
                        temp_PB_x=temp_PB_x/real(Plm_MHC_num)
                        temp_PB_y=temp_PB_y/real(Plm_MHC_num)
                        PB_x=cell_range_x*anint
     &                       (temp_PB_x/cell_range_x)
                        PB_y=cell_range_y*anint
     &                       (temp_PB_y/cell_range_y)
                        do i=1,Plm_MHC_num
                           do j=1,Plm_res_num
                              Plm_MHC_xo_new0(Plm_index,i,j)=
     &                             Plm_MHC_xo_new0(Plm_index,i,j)-PB_x     
                              Plm_MHC_yo_new0(Plm_index,i,j)=
     &                             Plm_MHC_yo_new0(Plm_index,i,j)-PB_y
                           enddo
                        enddo
                        do i=1,Plm_B7_num
                           do j=1,Plm_res_num
                              Plm_B7_xo_new0(Plm_index,i,j)=
     &                             Plm_B7_xo_new0(Plm_index,i,j)-PB_x     
                              Plm_B7_yo_new0(Plm_index,i,j)=
     &                             Plm_B7_yo_new0(Plm_index,i,j)-PB_y
                           enddo
                        enddo
                     else
                        do i=1,Plm_MHC_num
                           do j=1,Plm_res_num
                              Plm_MHC_xo_new0(Plm_index,i,j)=
     &                             Plm_MHC_xo(Plm_index,i,j)
                              Plm_MHC_yo_new0(Plm_index,i,j)=
     &                             Plm_MHC_yo(Plm_index,i,j)
                              Plm_MHC_zo_new0(Plm_index,i,j)=
     &                             Plm_MHC_zo(Plm_index,i,j)
                           enddo
                        enddo
                        do i=1,Plm_B7_num
                           do j=1,Plm_res_num
                              Plm_B7_xo_new0(Plm_index,i,j)=
     &                             Plm_B7_xo(Plm_index,i,j)
                              Plm_B7_yo_new0(Plm_index,i,j)=
     &                             Plm_B7_yo(Plm_index,i,j)
                              Plm_B7_zo_new0(Plm_index,i,j)=
     &                             Plm_B7_zo(Plm_index,i,j)
                           enddo
                        enddo
                     endif

c>>>>>>>>> Rotate the molecule
                     cm1_a_x=0
                     cm1_a_y=0
                     cm1_a_z=0
                     do j=1,Plm_MHC_num
                        do k=1,Plm_res_num
                           cm1_a_x=cm1_a_x+
     &                          Plm_MHC_xo_new0(Plm_index,j,k)
                           cm1_a_y=cm1_a_y+
     &                          Plm_MHC_yo_new0(Plm_index,j,k)
                           cm1_a_z=cm1_a_z+
     &                          Plm_MHC_zo_new0(Plm_index,j,k)
                        enddo
                     enddo
                     do j=1,Plm_B7_num
                        do k=1,Plm_res_num
                           cm1_a_x=cm1_a_x+Plm_B7_xo_new0(Plm_index,j,k)
                           cm1_a_y=cm1_a_y+Plm_B7_yo_new0(Plm_index,j,k)
                           cm1_a_z=cm1_a_z+Plm_B7_zo_new0(Plm_index,j,k)
                        enddo
                     enddo
                     cm1_a_x=cm1_a_x/real(Plm_MHC_num*Plm_res_num
     &                    +Plm_B7_num*Plm_res_num)
                     cm1_a_y=cm1_a_y/real(Plm_MHC_num*Plm_res_num
     &                    +Plm_B7_num*Plm_res_num)
                     cm1_a_z=cm1_a_z/real(Plm_MHC_num*Plm_res_num
     &                    +Plm_B7_num*Plm_res_num)
                     
                     theta=(2*rand3(r3)-1)*Plm_rot_D*pai/180.0
                     phi=(2*rand3(r3)-1)*Plm_rot_D*pai/180.0
                     psai=(2*rand3(r3)-1)*Plm_rot_D*pai/180.0
                     
                     t(1,1)=cos(psai)*cos(phi)-
     &                    cos(theta)*sin(phi)*sin(psai)
                     t(1,2)=-sin(psai)*cos(phi)-
     &                    cos(theta)*sin(phi)*cos(psai)
                     t(1,3)=sin(theta)*sin(phi)
                     
                     t(2,1)=cos(psai)*sin(phi)+
     &                    cos(theta)*cos(phi)*sin(psai)
                     t(2,2)=-sin(psai)*sin(phi)+
     &                    cos(theta)*cos(phi)*cos(psai)
                     t(2,3)=-sin(theta)*cos(phi)
                     
                     t(3,1)=sin(psai)*sin(theta)
                     t(3,2)=cos(psai)*sin(theta)
                     t(3,3)=cos(theta)

                     do j=1,Plm_MHC_num
                        do k=1,Plm_res_num
                           Plm_MHC_xo_new(Plm_index,j,k)=t(1,1)*
     &                          (Plm_MHC_xo_new0(Plm_index,j,k)
     &                          -cm1_a_x)+
     &                          t(1,2)*(Plm_MHC_yo_new0(Plm_index,j,k)
     &                          -cm1_a_y)
     &                          +t(1,3)*(Plm_MHC_zo_new0(Plm_index,j,k)
     &                          -cm1_a_z)
     &                          +cm1_a_x
                           Plm_MHC_yo_new(Plm_index,j,k)=t(2,1)*
     &                          (Plm_MHC_xo_new0(Plm_index,j,k)
     &                          -cm1_a_x)+
     &                          t(2,2)*(Plm_MHC_yo_new0(Plm_index,j,k)
     &                          -cm1_a_y)
     &                          +t(2,3)*(Plm_MHC_zo_new0(Plm_index,j,k)
     &                          -cm1_a_z)
     &                          +cm1_a_y
                           Plm_MHC_zo_new(Plm_index,j,k)=t(3,1)*
     &                          (Plm_MHC_xo_new0(Plm_index,j,k)
     &                          -cm1_a_x)+
     &                          t(3,2)*(Plm_MHC_yo_new0(Plm_index,j,k)
     &                          -cm1_a_y)
     &                          +t(3,3)*(Plm_MHC_zo_new0(Plm_index,j,k)
     &                          -cm1_a_z)
     &                          +cm1_a_z
                        enddo
                     enddo

                     do j=1,Plm_B7_num
                        do k=1,Plm_res_num
                           Plm_B7_xo_new(Plm_index,j,k)=t(1,1)*
     &                          (Plm_B7_xo_new0(Plm_index,j,k)
     &                          -cm1_a_x)+
     &                          t(1,2)*(Plm_B7_yo_new0(Plm_index,j,k)
     &                          -cm1_a_y)
     &                          +t(1,3)*(Plm_B7_zo_new0(Plm_index,j,k)
     &                          -cm1_a_z)
     &                          +cm1_a_x
                           Plm_B7_yo_new(Plm_index,j,k)=t(2,1)*
     &                          (Plm_B7_xo_new0(Plm_index,j,k)
     &                          -cm1_a_x)+
     &                          t(2,2)*(Plm_B7_yo_new0(Plm_index,j,k)
     &                          -cm1_a_y)
     &                          +t(2,3)*(Plm_B7_zo_new0(Plm_index,j,k)
     &                          -cm1_a_z)
     &                          +cm1_a_y
                           Plm_B7_zo_new(Plm_index,j,k)=t(3,1)*
     &                          (Plm_B7_xo_new0(Plm_index,j,k)
     &                          -cm1_a_x)+
     &                          t(3,2)*(Plm_B7_yo_new0(Plm_index,j,k)
     &                          -cm1_a_y)
     &                          +t(3,3)*(Plm_B7_zo_new0(Plm_index,j,k)
     &                          -cm1_a_z)
     &                          +cm1_a_z
                        enddo
                     enddo

c>>>>>   random translational flexibility

                     do j=1,Plm_MHC_num

                        amplitude_binselect1=rand5(r5)
                        amplitude_binindex1=0
                        do k=1,nbin
                           if(k.eq.1)then
                              if(amplitude_binselect1.lt.
     &                             Tlculm_dist_dom1(k))then
                                 amplitude_binindex1=k
                              endif
                           elseif(k.gt.1)then
                              if((amplitude_binselect1.lt.
     &                             Tlculm_dist_dom1(k))
     &                             .and.
     &                             (amplitude_binselect1.ge.
     &                             Tlculm_dist_dom1(k-1))
     &                             )then
                                 amplitude_binindex1=k
                              endif           
                           endif
                        enddo       
                        
                        if(amplitude_binindex1.gt.1)then
                           Flx_Rg=(rand3(r3)*(Tlcutoff_dist_dom1
     &                          (amplitude_binindex1)-
     &                          Tlcutoff_dist_dom1
     &                          (amplitude_binindex1-1))
     &                          +Tlcutoff_dist_dom1
     &                          (amplitude_binindex1-1))*0.1
                        elseif(amplitude_binindex1.eq.1)then
                           Flx_Rg=rand3(r3)*Tlcutoff_dist_dom1
     &                          (amplitude_binindex1)*0.1
                        endif
                        theta=rand3(r3)*pai
                        phai=rand3(r3)*2*pai
                        do k=1,Plm_res_num
                           Plm_MHC_x_new0(Plm_index,j,k)=
     &                          Plm_MHC_xo_new(Plm_index,j,k)
     &                          +Flx_Rg*sin(theta)*cos(phai)
                           Plm_MHC_y_new0(Plm_index,j,k)=
     &                          Plm_MHC_yo_new(Plm_index,j,k)
     &                          +Flx_Rg*sin(theta)*sin(phai)
                           Plm_MHC_z_new0(Plm_index,j,k)=
     &                          Plm_MHC_zo_new(Plm_index,j,k)
     &                          +Flx_Rg*cos(theta)
                        enddo
                     enddo
                     
                     do j=1,Plm_B7_num

                        amplitude_binselect2=rand5(r5)
                        amplitude_binindex2=0
                        do k=1,nbin
                           if(k.eq.1)then
                              if(amplitude_binselect2.lt.
     &                             Tlculm_dist_dom2(k))then
                                 amplitude_binindex2=k
                              endif
                           elseif(k.gt.1)then
                              if((amplitude_binselect2.lt.
     &                             Tlculm_dist_dom2(k))
     &                             .and.
     &                             (amplitude_binselect2.ge.
     &                             Tlculm_dist_dom2(k-1))
     &                             )then
                                 amplitude_binindex2=k
                              endif           
                           endif
                        enddo      
                        
                        if(amplitude_binindex2.gt.1)then
                           Flx_Rg=(rand3(r3)*(Tlcutoff_dist_dom2
     &                          (amplitude_binindex2)-
     &                          Tlcutoff_dist_dom2
     &                          (amplitude_binindex2-1))
     &                          +Tlcutoff_dist_dom2
     &                          (amplitude_binindex2-1))*0.1
                        elseif(amplitude_binindex2.eq.1)then
                           Flx_Rg=rand3(r3)*Tlcutoff_dist_dom2
     &                          (amplitude_binindex2)*0.1
                        endif
                        theta=rand3(r3)*pai
                        phai=rand3(r3)*2*pai
                        do k=1,Plm_res_num
                           Plm_B7_x_new0(Plm_index,j,k)=
     &                          Plm_B7_xo_new(Plm_index,j,k)
     &                          +Flx_Rg*sin(theta)*cos(phai)
                           Plm_B7_y_new0(Plm_index,j,k)=
     &                          Plm_B7_yo_new(Plm_index,j,k)
     &                          +Flx_Rg*sin(theta)*sin(phai)
                           Plm_B7_z_new0(Plm_index,j,k)=
     &                          Plm_B7_zo_new(Plm_index,j,k)
     &                          +Flx_Rg*cos(theta)
                        enddo
                     enddo

c>>>>>   random rotational flexibility

                     do j=1,Plm_MHC_num
                        
                        Plm_MHC_x_new(Plm_index,j,1)=
     &                       Plm_MHC_x_new0(Plm_index,j,1)
                        Plm_MHC_y_new(Plm_index,j,1)=
     &                       Plm_MHC_y_new0(Plm_index,j,1)
                        Plm_MHC_z_new(Plm_index,j,1)=
     &                       Plm_MHC_z_new0(Plm_index,j,1)

                        amplitude_binindex1=0
                        do k=1,nbin
                           if(k.eq.1)then
                              if(amplitude_binselect1.lt.
     &                             Rtculm_dist_dom1(k))then
                                 amplitude_binindex1=k
                              endif
                           elseif(k.gt.1)then
                              if((amplitude_binselect1.lt.
     &                             Rtculm_dist_dom1(k))
     &                             .and.
     &                             (amplitude_binselect1.ge.
     &                             Rtculm_dist_dom1(k-1))
     &                             )then
                                 amplitude_binindex1=k
                              endif           
                           endif
                        enddo       
                        
                        if(amplitude_binindex1.gt.1)then
                           Plm_Flx_RtRg=rand3(r3)*(Rtcutoff_dist_dom1
     &                          (amplitude_binindex1)-
     &                          Rtcutoff_dist_dom1
     &                          (amplitude_binindex1-1))
     &                          +Rtcutoff_dist_dom1
     &                          (amplitude_binindex1-1)
                        elseif(amplitude_binindex1.eq.1)then
                           Plm_Flx_RtRg=rand3(r3)*Rtcutoff_dist_dom1
     &                          (amplitude_binindex1)
                        endif

                        if(rand5(r5).lt.0.5)then
                           Plm_Flx_RtRg=-Plm_Flx_RtRg
                        endif

                        theta=Plm_Flx_RtRg*pai/180.0
                        phi=Plm_Flx_RtRg*pai/180.0
                        psai=Plm_Flx_RtRg*pai/180.0
                        
                        t(1,1)=cos(psai)*cos(phi)-
     &                       cos(theta)*sin(phi)*sin(psai)
                        t(1,2)=-sin(psai)*cos(phi)-
     &                       cos(theta)*sin(phi)*cos(psai)
                        t(1,3)=sin(theta)*sin(phi)
                        
                        t(2,1)=cos(psai)*sin(phi)+
     &                       cos(theta)*cos(phi)*sin(psai)
                        t(2,2)=-sin(psai)*sin(phi)+
     &                       cos(theta)*cos(phi)*cos(psai)
                        t(2,3)=-sin(theta)*cos(phi)
                        
                        t(3,1)=sin(psai)*sin(theta)
                        t(3,2)=cos(psai)*sin(theta)
                        t(3,3)=cos(theta)
                        
                        do k=2,Plm_res_num
                           Plm_MHC_x_new(Plm_index,j,k)=t(1,1)*
     &                          (Plm_MHC_x_new0(Plm_index,j,k)
     &                          -Plm_MHC_x_new(Plm_index,j,1))+
     &                          t(1,2)*(Plm_MHC_y_new0(Plm_index,j,k)
     &                          -Plm_MHC_y_new(Plm_index,j,1))
     &                          +t(1,3)*(Plm_MHC_z_new0(Plm_index,j,k)
     &                          -Plm_MHC_z_new(Plm_index,j,1))
     &                          +Plm_MHC_x_new(Plm_index,j,1)
                           Plm_MHC_y_new(Plm_index,j,k)=t(2,1)*
     &                          (Plm_MHC_x_new0(Plm_index,j,k)
     &                          -Plm_MHC_x_new(Plm_index,j,1))+
     &                          t(2,2)*(Plm_MHC_y_new0(Plm_index,j,k)
     &                          -Plm_MHC_y_new(Plm_index,j,1))
     &                          +t(2,3)*(Plm_MHC_z_new0(Plm_index,j,k)
     &                          -Plm_MHC_z_new(Plm_index,j,1))
     &                          +Plm_MHC_y_new(Plm_index,j,1)
                           Plm_MHC_z_new(Plm_index,j,k)=t(3,1)*
     &                          (Plm_MHC_x_new0(Plm_index,j,k)
     &                          -Plm_MHC_x_new(Plm_index,j,1))+
     &                          t(3,2)*(Plm_MHC_y_new0(Plm_index,j,k)
     &                          -Plm_MHC_y_new(Plm_index,j,1))
     &                          +t(3,3)*(Plm_MHC_z_new0(Plm_index,j,k)
     &                          -Plm_MHC_z_new(Plm_index,j,1))
     &                          +Plm_MHC_z_new(Plm_index,j,1)
                        enddo
                        
                     enddo

                     do j=1,Plm_B7_num
                        
                        Plm_B7_x_new(Plm_index,j,1)=
     &                       Plm_B7_x_new0(Plm_index,j,1)
                        Plm_B7_y_new(Plm_index,j,1)=
     &                       Plm_B7_y_new0(Plm_index,j,1)
                        Plm_B7_z_new(Plm_index,j,1)=
     &                       Plm_B7_z_new0(Plm_index,j,1)

                        amplitude_binselect2=rand5(r5)
                        amplitude_binindex2=0
                        do k=1,nbin
                           if(k.eq.1)then
                              if(amplitude_binselect2.lt.
     &                             Tlculm_dist_dom2(k))then
                                 amplitude_binindex2=k
                              endif
                           elseif(k.gt.1)then
                              if((amplitude_binselect2.lt.
     &                             Tlculm_dist_dom2(k))
     &                             .and.
     &                             (amplitude_binselect2.ge.
     &                             Tlculm_dist_dom2(k-1))
     &                             )then
                                 amplitude_binindex2=k
                              endif           
                           endif
                        enddo      
                        
                        if(amplitude_binindex2.gt.1)then
                           Plm_Flx_RtRg=rand3(r3)*(Rtcutoff_dist_dom2
     &                          (amplitude_binindex2)-
     &                          Rtcutoff_dist_dom2
     &                          (amplitude_binindex2-1))
     &                          +Rtcutoff_dist_dom2
     &                          (amplitude_binindex2-1)
                        elseif(amplitude_binindex2.eq.1)then
                           Plm_Flx_RtRg=rand3(r3)*Rtcutoff_dist_dom2
     &                          (amplitude_binindex2)
                        endif
                        
                        if(rand5(r5).lt.0.5)then
                           Plm_Flx_RtRg=-Plm_Flx_RtRg
                        endif

                        theta=Plm_Flx_RtRg*pai/180.0
                        phi=Plm_Flx_RtRg*pai/180.0
                        psai=Plm_Flx_RtRg*pai/180.0
                        
                        t(1,1)=cos(psai)*cos(phi)-
     &                       cos(theta)*sin(phi)*sin(psai)
                        t(1,2)=-sin(psai)*cos(phi)-
     &                       cos(theta)*sin(phi)*cos(psai)
                        t(1,3)=sin(theta)*sin(phi)
                        
                        t(2,1)=cos(psai)*sin(phi)+
     &                       cos(theta)*cos(phi)*sin(psai)
                        t(2,2)=-sin(psai)*sin(phi)+
     &                       cos(theta)*cos(phi)*cos(psai)
                        t(2,3)=-sin(theta)*cos(phi)
                        
                        t(3,1)=sin(psai)*sin(theta)
                        t(3,2)=cos(psai)*sin(theta)
                        t(3,3)=cos(theta)
                        
                        do k=2,Plm_res_num
                           Plm_B7_x_new(Plm_index,j,k)=t(1,1)*
     &                          (Plm_B7_x_new0(Plm_index,j,k)
     &                          -Plm_B7_x_new(Plm_index,j,1))+
     &                          t(1,2)*(Plm_B7_y_new0(Plm_index,j,k)
     &                          -Plm_B7_y_new(Plm_index,j,1))
     &                          +t(1,3)*(Plm_B7_z_new0(Plm_index,j,k)
     &                          -Plm_B7_z_new(Plm_index,j,1))
     &                          +Plm_B7_x_new(Plm_index,j,1)
                           Plm_B7_y_new(Plm_index,j,k)=t(2,1)*
     &                          (Plm_B7_x_new0(Plm_index,j,k)
     &                          -Plm_B7_x_new(Plm_index,j,1))+
     &                          t(2,2)*(Plm_B7_y_new0(Plm_index,j,k)
     &                          -Plm_B7_y_new(Plm_index,j,1))
     &                          +t(2,3)*(Plm_B7_z_new0(Plm_index,j,k)
     &                          -Plm_B7_z_new(Plm_index,j,1))
     &                          +Plm_B7_y_new(Plm_index,j,1)
                           Plm_B7_z_new(Plm_index,j,k)=t(3,1)*
     &                          (Plm_B7_x_new0(Plm_index,j,k)
     &                          -Plm_B7_x_new(Plm_index,j,1))+
     &                          t(3,2)*(Plm_B7_y_new0(Plm_index,j,k)
     &                          -Plm_B7_y_new(Plm_index,j,1))
     &                          +t(3,3)*(Plm_B7_z_new0(Plm_index,j,k)
     &                          -Plm_B7_z_new(Plm_index,j,1))
     &                          +Plm_B7_z_new(Plm_index,j,1)
                        enddo
                        
                     enddo

ccc>>>>   Check Collisions
                     collision_flag=0
                     do i=1,Plm_tot_num
                        if(Plm_index.ne.i)then

c                           cm1_a_x=0
c                           cm1_a_y=0
c                           cm1_a_z=0
c                           do k=1,Plm_MHC_num
c                              cm1_a_x=cm1_a_x+Plm_MHC_x_new
c     &                             (Plm_index,k,1)
c                              cm1_a_y=cm1_a_y+Plm_MHC_y_new
c     &                             (Plm_index,k,1)
c                              cm1_a_z=cm1_a_z+Plm_MHC_z_new
c     &                             (Plm_index,k,1)
c                           enddo
c                           cm1_a_x=cm1_a_x/real(Plm_MHC_num)
c                           cm1_a_y=cm1_a_y/real(Plm_MHC_num)
c                           cm1_a_z=cm1_a_z/real(Plm_MHC_num)
c                           cm1_b_x=0
c                           cm1_b_y=0
c                           cm1_b_z=0
c                           do k=1,Plm_MHC_num
c                              cm1_b_x=cm1_b_x+Plm_MHC_x_new
c     &                             (i,k,1)
c                              cm1_b_y=cm1_b_y+Plm_MHC_y_new
c     &                             (i,k,1)
c                              cm1_b_z=cm1_b_z+Plm_MHC_z_new
c     &                             (i,k,1)
c                           enddo
c                           cm1_b_x=cm1_b_x/real(Plm_MHC_num)
c                           cm1_b_y=cm1_b_y/real(Plm_MHC_num)
c                           cm1_b_z=cm1_b_z/real(Plm_MHC_num)
c                           dist=sqrt((cm1_b_x-cm1_a_x)**2+
c     &                          (cm1_b_y-cm1_a_y)**2+
c     &                          (cm1_b_x-cm1_a_z)**2)
c                           if(dist.le.
c     &                          (2*Plm_radius))then
c                              collision_flag=1
c                           endif

                           do j=1,Plm_MHC_num
                              do k=1,Plm_MHC_num
                                 dist=sqrt((Plm_MHC_x_new(i,j,1)-
     &                                Plm_MHC_x_new
     &                                (Plm_index,k,1))**2+
     &                                (Plm_MHC_y_new(i,j,1)-
     &                                Plm_MHC_y_new
     &                                (Plm_index,j,1))**2+
     &                                (Plm_MHC_z_new(i,j,1)-
     &                                Plm_MHC_z_new
     &                                (Plm_index,j,1))**2)
                                 if(dist.lt.(MHC_radius+MHC_radius))then
                                    collision_flag=1
                                 endif
                              enddo
                           enddo
                           do j=1,Plm_MHC_num
                              do k=1,Plm_B7_num
                                 dist=sqrt((Plm_B7_x_new(i,j,1)-
     &                                Plm_MHC_x_new
     &                                (Plm_index,k,1))**2+
     &                                (Plm_B7_y_new(i,j,1)-
     &                                Plm_MHC_y_new
     &                                (Plm_index,j,1))**2+
     &                                (Plm_B7_z_new(i,j,1)-
     &                                Plm_MHC_z_new
     &                                (Plm_index,j,1))**2)
                                 if(dist.lt.(B7_radius+MHC_radius))then
                                    collision_flag=1
                                 endif
                              enddo
                           enddo
                           do j=1,Plm_B7_num
                              do k=1,Plm_B7_num
                                 dist=sqrt((Plm_B7_x_new(i,j,1)-
     &                                Plm_B7_x_new
     &                                (Plm_index,k,1))**2+
     &                                (Plm_B7_y_new(i,j,1)-
     &                                Plm_B7_y_new
     &                                (Plm_index,j,1))**2+
     &                                (Plm_B7_z_new(i,j,1)-
     &                                Plm_B7_z_new
     &                                (Plm_index,j,1))**2)
                                 if(dist.lt.(B7_radius+B7_radius))then
                                    collision_flag=1
                                 endif
                              enddo
                           enddo
                           do j=1,Plm_B7_num
                              do k=1,Plm_MHC_num
                                 dist=sqrt((Plm_MHC_x_new(i,j,1)-
     &                                Plm_B7_x_new
     &                                (Plm_index,k,1))**2+
     &                                (Plm_MHC_y_new(i,j,1)-
     &                                Plm_B7_y_new
     &                                (Plm_index,j,1))**2+
     &                                (Plm_MHC_z_new(i,j,1)-
     &                                Plm_B7_z_new
     &                                (Plm_index,j,1))**2)
                                 if(dist.lt.(MHC_radius+B7_radius))then
                                    collision_flag=1
                                 endif
                              enddo
                           enddo
                        endif
                     enddo

                     do i=1,TCR_tot_num
                        do j=1,Plm_MHC_num
                           dist=sqrt((TCR_x_new(i,1)-
     &                          Plm_MHC_x_new(Plm_index,j,1))**2+
     &                          (TCR_y_new(i,1)-
     &                          Plm_MHC_y_new(Plm_index,j,1))**2+
     &                          (TCR_z_new(i,1)-
     &                          Plm_MHC_z_new(Plm_index,j,1))**2)
                           if(dist.lt.(TCR_radius+MHC_radius)
     &                          )then
                              collision_flag=1
                           endif
                        enddo
                        do j=1,Plm_B7_num
                           dist=sqrt((TCR_x_new(i,1)-
     &                          Plm_B7_x_new(Plm_index,j,1))**2+
     &                          (TCR_y_new(i,1)-
     &                          Plm_B7_y_new(Plm_index,j,1))**2+
     &                          (TCR_z_new(i,1)-
     &                          Plm_B7_z_new(Plm_index,j,1))**2)
                           if(dist.lt.(TCR_radius+B7_radius)
     &                          )then
                              collision_flag=1
                           endif
                        enddo
                     enddo


                     do i=1,CD_tot_num
                        do j=1,Plm_MHC_num
                           dist=sqrt((CD_x_new(i,1)-
     &                          Plm_MHC_x_new(Plm_index,j,1))**2+
     &                          (CD_y_new(i,1)-
     &                          Plm_MHC_y_new(Plm_index,j,1))**2+
     &                          (CD_z_new(i,1)-
     &                          Plm_MHC_z_new(Plm_index,j,1))**2)
                           if(dist.lt.(CD_radius+MHC_radius)
     &                          )then
                              collision_flag=1
                           endif
                        enddo
                        do j=1,Plm_B7_num
                           dist=sqrt((CD_x_new(i,1)-
     &                          Plm_B7_x_new(Plm_index,j,1))**2+
     &                          (CD_y_new(i,1)-
     &                          Plm_B7_y_new(Plm_index,j,1))**2+
     &                          (CD_z_new(i,1)-
     &                          Plm_B7_z_new(Plm_index,j,1))**2)
                           if(dist.lt.(CD_radius+B7_radius)
     &                          )then
                              collision_flag=1
                           endif
                        enddo
                     enddo

                     if(collision_flag.eq.1)then
                        do i=1,Plm_MHC_num
                           do j=1,Plm_res_num
                              Plm_MHC_x_new(Plm_index,i,j)=
     &                             Plm_MHC_x(Plm_index,i,j)
                              Plm_MHC_y_new(Plm_index,i,j)=
     &                             Plm_MHC_y(Plm_index,i,j)
                              Plm_MHC_z_new(Plm_index,i,j)=
     &                             Plm_MHC_z(Plm_index,i,j)
                              Plm_MHC_xo_new(Plm_index,i,j)=
     &                             Plm_MHC_xo(Plm_index,i,j)
                              Plm_MHC_yo_new(Plm_index,i,j)=
     &                             Plm_MHC_yo(Plm_index,i,j)
                              Plm_MHC_zo_new(Plm_index,i,j)=
     &                             Plm_MHC_zo(Plm_index,i,j)
                           enddo
                        enddo
                        do i=1,Plm_B7_num
                           do j=1,Plm_res_num
                              Plm_B7_x_new(Plm_index,i,j)=
     &                             Plm_B7_x(Plm_index,i,j)
                              Plm_B7_y_new(Plm_index,i,j)=
     &                             Plm_B7_y(Plm_index,i,j)
                              Plm_B7_z_new(Plm_index,i,j)=
     &                             Plm_B7_z(Plm_index,i,j)
                              Plm_B7_xo_new(Plm_index,i,j)=
     &                             Plm_B7_xo(Plm_index,i,j)
                              Plm_B7_yo_new(Plm_index,i,j)=
     &                             Plm_B7_yo(Plm_index,i,j)
                              Plm_B7_zo_new(Plm_index,i,j)=
     &                             Plm_B7_zo(Plm_index,i,j)
                           enddo
                        enddo
                     endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc>>>>   elseif it is in the complex with TRC or CD28
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

                  elseif(Plm_bound_flag.eq.1)then

                     Prob_Diff=(4*complex_D*time_step)/
     &                    (distance_step**2)
 
                     if(rand3(r3).lt.Prob_Diff)then

                        phai=rand3(r3)*2*pai
                        distance_amp=2*rand4(r4)*distance_step
                        do i=1,Plm_MHC_num
                           do j=1,Plm_res_num
                              Plm_MHC_xo_new0(Plm_index,i,j)=
     &                             Plm_MHC_xo(Plm_index,i,j)+
     &                             distance_amp*cos(phai)
                              Plm_MHC_yo_new0(Plm_index,i,j)=
     &                             Plm_MHC_yo(Plm_index,i,j)+
     &                             distance_amp*sin(phai)
                              Plm_MHC_zo_new0(Plm_index,i,j)=
     &                             Plm_MHC_zo(Plm_index,i,j)
                           enddo
                           do j=1,Plm_res_num
                              Plm_MHC_x_new0(Plm_index,i,j)=
     &                             Plm_MHC_x(Plm_index,i,j)+
     &                             distance_amp*cos(phai)
                              Plm_MHC_y_new0(Plm_index,i,j)=
     &                             Plm_MHC_y(Plm_index,i,j)+
     &                             distance_amp*sin(phai)
                              Plm_MHC_z_new0(Plm_index,i,j)=
     &                             Plm_MHC_z(Plm_index,i,j)
                           enddo
                           if(Plm_MHC_status_new(Plm_index,i).eq.1)then
                              do j=1,TCR_res_num
                                 TCR_x_new0(Plm_MHC_idx_new
     &                                (Plm_index,i),j)=
     &                                TCR_x(Plm_MHC_idx_new
     &                                (Plm_index,i),j)+
     &                                distance_amp*cos(phai)
                                 TCR_y_new0(Plm_MHC_idx_new
     &                                (Plm_index,i),j)=
     &                                TCR_y(Plm_MHC_idx_new
     &                                (Plm_index,i),j)+
     &                                distance_amp*sin(phai)
                                 TCR_z_new0(Plm_MHC_idx_new
     &                                (Plm_index,i),j)=
     &                                TCR_z(Plm_MHC_idx_new
     &                                (Plm_index,i),j)
                              enddo
                           endif
                        enddo
                        do i=1,Plm_B7_num
                           do j=1,Plm_res_num
                              Plm_B7_xo_new0(Plm_index,i,j)=
     &                             Plm_B7_xo(Plm_index,i,j)+
     &                             distance_amp*cos(phai)
                              Plm_B7_yo_new0(Plm_index,i,j)=
     &                             Plm_B7_yo(Plm_index,i,j)+
     &                             distance_amp*sin(phai)
                              Plm_B7_zo_new0(Plm_index,i,j)=
     &                             Plm_B7_zo(Plm_index,i,j)
                           enddo
                           do j=1,Plm_res_num
                              Plm_B7_x_new0(Plm_index,i,j)=
     &                             Plm_B7_x(Plm_index,i,j)+
     &                             distance_amp*cos(phai)
                              Plm_B7_y_new0(Plm_index,i,j)=
     &                             Plm_B7_y(Plm_index,i,j)+
     &                             distance_amp*sin(phai)
                              Plm_B7_z_new0(Plm_index,i,j)=
     &                             Plm_B7_z(Plm_index,i,j)
                           enddo
                           if(Plm_B7_status_new(Plm_index,i).eq.1)then
                              do j=1,CD_res_num
                                 CD_x_new0(Plm_B7_idx_new
     &                                (Plm_index,i),j)=
     &                                CD_x(Plm_B7_idx_new
     &                                (Plm_index,i),j)+
     &                                distance_amp*cos(phai)
                                 CD_y_new0(Plm_B7_idx_new
     &                                (Plm_index,i),j)=
     &                                CD_y(Plm_B7_idx_new
     &                                (Plm_index,i),j)+
     &                                distance_amp*sin(phai)
                                 CD_z_new0(Plm_B7_idx_new
     &                                (Plm_index,i),j)=
     &                                CD_z(Plm_B7_idx_new
     &                                (Plm_index,i),j)
                              enddo
                           endif
                        enddo                        
                        temp_PB_x=0
                        temp_PB_y=0
                        do i=1,Plm_MHC_num
                           temp_PB_x=temp_PB_x+
     &                          Plm_MHC_xo_new0(Plm_index,i,1)
                           temp_PB_y=temp_PB_y+
     &                          Plm_MHC_yo_new0(Plm_index,i,1)
                        enddo
                        temp_PB_x=temp_PB_x/real(Plm_MHC_num)
                        temp_PB_y=temp_PB_y/real(Plm_MHC_num)
                        PB_x=cell_range_x*anint
     &                       (temp_PB_x/cell_range_x)
                        PB_y=cell_range_y*anint
     &                       (temp_PB_y/cell_range_y)
c                        if(Plm_index.eq.18)then
c                           print*,'PB',Plm_index,
c     &                          Plm_MHC_xo_new0(Plm_index,1,1),
c     &                          Plm_MHC_x_new0(Plm_index,1,1),PB_x
c                           print*,'PB',Plm_index,
c     &                          Plm_MHC_yo_new0(Plm_index,1,1),
c     &                          Plm_MHC_y_new0(Plm_index,1,1),PB_y
c                        endif
                        do i=1,Plm_MHC_num
                           do j=1,Plm_res_num
                              Plm_MHC_xo_new0(Plm_index,i,j)=
     &                             Plm_MHC_xo_new0(Plm_index,i,j)-PB_x     
                              Plm_MHC_yo_new0(Plm_index,i,j)=
     &                             Plm_MHC_yo_new0(Plm_index,i,j)-PB_y
                           enddo
                           do j=1,Plm_res_num
                              Plm_MHC_x_new0(Plm_index,i,j)=
     &                             Plm_MHC_x_new0(Plm_index,i,j)-PB_x     
                              Plm_MHC_y_new0(Plm_index,i,j)=
     &                             Plm_MHC_y_new0(Plm_index,i,j)-PB_y
                           enddo
                           if(Plm_MHC_status_new(Plm_index,i).eq.1)then
                              do j=1,TCR_res_num
                                 TCR_x_new0(Plm_MHC_idx_new
     &                                (Plm_index,i),j)=
     &                                TCR_x_new0(Plm_MHC_idx_new
     &                                (Plm_index,i),j)-PB_x
                                 TCR_y_new0(Plm_MHC_idx_new
     &                                (Plm_index,i),j)=
     &                                TCR_y_new0(Plm_MHC_idx_new
     &                                (Plm_index,i),j)-PB_y
                              enddo
                           endif
                        enddo
                        do i=1,Plm_B7_num
                           do j=1,Plm_res_num
                              Plm_B7_xo_new0(Plm_index,i,j)=
     &                             Plm_B7_xo_new0(Plm_index,i,j)-PB_x     
                              Plm_B7_yo_new0(Plm_index,i,j)=
     &                             Plm_B7_yo_new0(Plm_index,i,j)-PB_y
                           enddo
                           do j=1,Plm_res_num
                              Plm_B7_x_new0(Plm_index,i,j)=
     &                             Plm_B7_x_new0(Plm_index,i,j)-PB_x     
                              Plm_B7_y_new0(Plm_index,i,j)=
     &                             Plm_B7_y_new0(Plm_index,i,j)-PB_y
                           enddo
                           if(Plm_B7_status_new(Plm_index,i).eq.1)then
                              do j=1,CD_res_num
                                 CD_x_new0(Plm_B7_idx_new
     &                                (Plm_index,i),j)=
     &                                CD_x_new0(Plm_B7_idx_new
     &                                (Plm_index,i),j)-PB_x
                                 CD_y_new0(Plm_B7_idx_new
     &                                (Plm_index,i),j)=
     &                                CD_y_new0(Plm_B7_idx_new
     &                                (Plm_index,i),j)-PB_y
                              enddo
                           endif
                        enddo
                     else
                        do i=1,Plm_MHC_num
                           do j=1,Plm_res_num
                              Plm_MHC_xo_new0(Plm_index,i,j)=
     &                             Plm_MHC_xo(Plm_index,i,j)
                              Plm_MHC_yo_new0(Plm_index,i,j)=
     &                             Plm_MHC_yo(Plm_index,i,j)
                              Plm_MHC_zo_new0(Plm_index,i,j)=
     &                             Plm_MHC_zo(Plm_index,i,j)
                           enddo
                           do j=1,Plm_res_num
                              Plm_MHC_x_new0(Plm_index,i,j)=
     &                             Plm_MHC_x(Plm_index,i,j)
                              Plm_MHC_y_new0(Plm_index,i,j)=
     &                             Plm_MHC_y(Plm_index,i,j)
                              Plm_MHC_z_new0(Plm_index,i,j)=
     &                             Plm_MHC_z(Plm_index,i,j)
                           enddo
                           if(Plm_MHC_status_new(Plm_index,i).eq.1)then
                              do j=1,TCR_res_num
                                 TCR_x_new0(Plm_MHC_idx_new
     &                                (Plm_index,i),j)=
     &                                TCR_x(Plm_MHC_idx_new
     &                                (Plm_index,i),j)
                                 TCR_y_new0(Plm_MHC_idx_new
     &                                (Plm_index,i),j)=
     &                                TCR_y(Plm_MHC_idx_new
     &                                (Plm_index,i),j)
                                 TCR_z_new0(Plm_MHC_idx_new
     &                                (Plm_index,i),j)=
     &                                TCR_z(Plm_MHC_idx_new
     &                                (Plm_index,i),j)
                              enddo
                           endif
                        enddo
                        do i=1,Plm_B7_num
                           do j=1,Plm_res_num
                              Plm_B7_xo_new0(Plm_index,i,j)=
     &                             Plm_B7_xo(Plm_index,i,j)
                              Plm_B7_yo_new0(Plm_index,i,j)=
     &                             Plm_B7_yo(Plm_index,i,j)
                              Plm_B7_zo_new0(Plm_index,i,j)=
     &                             Plm_B7_zo(Plm_index,i,j)
                           enddo
                           do j=1,Plm_res_num
                              Plm_B7_x_new0(Plm_index,i,j)=
     &                             Plm_B7_x(Plm_index,i,j)
                              Plm_B7_y_new0(Plm_index,i,j)=
     &                             Plm_B7_y(Plm_index,i,j)
                              Plm_B7_z_new0(Plm_index,i,j)=
     &                             Plm_B7_z(Plm_index,i,j)
                           enddo
                           if(Plm_B7_status_new(Plm_index,i).eq.1)then
                              do j=1,CD_res_num
                                 CD_x_new0(Plm_B7_idx_new
     &                                (Plm_index,i),j)=
     &                                CD_x(Plm_B7_idx_new
     &                                (Plm_index,i),j)
                                 CD_y_new0(Plm_B7_idx_new
     &                                (Plm_index,i),j)=
     &                                CD_y(Plm_B7_idx_new
     &                                (Plm_index,i),j)
                                 CD_z_new0(Plm_B7_idx_new
     &                                (Plm_index,i),j)=
     &                                CD_z(Plm_B7_idx_new
     &                                (Plm_index,i),j)
                              enddo
                           endif
                        enddo
                        temp_PB_x=0
                        temp_PB_y=0
                        do i=1,Plm_MHC_num
                           temp_PB_x=temp_PB_x+
     &                          Plm_MHC_xo_new0(Plm_index,i,1)
                           temp_PB_y=temp_PB_y+
     &                          Plm_MHC_yo_new0(Plm_index,i,1)
                        enddo
                        temp_PB_x=temp_PB_x/real(Plm_MHC_num)
                        temp_PB_y=temp_PB_y/real(Plm_MHC_num)
                        PB_x=cell_range_x*anint
     &                       (temp_PB_x/cell_range_x)
                        PB_y=cell_range_y*anint
     &                       (temp_PB_y/cell_range_y)
                        do i=1,Plm_MHC_num
                           do j=1,Plm_res_num
                              Plm_MHC_xo_new0(Plm_index,i,j)=
     &                             Plm_MHC_xo_new0(Plm_index,i,j)-PB_x     
                              Plm_MHC_yo_new0(Plm_index,i,j)=
     &                             Plm_MHC_yo_new0(Plm_index,i,j)-PB_y
                           enddo
                           do j=1,Plm_res_num
                              Plm_MHC_x_new0(Plm_index,i,j)=
     &                             Plm_MHC_x_new0(Plm_index,i,j)-PB_x     
                              Plm_MHC_y_new0(Plm_index,i,j)=
     &                             Plm_MHC_y_new0(Plm_index,i,j)-PB_y
                           enddo
                           if(Plm_MHC_status_new(Plm_index,i).eq.1)then
                              do j=1,TCR_res_num
                                 TCR_x_new0(Plm_MHC_idx_new
     &                                (Plm_index,i),j)=
     &                                TCR_x_new0(Plm_MHC_idx_new
     &                                (Plm_index,i),j)-PB_x
                                 TCR_y_new0(Plm_MHC_idx_new
     &                                (Plm_index,i),j)=
     &                                TCR_y_new0(Plm_MHC_idx_new
     &                                (Plm_index,i),j)-PB_y
                              enddo
                           endif
                        enddo
                        do i=1,Plm_B7_num
                           do j=1,Plm_res_num
                              Plm_B7_xo_new0(Plm_index,i,j)=
     &                             Plm_B7_xo_new0(Plm_index,i,j)-PB_x     
                              Plm_B7_yo_new0(Plm_index,i,j)=
     &                             Plm_B7_yo_new0(Plm_index,i,j)-PB_y
                           enddo
                           do j=1,Plm_res_num
                              Plm_B7_x_new0(Plm_index,i,j)=
     &                             Plm_B7_x_new0(Plm_index,i,j)-PB_x     
                              Plm_B7_y_new0(Plm_index,i,j)=
     &                             Plm_B7_y_new0(Plm_index,i,j)-PB_y
                           enddo
                           if(Plm_B7_status_new(Plm_index,i).eq.1)then
                              do j=1,CD_res_num
                                 CD_x_new0(Plm_B7_idx_new
     &                                (Plm_index,i),j)=
     &                                CD_x_new0(Plm_B7_idx_new
     &                                (Plm_index,i),j)-PB_x
                                 CD_y_new0(Plm_B7_idx_new
     &                                (Plm_index,i),j)=
     &                                CD_y_new0(Plm_B7_idx_new
     &                                (Plm_index,i),j)-PB_y
                              enddo
                           endif
                        enddo
                     endif
                        
c>>>>>>>>> Rotate the complex

                     rot_angle=2*rand3(r3)*Complex_rot_D
     &                    -Complex_rot_D
                     x_axis_bg=0
                     y_axis_bg=0
                     x_axis_ed=0
                     y_axis_ed=0
                     do i=1,Plm_MHC_num
                        x_axis_bg=x_axis_bg+
     &                       Plm_MHC_xo_new0(Plm_index,i,1)
                        y_axis_bg=y_axis_bg+
     &                       Plm_MHC_yo_new0(Plm_index,i,1)
                        x_axis_ed=x_axis_ed+
     &                       Plm_MHC_xo_new0(Plm_index,i,1)
                        y_axis_ed=y_axis_ed+
     &                       Plm_MHC_yo_new0(Plm_index,i,1)
                     enddo
                     do i=1,Plm_B7_num
                        x_axis_bg=x_axis_bg+
     &                       Plm_B7_xo_new0(Plm_index,i,1)
                        y_axis_bg=y_axis_bg+
     &                       Plm_B7_yo_new0(Plm_index,i,1)
                        x_axis_ed=x_axis_ed+
     &                       Plm_B7_xo_new0(Plm_index,i,1)
                        y_axis_ed=y_axis_ed+
     &                       Plm_B7_yo_new0(Plm_index,i,1)
                     enddo
                     x_axis_bg=x_axis_bg/real(Plm_MHC_num+Plm_B7_num)
                     y_axis_bg=y_axis_bg/real(Plm_MHC_num+Plm_B7_num)
                     x_axis_ed=x_axis_ed/real(Plm_MHC_num+Plm_B7_num)
                     y_axis_ed=y_axis_ed/real(Plm_MHC_num+Plm_B7_num)
                     z_axis_bg=-100.0
                     z_axis_ed=100.0
                     do i=1,Plm_MHC_num
                        do j=1,Plm_res_num
                           x_bf_rot=Plm_MHC_xo_new0(Plm_index,i,j)
                           y_bf_rot=Plm_MHC_yo_new0(Plm_index,i,j)
                           z_bf_rot=Plm_MHC_zo_new0(Plm_index,i,j)
                           x_af_rot=0
                           y_af_rot=0
                           z_af_rot=0
                           call rot_along_axis(rot_angle,
     &                          x_axis_bg,y_axis_bg,z_axis_bg,
     &                          x_axis_ed,y_axis_ed,z_axis_ed,
     &                          x_bf_rot,y_bf_rot,z_bf_rot,
     &                          x_af_rot,y_af_rot,z_af_rot)
                           Plm_MHC_xo_new(Plm_index,i,j)=x_af_rot
                           Plm_MHC_yo_new(Plm_index,i,j)=y_af_rot
                           Plm_MHC_zo_new(Plm_index,i,j)=z_af_rot
                        enddo
                        do j=1,Plm_res_num
                           x_bf_rot=Plm_MHC_x_new0(Plm_index,i,j)
                           y_bf_rot=Plm_MHC_y_new0(Plm_index,i,j)
                           z_bf_rot=Plm_MHC_z_new0(Plm_index,i,j)
                           x_af_rot=0
                           y_af_rot=0
                           z_af_rot=0
                           call rot_along_axis(rot_angle,
     &                          x_axis_bg,y_axis_bg,z_axis_bg,
     &                          x_axis_ed,y_axis_ed,z_axis_ed,
     &                          x_bf_rot,y_bf_rot,z_bf_rot,
     &                          x_af_rot,y_af_rot,z_af_rot)
                           Plm_MHC_x_new(Plm_index,i,j)=x_af_rot
                           Plm_MHC_y_new(Plm_index,i,j)=y_af_rot
                           Plm_MHC_z_new(Plm_index,i,j)=z_af_rot
                        enddo
                        if(Plm_MHC_status_new(Plm_index,i).eq.1)then
                           do j=1,TCR_res_num
                              x_bf_rot=TCR_x_new0(Plm_MHC_idx_new
     &                             (Plm_index,i),j)
                              y_bf_rot=TCR_y_new0(Plm_MHC_idx_new
     &                             (Plm_index,i),j)
                              z_bf_rot=TCR_z_new0(Plm_MHC_idx_new
     &                             (Plm_index,i),j)
                              call rot_along_axis(rot_angle,
     &                             x_axis_bg,y_axis_bg,z_axis_bg,
     &                             x_axis_ed,y_axis_ed,z_axis_ed,
     &                             x_bf_rot,y_bf_rot,z_bf_rot,
     &                             x_af_rot,y_af_rot,z_af_rot)
                              TCR_x_new(Plm_MHC_idx_new
     &                             (Plm_index,i),j)=x_af_rot
                              TCR_y_new(Plm_MHC_idx_new
     &                             (Plm_index,i),j)=y_af_rot
                              TCR_z_new(Plm_MHC_idx_new
     &                             (Plm_index,i),j)=z_af_rot
                           enddo
                        endif
                     enddo
                     do i=1,Plm_B7_num
                        do j=1,Plm_res_num
                           x_bf_rot=Plm_B7_xo_new0(Plm_index,i,j)
                           y_bf_rot=Plm_B7_yo_new0(Plm_index,i,j)
                           z_bf_rot=Plm_B7_zo_new0(Plm_index,i,j)
                           x_af_rot=0
                           y_af_rot=0
                           z_af_rot=0
                           call rot_along_axis(rot_angle,
     &                          x_axis_bg,y_axis_bg,z_axis_bg,
     &                          x_axis_ed,y_axis_ed,z_axis_ed,
     &                          x_bf_rot,y_bf_rot,z_bf_rot,
     &                          x_af_rot,y_af_rot,z_af_rot)
                           Plm_B7_xo_new(Plm_index,i,j)=x_af_rot
                           Plm_B7_yo_new(Plm_index,i,j)=y_af_rot
                           Plm_B7_zo_new(Plm_index,i,j)=z_af_rot
                        enddo
                        do j=1,Plm_res_num
                           x_bf_rot=Plm_B7_x_new0(Plm_index,i,j)
                           y_bf_rot=Plm_B7_y_new0(Plm_index,i,j)
                           z_bf_rot=Plm_B7_z_new0(Plm_index,i,j)
                           x_af_rot=0
                           y_af_rot=0
                           z_af_rot=0
                           call rot_along_axis(rot_angle,
     &                          x_axis_bg,y_axis_bg,z_axis_bg,
     &                          x_axis_ed,y_axis_ed,z_axis_ed,
     &                          x_bf_rot,y_bf_rot,z_bf_rot,
     &                          x_af_rot,y_af_rot,z_af_rot)
                           Plm_B7_x_new(Plm_index,i,j)=x_af_rot
                           Plm_B7_y_new(Plm_index,i,j)=y_af_rot
                           Plm_B7_z_new(Plm_index,i,j)=z_af_rot
                        enddo
                        if(Plm_B7_status_new(Plm_index,i).eq.1)then
                           do j=1,CD_res_num
                              x_bf_rot=CD_x_new0(Plm_B7_idx_new
     &                             (Plm_index,i),j)
                              y_bf_rot=CD_y_new0(Plm_B7_idx_new
     &                             (Plm_index,i),j)
                              z_bf_rot=CD_z_new0(Plm_B7_idx_new
     &                             (Plm_index,i),j)
                              call rot_along_axis(rot_angle,
     &                             x_axis_bg,y_axis_bg,z_axis_bg,
     &                             x_axis_ed,y_axis_ed,z_axis_ed,
     &                             x_bf_rot,y_bf_rot,z_bf_rot,
     &                             x_af_rot,y_af_rot,z_af_rot)
                              CD_x_new(Plm_B7_idx_new
     &                             (Plm_index,i),j)=x_af_rot
                              CD_y_new(Plm_B7_idx_new
     &                             (Plm_index,i),j)=y_af_rot
                              CD_z_new(Plm_B7_idx_new
     &                             (Plm_index,i),j)=z_af_rot
                           enddo
                        endif
                     enddo

c>>>>>   random translational flexibility

                     do j=1,Plm_MHC_num
                        if(Plm_MHC_status_new(Plm_index,j).eq.0)then

                           amplitude_binselect1=rand5(r5)
                           amplitude_binindex1=0
                           do k=1,nbin
                              if(k.eq.1)then
                                 if(amplitude_binselect1.lt.
     &                                Tlculm_dist_dom1(k))then
                                    amplitude_binindex1=k
                                 endif
                              elseif(k.gt.1)then
                                 if((amplitude_binselect1.lt.
     &                                Tlculm_dist_dom1(k))
     &                                .and.
     &                                (amplitude_binselect1.ge.
     &                                Tlculm_dist_dom1(k-1))
     &                                )then
                                    amplitude_binindex1=k
                                 endif           
                              endif
                           enddo       
                           
                           if(amplitude_binindex1.gt.1)then
                              Flx_Rg=(rand3(r3)*(Tlcutoff_dist_dom1
     &                             (amplitude_binindex1)-
     &                             Tlcutoff_dist_dom1
     &                             (amplitude_binindex1-1))
     &                             +Tlcutoff_dist_dom1
     &                             (amplitude_binindex1-1))*0.1
                           elseif(amplitude_binindex1.eq.1)then
                              Flx_Rg=rand3(r3)*Tlcutoff_dist_dom1
     &                             (amplitude_binindex1)*0.1
                           endif
                           theta=rand3(r3)*pai
                           phai=rand3(r3)*2*pai

                           do k=1,Plm_res_num
                              Plm_MHC_x_new0(Plm_index,j,k)=
     &                             Plm_MHC_xo_new(Plm_index,j,k)
     &                             +Flx_Rg*sin(theta)*cos(phai)
                              Plm_MHC_y_new0(Plm_index,j,k)=
     &                             Plm_MHC_yo_new(Plm_index,j,k)
     &                             +Flx_Rg*sin(theta)*sin(phai)
                              Plm_MHC_z_new0(Plm_index,j,k)=
     &                             Plm_MHC_zo_new(Plm_index,j,k)
     &                             +Flx_Rg*cos(theta)
                           enddo
                        endif
                     enddo
                     
                     do j=1,Plm_B7_num
                        if(Plm_B7_status_new(Plm_index,j).eq.0)then

                           amplitude_binselect2=rand5(r5)
                           amplitude_binindex2=0
                           do k=1,nbin
                              if(k.eq.1)then
                                 if(amplitude_binselect2.lt.
     &                                Tlculm_dist_dom2(k))then
                                    amplitude_binindex2=k
                                 endif
                              elseif(k.gt.1)then
                                 if((amplitude_binselect2.lt.
     &                                Tlculm_dist_dom2(k))
     &                                .and.
     &                                (amplitude_binselect2.ge.
     &                                Tlculm_dist_dom2(k-1))
     &                                )then
                                    amplitude_binindex2=k
                                 endif           
                              endif
                           enddo      
                           
                           if(amplitude_binindex2.gt.1)then
                              Flx_Rg=(rand3(r3)*(Tlcutoff_dist_dom2
     &                             (amplitude_binindex2)-
     &                             Tlcutoff_dist_dom2
     &                             (amplitude_binindex2-1))
     &                             +Tlcutoff_dist_dom2
     &                             (amplitude_binindex2-1))*0.1
                           elseif(amplitude_binindex2.eq.1)then
                              Flx_Rg=rand3(r3)*Tlcutoff_dist_dom2
     &                             (amplitude_binindex2)*0.1
                           endif
                           theta=rand3(r3)*pai
                           phai=rand3(r3)*2*pai

                           do k=1,Plm_res_num
                              Plm_B7_x_new0(Plm_index,j,k)=
     &                             Plm_B7_xo_new(Plm_index,j,k)
     &                             +Flx_Rg*sin(theta)*cos(phai)
                              Plm_B7_y_new0(Plm_index,j,k)=
     &                             Plm_B7_yo_new(Plm_index,j,k)
     &                             +Flx_Rg*sin(theta)*sin(phai)
                              Plm_B7_z_new0(Plm_index,j,k)=
     &                             Plm_B7_zo_new(Plm_index,j,k)
     &                             +Flx_Rg*cos(theta)
                           enddo
                        endif
                     enddo

c>>>>>   random rotational flexibility

                     do j=1,Plm_MHC_num
                        
                        if(Plm_MHC_status_new(Plm_index,j).eq.0)then

                           Plm_MHC_x_new(Plm_index,j,1)=
     &                          Plm_MHC_x_new0(Plm_index,j,1)
                           Plm_MHC_y_new(Plm_index,j,1)=
     &                          Plm_MHC_y_new0(Plm_index,j,1)
                           Plm_MHC_z_new(Plm_index,j,1)=
     &                          Plm_MHC_z_new0(Plm_index,j,1)

                           amplitude_binindex1=0
                           do k=1,nbin
                              if(k.eq.1)then
                                 if(amplitude_binselect1.lt.
     &                                Rtculm_dist_dom1(k))then
                                    amplitude_binindex1=k
                                 endif
                              elseif(k.gt.1)then
                                 if((amplitude_binselect1.lt.
     &                                Rtculm_dist_dom1(k))
     &                                .and.
     &                                (amplitude_binselect1.ge.
     &                                Rtculm_dist_dom1(k-1))
     &                                )then
                                    amplitude_binindex1=k
                                 endif           
                              endif
                           enddo       
                           
                           if(amplitude_binindex1.gt.1)then
                              Plm_Flx_RtRg=rand3(r3)*(Rtcutoff_dist_dom1
     &                             (amplitude_binindex1)-
     &                             Rtcutoff_dist_dom1
     &                             (amplitude_binindex1-1))
     &                             +Rtcutoff_dist_dom1
     &                             (amplitude_binindex1-1)
                           elseif(amplitude_binindex1.eq.1)then
                              Plm_Flx_RtRg=rand3(r3)*Rtcutoff_dist_dom1
     &                             (amplitude_binindex1)
                           endif
                           
                           if(rand5(r5).lt.0.5)then
                              Plm_Flx_RtRg=-Plm_Flx_RtRg
                           endif

                           
                           theta=Plm_Flx_RtRg*pai/180.0
                           phi=Plm_Flx_RtRg*pai/180.0
                           psai=Plm_Flx_RtRg*pai/180.0
                           
                           t(1,1)=cos(psai)*cos(phi)-
     &                          cos(theta)*sin(phi)*sin(psai)
                           t(1,2)=-sin(psai)*cos(phi)-
     &                          cos(theta)*sin(phi)*cos(psai)
                           t(1,3)=sin(theta)*sin(phi)
                           
                           t(2,1)=cos(psai)*sin(phi)+
     &                          cos(theta)*cos(phi)*sin(psai)
                           t(2,2)=-sin(psai)*sin(phi)+
     &                          cos(theta)*cos(phi)*cos(psai)
                           t(2,3)=-sin(theta)*cos(phi)
                           
                           t(3,1)=sin(psai)*sin(theta)
                           t(3,2)=cos(psai)*sin(theta)
                           t(3,3)=cos(theta)
                        
                           do k=2,Plm_res_num
                              Plm_MHC_x_new(Plm_index,j,k)=t(1,1)*
     &                             (Plm_MHC_x_new0(Plm_index,j,k)
     &                             -Plm_MHC_x_new(Plm_index,j,1))+
     &                             t(1,2)*(Plm_MHC_y_new0(Plm_index,j,k)
     &                             -Plm_MHC_y_new(Plm_index,j,1))
     &                             +t(1,3)*(Plm_MHC_z_new0
     &                             (Plm_index,j,k)
     &                             -Plm_MHC_z_new(Plm_index,j,1))
     &                             +Plm_MHC_x_new(Plm_index,j,1)
                              Plm_MHC_y_new(Plm_index,j,k)=t(2,1)*
     &                             (Plm_MHC_x_new0(Plm_index,j,k)
     &                             -Plm_MHC_x_new(Plm_index,j,1))+
     &                             t(2,2)*(Plm_MHC_y_new0(Plm_index,j,k)
     &                             -Plm_MHC_y_new(Plm_index,j,1))
     &                             +t(2,3)*(Plm_MHC_z_new0
     &                             (Plm_index,j,k)
     &                             -Plm_MHC_z_new(Plm_index,j,1))
     &                             +Plm_MHC_y_new(Plm_index,j,1)
                              Plm_MHC_z_new(Plm_index,j,k)=t(3,1)*
     &                             (Plm_MHC_x_new0(Plm_index,j,k)
     &                             -Plm_MHC_x_new(Plm_index,j,1))+
     &                             t(3,2)*(Plm_MHC_y_new0(Plm_index,j,k)
     &                             -Plm_MHC_y_new(Plm_index,j,1))
     &                             +t(3,3)*(Plm_MHC_z_new0
     &                             (Plm_index,j,k)
     &                             -Plm_MHC_z_new(Plm_index,j,1))
     &                             +Plm_MHC_z_new(Plm_index,j,1)
                           enddo

                        endif
                        
                     enddo

                     do j=1,Plm_B7_num
                        
                        if(Plm_B7_status_new(Plm_index,j).eq.0)then

                           Plm_B7_x_new(Plm_index,j,1)=
     &                          Plm_B7_x_new0(Plm_index,j,1)
                           Plm_B7_y_new(Plm_index,j,1)=
     &                          Plm_B7_y_new0(Plm_index,j,1)
                           Plm_B7_z_new(Plm_index,j,1)=
     &                          Plm_B7_z_new0(Plm_index,j,1)

                           amplitude_binselect2=rand5(r5)
                           amplitude_binindex2=0
                           do k=1,nbin
                              if(k.eq.1)then
                                 if(amplitude_binselect2.lt.
     &                                Tlculm_dist_dom2(k))then
                                    amplitude_binindex2=k
                                 endif
                              elseif(k.gt.1)then
                                 if((amplitude_binselect2.lt.
     &                                Tlculm_dist_dom2(k))
     &                                .and.
     &                                (amplitude_binselect2.ge.
     &                                Tlculm_dist_dom2(k-1))
     &                                )then
                                    amplitude_binindex2=k
                                 endif           
                              endif
                           enddo      
                           
                           if(amplitude_binindex2.gt.1)then
                              Plm_Flx_RtRg=rand3(r3)*(Rtcutoff_dist_dom2
     &                             (amplitude_binindex2)-
     &                             Rtcutoff_dist_dom2
     &                             (amplitude_binindex2-1))
     &                             +Rtcutoff_dist_dom2
     &                             (amplitude_binindex2-1)
                           elseif(amplitude_binindex2.eq.1)then
                              Plm_Flx_RtRg=rand3(r3)*Rtcutoff_dist_dom2
     &                             (amplitude_binindex2)
                           endif
                           
                           if(rand5(r5).lt.0.5)then
                              Plm_Flx_RtRg=-Plm_Flx_RtRg
                           endif
                           
                           theta=Plm_Flx_RtRg*pai/180.0
                           phi=Plm_Flx_RtRg*pai/180.0
                           psai=Plm_Flx_RtRg*pai/180.0
                           
                           t(1,1)=cos(psai)*cos(phi)-
     &                          cos(theta)*sin(phi)*sin(psai)
                           t(1,2)=-sin(psai)*cos(phi)-
     &                          cos(theta)*sin(phi)*cos(psai)
                           t(1,3)=sin(theta)*sin(phi)
                           
                           t(2,1)=cos(psai)*sin(phi)+
     &                          cos(theta)*cos(phi)*sin(psai)
                           t(2,2)=-sin(psai)*sin(phi)+
     &                          cos(theta)*cos(phi)*cos(psai)
                           t(2,3)=-sin(theta)*cos(phi)
                           
                           t(3,1)=sin(psai)*sin(theta)
                           t(3,2)=cos(psai)*sin(theta)
                           t(3,3)=cos(theta)
                           
                           do k=2,Plm_res_num
                              Plm_B7_x_new(Plm_index,j,k)=t(1,1)*
     &                             (Plm_B7_x_new0(Plm_index,j,k)
     &                             -Plm_B7_x_new(Plm_index,j,1))+
     &                             t(1,2)*(Plm_B7_y_new0(Plm_index,j,k)
     &                             -Plm_B7_y_new(Plm_index,j,1))
     &                             +t(1,3)*(Plm_B7_z_new0(Plm_index,j,k)
     &                             -Plm_B7_z_new(Plm_index,j,1))
     &                             +Plm_B7_x_new(Plm_index,j,1)
                              Plm_B7_y_new(Plm_index,j,k)=t(2,1)*
     &                             (Plm_B7_x_new0(Plm_index,j,k)
     &                             -Plm_B7_x_new(Plm_index,j,1))+
     &                             t(2,2)*(Plm_B7_y_new0(Plm_index,j,k)
     &                             -Plm_B7_y_new(Plm_index,j,1))
     &                             +t(2,3)*(Plm_B7_z_new0(Plm_index,j,k)
     &                             -Plm_B7_z_new(Plm_index,j,1))
     &                             +Plm_B7_y_new(Plm_index,j,1)
                              Plm_B7_z_new(Plm_index,j,k)=t(3,1)*
     &                             (Plm_B7_x_new0(Plm_index,j,k)
     &                             -Plm_B7_x_new(Plm_index,j,1))+
     &                             t(3,2)*(Plm_B7_y_new0(Plm_index,j,k)
     &                             -Plm_B7_y_new(Plm_index,j,1))
     &                             +t(3,3)*(Plm_B7_z_new0(Plm_index,j,k)
     &                             -Plm_B7_z_new(Plm_index,j,1))
     &                             +Plm_B7_z_new(Plm_index,j,1)
                           enddo

                        endif
                        
                     enddo

ccc>>>>   Check Collisions
                     collision_flag=0
                     do i=1,Plm_tot_num
                        if(Plm_index.ne.i)then
c                           cm1_a_x=0
c                           cm1_a_y=0
c                           cm1_a_z=0
c                           do k=1,Plm_MHC_num
c                              cm1_a_x=cm1_a_x+Plm_MHC_x_new
c     &                             (Plm_index,k,1)
c                              cm1_a_y=cm1_a_y+Plm_MHC_y_new
c     &                             (Plm_index,k,1)
c                              cm1_a_z=cm1_a_z+Plm_MHC_z_new
c     &                             (Plm_index,k,1)
c                           enddo
c                           cm1_a_x=cm1_a_x/real(Plm_MHC_num)
c                           cm1_a_y=cm1_a_y/real(Plm_MHC_num)
c                           cm1_a_z=cm1_a_z/real(Plm_MHC_num)
c                           cm1_b_x=0
c                           cm1_b_y=0
c                           cm1_b_z=0
c                           do k=1,Plm_MHC_num
c                              cm1_b_x=cm1_b_x+Plm_MHC_x_new
c     &                             (i,k,1)
c                              cm1_b_y=cm1_b_y+Plm_MHC_y_new
c     &                             (i,k,1)
c                              cm1_b_z=cm1_b_z+Plm_MHC_z_new
c     &                             (i,k,1)
c                           enddo
c                           cm1_b_x=cm1_b_x/real(Plm_MHC_num)
c                           cm1_b_y=cm1_b_y/real(Plm_MHC_num)
c                           cm1_b_z=cm1_b_z/real(Plm_MHC_num)
c                           dist=sqrt((cm1_b_x-cm1_a_x)**2+
c     &                          (cm1_b_y-cm1_a_y)**2+
c     &                          (cm1_b_x-cm1_a_z)**2)
c                           if(dist.le.
c     &                          (2*Plm_radius))then
c                              collision_flag=1
c                           endif

                           do j=1,Plm_MHC_num
                              do k=1,Plm_MHC_num
                                 dist=sqrt((Plm_MHC_x_new(i,j,1)-
     &                                Plm_MHC_x_new
     &                                (Plm_index,k,1))**2+
     &                                (Plm_MHC_y_new(i,j,1)-
     &                                Plm_MHC_y_new
     &                                (Plm_index,j,1))**2+
     &                                (Plm_MHC_z_new(i,j,1)-
     &                                Plm_MHC_z_new
     &                                (Plm_index,j,1))**2)
                                 if(dist.lt.(MHC_radius+MHC_radius))then
                                    collision_flag=1
                                 endif
                              enddo
                           enddo
                           do j=1,Plm_MHC_num
                              do k=1,Plm_B7_num
                                 dist=sqrt((Plm_B7_x_new(i,j,1)-
     &                                Plm_MHC_x_new
     &                                (Plm_index,k,1))**2+
     &                                (Plm_B7_y_new(i,j,1)-
     &                                Plm_MHC_y_new
     &                                (Plm_index,j,1))**2+
     &                                (Plm_B7_z_new(i,j,1)-
     &                                Plm_MHC_z_new
     &                                (Plm_index,j,1))**2)
                                 if(dist.lt.(B7_radius+MHC_radius))then
                                    collision_flag=1
                                 endif
                              enddo
                           enddo
                           do j=1,Plm_B7_num
                              do k=1,Plm_B7_num
                                 dist=sqrt((Plm_B7_x_new(i,j,1)-
     &                                Plm_B7_x_new
     &                                (Plm_index,k,1))**2+
     &                                (Plm_B7_y_new(i,j,1)-
     &                                Plm_B7_y_new
     &                                (Plm_index,j,1))**2+
     &                                (Plm_B7_z_new(i,j,1)-
     &                                Plm_B7_z_new
     &                                (Plm_index,j,1))**2)
                                 if(dist.lt.(B7_radius+B7_radius))then
                                    collision_flag=1
                                 endif
                              enddo
                           enddo
                           do j=1,Plm_B7_num
                              do k=1,Plm_MHC_num
                                 dist=sqrt((Plm_MHC_x_new(i,j,1)-
     &                                Plm_B7_x_new
     &                                (Plm_index,k,1))**2+
     &                                (Plm_MHC_y_new(i,j,1)-
     &                                Plm_B7_y_new
     &                                (Plm_index,j,1))**2+
     &                                (Plm_MHC_z_new(i,j,1)-
     &                                Plm_B7_z_new
     &                                (Plm_index,j,1))**2)
                                 if(dist.lt.(MHC_radius+B7_radius))then
                                    collision_flag=1
                                 endif
                              enddo
                           enddo
                        endif
                     enddo

                     do i=1,TCR_tot_num
                        do j=1,Plm_MHC_num
                           dist=sqrt((TCR_x_new(i,1)-
     &                          Plm_MHC_x_new(Plm_index,j,1))**2+
     &                          (TCR_y_new(i,1)-
     &                          Plm_MHC_y_new(Plm_index,j,1))**2+
     &                          (TCR_z_new(i,1)-
     &                          Plm_MHC_z_new(Plm_index,j,1))**2)
                           if(dist.lt.(TCR_radius+MHC_radius)
     &                          )then
                              collision_flag=1
                           endif
                        enddo
                        do j=1,Plm_B7_num
                           dist=sqrt((TCR_x_new(i,1)-
     &                          Plm_B7_x_new(Plm_index,j,1))**2+
     &                          (TCR_y_new(i,1)-
     &                          Plm_B7_y_new(Plm_index,j,1))**2+
     &                          (TCR_z_new(i,1)-
     &                          Plm_B7_z_new(Plm_index,j,1))**2)
                           if(dist.lt.(TCR_radius+B7_radius)
     &                          )then
                              collision_flag=1
                           endif
                        enddo
                     enddo

                     do i=1,CD_tot_num
                        do j=1,Plm_MHC_num
                           dist=sqrt((CD_x_new(i,1)-
     &                          Plm_MHC_x_new(Plm_index,j,1))**2+
     &                          (CD_y_new(i,1)-
     &                          Plm_MHC_y_new(Plm_index,j,1))**2+
     &                          (CD_z_new(i,1)-
     &                          Plm_MHC_z_new(Plm_index,j,1))**2)
                           if(dist.lt.(CD_radius+MHC_radius)
     &                          )then
                              collision_flag=1
                           endif
                        enddo
                        do j=1,Plm_B7_num
                           dist=sqrt((CD_x_new(i,1)-
     &                          Plm_B7_x_new(Plm_index,j,1))**2+
     &                          (CD_y_new(i,1)-
     &                          Plm_B7_y_new(Plm_index,j,1))**2+
     &                          (CD_z_new(i,1)-
     &                          Plm_B7_z_new(Plm_index,j,1))**2)
                           if(dist.lt.(CD_radius+B7_radius)
     &                          )then
                              collision_flag=1
                           endif
                        enddo
                     enddo

                     do i=1,Plm_MHC_num
                        if(Plm_MHC_status_new(Plm_index,i).eq.1)then
                           do j=1,TCR_tot_num
                              if(Plm_MHC_idx_new(Plm_index,i).ne.j)then
                                 dist=sqrt((TCR_x_new(Plm_MHC_idx_new
     &                                (Plm_index,i),1)-
     &                                TCR_x_new(j,1))**2+
     &                                (TCR_y_new(Plm_MHC_idx_new
     &                                (Plm_index,i),1)-
     &                                TCR_y_new(j,1))**2+
     &                                (TCR_z_new(Plm_MHC_idx_new
     &                                (Plm_index,i),1)-
     &                                TCR_z_new(j,1))**2)
                                 if(dist.lt.(TCR_radius+TCR_radius)
     &                                )then
                                    collision_flag=1
                                 endif
                              endif
                           enddo
                           do j=1,CD_tot_num
                              dist=sqrt((TCR_x_new(Plm_MHC_idx_new
     &                             (Plm_index,i),1)-
     &                             CD_x_new(j,1))**2+
     &                             (TCR_y_new(Plm_MHC_idx_new
     &                             (Plm_index,i),1)-
     &                             CD_y_new(j,1))**2+
     &                             (TCR_z_new(Plm_MHC_idx_new
     &                             (Plm_index,i),1)-
     &                             CD_z_new(j,1))**2)
                              if(dist.lt.(TCR_radius+CD_radius)
     &                             )then
                                 collision_flag=1
                              endif
                           enddo
                        endif
                     enddo


                     do i=1,Plm_B7_num
                        if(Plm_B7_status_new(Plm_index,i).eq.1)then
                           do j=1,TCR_tot_num
                              dist=sqrt((CD_x_new(Plm_B7_idx_new
     &                             (Plm_index,i),1)-
     &                             TCR_x_new(j,1))**2+
     &                             (CD_y_new(Plm_B7_idx_new
     &                             (Plm_index,i),1)-
     &                             TCR_y_new(j,1))**2+
     &                             (CD_z_new(Plm_B7_idx_new
     &                             (Plm_index,i),1)-
     &                             TCR_z_new(j,1))**2)
                              if(dist.lt.(TCR_radius+CD_radius)
     &                             )then
                                 collision_flag=1
                              endif
                           enddo
                           do j=1,CD_tot_num
                              if(Plm_B7_idx_new(Plm_index,i).ne.j)then
                                 dist=sqrt((CD_x_new(Plm_B7_idx_new
     &                                (Plm_index,i),1)-
     &                                CD_x_new(j,1))**2+
     &                                (CD_y_new(Plm_B7_idx_new
     &                                (Plm_index,i),1)-
     &                                CD_y_new(j,1))**2+
     &                                (CD_z_new(Plm_B7_idx_new
     &                                (Plm_index,i),1)-
     &                                CD_z_new(j,1))**2)
                                 if(dist.lt.(CD_radius+CD_radius)
     &                                )then
                                    collision_flag=1
                                 endif
                              endif
                           enddo
                        endif
                     enddo

                     if(collision_flag.eq.1)then
                        do i=1,Plm_MHC_num
                           do j=1,Plm_res_num
                              Plm_MHC_x_new(Plm_index,i,j)=
     &                             Plm_MHC_x(Plm_index,i,j)
                              Plm_MHC_y_new(Plm_index,i,j)=
     &                             Plm_MHC_y(Plm_index,i,j)
                              Plm_MHC_z_new(Plm_index,i,j)=
     &                             Plm_MHC_z(Plm_index,i,j)
                              Plm_MHC_xo_new(Plm_index,i,j)=
     &                             Plm_MHC_xo(Plm_index,i,j)
                              Plm_MHC_yo_new(Plm_index,i,j)=
     &                             Plm_MHC_yo(Plm_index,i,j)
                              Plm_MHC_zo_new(Plm_index,i,j)=
     &                             Plm_MHC_zo(Plm_index,i,j)
                           enddo
                           if(Plm_MHC_status_new(Plm_index,i).eq.1)then
                              do j=1,TCR_res_num
                                 TCR_x_new(Plm_MHC_idx_new
     &                                (Plm_index,i),j)=
     &                                TCR_x(Plm_MHC_idx_new
     &                                (Plm_index,i),j)
                                 TCR_y_new(Plm_MHC_idx_new
     &                                (Plm_index,i),j)=
     &                                TCR_y(Plm_MHC_idx_new
     &                                (Plm_index,i),j)
                                 TCR_z_new(Plm_MHC_idx_new
     &                                (Plm_index,i),j)=
     &                                TCR_z(Plm_MHC_idx_new
     &                                (Plm_index,i),j)
                              enddo
                           endif
                        enddo
                        do i=1,Plm_B7_num
                           do j=1,Plm_res_num
                              Plm_B7_x_new(Plm_index,i,j)=
     &                             Plm_B7_x(Plm_index,i,j)
                              Plm_B7_y_new(Plm_index,i,j)=
     &                             Plm_B7_y(Plm_index,i,j)
                              Plm_B7_z_new(Plm_index,i,j)=
     &                             Plm_B7_z(Plm_index,i,j)
                              Plm_B7_xo_new(Plm_index,i,j)=
     &                             Plm_B7_xo(Plm_index,i,j)
                              Plm_B7_yo_new(Plm_index,i,j)=
     &                             Plm_B7_yo(Plm_index,i,j)
                              Plm_B7_zo_new(Plm_index,i,j)=
     &                             Plm_B7_zo(Plm_index,i,j)
                           enddo
                           if(Plm_B7_status_new(Plm_index,i).eq.1)then
                              do j=1,CD_res_num
                                 CD_x_new(Plm_B7_idx_new
     &                                (Plm_index,i),j)=
     &                                CD_x(Plm_B7_idx_new
     &                                (Plm_index,i),j)
                                 CD_y_new(Plm_B7_idx_new
     &                                (Plm_index,i),j)=
     &                                CD_y(Plm_B7_idx_new
     &                                (Plm_index,i),j)
                                 CD_z_new(Plm_B7_idx_new
     &                                (Plm_index,i),j)=
     &                                CD_z(Plm_B7_idx_new
     &                                (Plm_index,i),j)
                              enddo
                           endif
                        enddo

                     endif



                  endif 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   if the selected molecule is an TCR on 2D

               elseif((selecting_mole_index.gt.Plm_tot_num).AND.
     &                 (selecting_mole_index.le.
     &                 (Plm_tot_num+TCR_tot_num)))then

c                  TCR_index=int(rand3(r3)*TCR_tot_num)+1
                  TCR_index=selecting_mole_index-Plm_tot_num
cc>>>>   make random diffusion for selected molecule
cc>>>>   if it is unbound
                  if(TCR_status_new(TCR_index).eq.0)then
                     Prob_Diff=(4*TCR_D*time_step)/(distance_step**2)
                     if(rand3(r3).lt.Prob_Diff)then
cc>>     if smaller than the diffusion probability, move alone
                        phai=rand3(r3)*2*pai
                        distance_amp=2*rand4(r4)*distance_step
                        do i=1,TCR_res_num
                           TCR_x_new0(TCR_index,i)=
     &                          TCR_x(TCR_index,i)+
     &                          distance_amp*cos(phai)
                           TCR_y_new0(TCR_index,i)=
     &                          TCR_y(TCR_index,i)+
     &                          distance_amp*sin(phai)
                           TCR_z_new0(TCR_index,i)=
     &                          TCR_z(TCR_index,i)
                        enddo
                        PB_x=cell_range_x*anint
     &                          (TCR_x_new0(TCR_index,1)/cell_range_x)
                        PB_y=cell_range_y*anint
     &                          (TCR_y_new0(TCR_index,1)/cell_range_y)
                        do i=1,TCR_res_num
                           TCR_x_new0(TCR_index,i)=
     &                          TCR_x_new0(TCR_index,i)-PB_x     
                           TCR_y_new0(TCR_index,i)=
     &                          TCR_y_new0(TCR_index,i)-PB_y
                        enddo
                     else
                        do i=1,TCR_res_num
                           TCR_x_new0(TCR_index,i)=
     &                          TCR_x(TCR_index,i)
                           TCR_y_new0(TCR_index,i)=
     &                          TCR_y(TCR_index,i)
                           TCR_z_new0(TCR_index,i)=
     &                          TCR_z(TCR_index,i)
                        enddo
                     endif
c>>>>>>>>> Rotate the molecule

                     TCR_x_new(TCR_index,1)=TCR_x_new0(TCR_index,1)
                     TCR_y_new(TCR_index,1)=TCR_y_new0(TCR_index,1)
                     TCR_z_new(TCR_index,1)=TCR_z_new0(TCR_index,1)

                     TCR_x_new(TCR_index,TCR_res_num)=
     &                    TCR_x_new0(TCR_index,TCR_res_num)
                     TCR_y_new(TCR_index,TCR_res_num)=
     &                    TCR_y_new0(TCR_index,TCR_res_num)
                     TCR_z_new(TCR_index,TCR_res_num)=
     &                    TCR_z_new0(TCR_index,TCR_res_num)

                     rot_angle=2*rand3(r3)*TCR_rot_D-TCR_rot_D
                     x_axis_bg=TCR_x_new(TCR_index,1)
                     y_axis_bg=TCR_y_new(TCR_index,1)
                     z_axis_bg=-100.0
                     x_axis_ed=TCR_x_new(TCR_index,TCR_res_num)
                     y_axis_ed=TCR_y_new(TCR_index,TCR_res_num)
                     z_axis_ed=100.0
                     do j=2,TCR_res_num-1
                        x_bf_rot=TCR_x_new0(TCR_index,j)
                        y_bf_rot=TCR_y_new0(TCR_index,j)
                        z_bf_rot=TCR_z_new0(TCR_index,j)
                        x_af_rot=0
                        y_af_rot=0
                        z_af_rot=0
                        call rot_along_axis(rot_angle,
     &                       x_axis_bg,y_axis_bg,z_axis_bg,
     &                       x_axis_ed,y_axis_ed,z_axis_ed,
     &                       x_bf_rot,y_bf_rot,z_bf_rot,
     &                       x_af_rot,y_af_rot,z_af_rot)
                        TCR_x_new(TCR_index,j)=x_af_rot
                        TCR_y_new(TCR_index,j)=y_af_rot
                        TCR_z_new(TCR_index,j)=z_af_rot
                     enddo

ccc>>>>   Check Collisions

                     collision_flag=0
                     do i=1,TCR_tot_num
                        if(TCR_index.ne.i)then
                           dist=sqrt((TCR_x_new(i,1)-
     &                          TCR_x_new(TCR_index,1))**2+
     &                          (TCR_y_new(i,1)-
     &                          TCR_y_new(TCR_index,1))**2+
     &                          (TCR_z_new(i,1)-
     &                          TCR_z_new(TCR_index,1))**2)
                           if(dist.lt.(TCR_radius+TCR_radius))then
                              collision_flag=1
                           endif
                        endif
                     enddo

                     do i=1,CD_tot_num
                        dist=sqrt((CD_x_new(i,1)-
     &                       TCR_x_new(TCR_index,1))**2+
     &                       (CD_y_new(i,1)-
     &                       TCR_y_new(TCR_index,1))**2+
     &                       (CD_z_new(i,1)-
     &                       TCR_z_new(TCR_index,1))**2)
                        if(dist.lt.(TCR_radius+B7_radius)
     &                       )then
                           collision_flag=1
                        endif
                     enddo

                     do i=1,Plm_tot_num
                        do j=1,Plm_MHC_num
                           dist=sqrt((Plm_MHC_x_new(i,j,1)-
     &                          TCR_x_new(TCR_index,1))**2+
     &                          (Plm_MHC_y_new(i,j,1)-
     &                          TCR_y_new(TCR_index,1))**2+
     &                          (Plm_MHC_z_new(i,j,1)-
     &                          TCR_z_new(TCR_index,1))**2)
                           if(dist.lt.(MHC_radius+TCR_radius))then
                              collision_flag=1
                           endif
                        enddo
                        do j=1,Plm_B7_num
                           dist=sqrt((Plm_B7_x_new(i,j,1)-
     &                          TCR_x_new(TCR_index,1))**2+
     &                          (Plm_B7_y_new(i,j,1)-
     &                          TCR_y_new(TCR_index,1))**2+
     &                          (Plm_B7_z_new(i,j,1)-
     &                          TCR_z_new(TCR_index,1))**2)
                           if(dist.lt.(B7_radius+TCR_radius))then
                              collision_flag=1
                           endif
                        enddo
                     enddo

                     if(collision_flag.eq.1)then
                        do i=1,TCR_res_num
                           TCR_x_new(TCR_index,i)=
     &                          TCR_x(TCR_index,i)
                           TCR_y_new(TCR_index,i)=
     &                          TCR_y(TCR_index,i)
                           TCR_z_new(TCR_index,i)=
     &                          TCR_z(TCR_index,i)
                        enddo
                     endif




                  endif


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   if the selected molecule is an CD28 on 2D


               elseif((selecting_mole_index.gt.
     &                 (Plm_tot_num+TCR_tot_num)).AND.
     &                 (selecting_mole_index.le.
     &                 (Plm_tot_num+TCR_tot_num+CD_tot_num)))then


c                  CD_index=int(rand3(r3)*CD_tot_num)+1
                  CD_index=selecting_mole_index-Plm_tot_num-TCR_tot_num
cc>>>>   make random diffusion for selected molecule
cc>>>>   if it is unbound
                  if(CD_status_new(CD_index).eq.0)then
                     Prob_Diff=(4*CD_D*time_step)/(distance_step**2)
                     if(rand3(r3).lt.Prob_Diff)then
cc>>     if smaller than the diffusion probability, move alone
                        phai=rand3(r3)*2*pai
                        distance_amp=2*rand4(r4)*distance_step
                        do i=1,CD_res_num
                           CD_x_new0(CD_index,i)=
     &                          CD_x(CD_index,i)+
     &                          distance_amp*cos(phai)
                           CD_y_new0(CD_index,i)=
     &                          CD_y(CD_index,i)+
     &                          distance_amp*sin(phai)
                           CD_z_new0(CD_index,i)=
     &                          CD_z(CD_index,i)
                        enddo
                        PB_x=cell_range_x*anint
     &                          (CD_x_new0(CD_index,1)/cell_range_x)
                        PB_y=cell_range_y*anint
     &                          (CD_y_new0(CD_index,1)/cell_range_y)
                        do i=1,CD_res_num
                           CD_x_new0(CD_index,i)=
     &                          CD_x_new0(CD_index,i)-PB_x     
                           CD_y_new0(CD_index,i)=
     &                          CD_y_new0(CD_index,i)-PB_y
                        enddo
                     else
                        do i=1,CD_res_num
                           CD_x_new0(CD_index,i)=
     &                          CD_x(CD_index,i)
                           CD_y_new0(CD_index,i)=
     &                          CD_y(CD_index,i)
                           CD_z_new0(CD_index,i)=
     &                          CD_z(CD_index,i)
                        enddo
                     endif

c>>>>>>>>> Rotate the molecule

                     CD_x_new(CD_index,1)=CD_x_new0(CD_index,1)
                     CD_y_new(CD_index,1)=CD_y_new0(CD_index,1)
                     CD_z_new(CD_index,1)=CD_z_new0(CD_index,1)

                     CD_x_new(CD_index,CD_res_num)=
     &                    CD_x_new0(CD_index,CD_res_num)
                     CD_y_new(CD_index,CD_res_num)=
     &                    CD_y_new0(CD_index,CD_res_num)
                     CD_z_new(CD_index,CD_res_num)=
     &                    CD_z_new0(CD_index,CD_res_num)

                     rot_angle=2*rand3(r3)*CD_rot_D-CD_rot_D
                     x_axis_bg=CD_x_new(CD_index,1)
                     y_axis_bg=CD_y_new(CD_index,1)
                     z_axis_bg=-100.0
                     x_axis_ed=CD_x_new(CD_index,CD_res_num)
                     y_axis_ed=CD_y_new(CD_index,CD_res_num)
                     z_axis_ed=100.0
                     do j=2,TCR_res_num-1
                        x_bf_rot=CD_x_new0(CD_index,j)
                        y_bf_rot=CD_y_new0(CD_index,j)
                        z_bf_rot=CD_z_new0(CD_index,j)
                        x_af_rot=0
                        y_af_rot=0
                        z_af_rot=0
                        call rot_along_axis(rot_angle,
     &                       x_axis_bg,y_axis_bg,z_axis_bg,
     &                       x_axis_ed,y_axis_ed,z_axis_ed,
     &                       x_bf_rot,y_bf_rot,z_bf_rot,
     &                       x_af_rot,y_af_rot,z_af_rot)
                        CD_x_new(CD_index,j)=x_af_rot
                        CD_y_new(CD_index,j)=y_af_rot
                        CD_z_new(CD_index,j)=z_af_rot
                     enddo

ccc>>>>   Check Collisions

                     collision_flag=0
                     do i=1,TCR_tot_num
                        dist=sqrt((TCR_x_new(i,1)-
     &                       CD_x_new(CD_index,1))**2+
     &                       (TCR_y_new(i,1)-
     &                       CD_y_new(CD_index,1))**2+
     &                       (TCR_z_new(i,1)-
     &                       CD_z_new(CD_index,1))**2)
                        if(dist.lt.(TCR_radius+CD_radius))then
                           collision_flag=1
                        endif
                     enddo
 
                     do i=1,CD_tot_num
                        if(CD_index.ne.i)then
                           dist=sqrt((CD_x_new(i,1)-
     &                          CD_x_new(CD_index,1))**2+
     &                          (CD_y_new(i,1)-
     &                          CD_y_new(CD_index,1))**2+
     &                          (CD_z_new(i,1)-
     &                          CD_z_new(CD_index,1))**2)
                           if(dist.lt.(CD_radius+CD_radius))then
                              collision_flag=1
                           endif
                        endif
                     enddo

                     do i=1,Plm_tot_num
                        do j=1,Plm_MHC_num
                           dist=sqrt((Plm_MHC_x_new(i,j,1)-
     &                          CD_x_new(CD_index,1))**2+
     &                          (Plm_MHC_y_new(i,j,1)-
     &                          CD_y_new(CD_index,1))**2+
     &                          (Plm_MHC_z_new(i,j,1)-
     &                          CD_z_new(CD_index,1))**2)
                           if(dist.lt.(MHC_radius+TCR_radius))then
                              collision_flag=1
                           endif
                        enddo
                        do j=1,Plm_B7_num
                           dist=sqrt((Plm_B7_x_new(i,j,1)-
     &                          CD_x_new(CD_index,1))**2+
     &                          (Plm_B7_y_new(i,j,1)-
     &                          CD_y_new(CD_index,1))**2+
     &                          (Plm_B7_z_new(i,j,1)-
     &                          CD_z_new(CD_index,1))**2)
                           if(dist.lt.(B7_radius+TCR_radius))then
                              collision_flag=1
                           endif
                        enddo
                     enddo

                     if(collision_flag.eq.1)then
                        do i=1,CD_res_num
                           CD_x_new(CD_index,i)=
     &                          CD_x(CD_index,i)
                           CD_y_new(CD_index,i)=
     &                          CD_y(CD_index,i)
                           CD_z_new(CD_index,i)=
     &                          CD_z(CD_index,i)
                        enddo
                     endif


                  endif


               endif

            enddo


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   endo of diffusion for molecules
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   perform the reaction for molecules
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


ccc  1)  associate into a complex
            do i=1,Plm_tot_num
c>>   Between TCR and MHC
               do j=1,Plm_MHC_num
                  if(Plm_MHC_status_new(i,j).eq.0)then
                     do k=1,TCR_tot_num
                        if(TCR_status_new(k).eq.0)then
c>>>>>>>>>>>>>   calculate distance of two reaction sites between TCR and MHC
                           dist=sqrt((Plm_MHC_x_new(i,j,4)-
     &                          TCR_x_new(k,4))**2+
     &                          (Plm_MHC_y_new(i,j,4)-
     &                          TCR_y_new(k,4))**2+
     &                          (Plm_MHC_z_new(i,j,4)-
     &                          TCR_z_new(k,4))**2)
                           if(dist.lt.bond_dist_cutoff_TM)then
                         
c>>>>>>>>>>>> calculate two angle
c>> perpendicular angle
                              theta_pd=0
                              point_x(1)=TCR_x_new(k,1)-TCR_x_new(k,4)
                              point_y(1)=TCR_y_new(k,1)-TCR_y_new(k,4)
                              point_z(1)=TCR_z_new(k,1)-TCR_z_new(k,4)
                              point_x(2)=0
                              point_y(2)=0
                              point_z(2)=0
                              point_x(3)=Plm_MHC_x_new(i,j,1)
     &                             -Plm_MHC_x_new(i,j,4)
                              point_y(3)=Plm_MHC_y_new(i,j,1)
     &                             -Plm_MHC_y_new(i,j,4)
                              point_z(3)=Plm_MHC_z_new(i,j,1)
     &                             -Plm_MHC_z_new(i,j,4)
                              call gettheta(point_x,point_y,point_z,
     &                             theta_pd)
c>>   orientational angle
                              theta_ot=0
                              point_x(1)=TCR_x_new(k,1)-TCR_x_new(k,2)
                              point_y(1)=TCR_y_new(k,1)-TCR_y_new(k,2)
                              point_z(1)=TCR_z_new(k,1)-TCR_z_new(k,2)
                              point_x(2)=0
                              point_y(2)=0
                              point_z(2)=0
                              point_x(3)=Plm_MHC_x_new(i,j,1)
     &                             -Plm_MHC_x_new(i,j,2)
                              point_y(3)=Plm_MHC_y_new(i,j,1)
     &                             -Plm_MHC_y_new(i,j,2)
                              point_z(3)=Plm_MHC_z_new(i,j,1)
     &                             -Plm_MHC_z_new(i,j,2)
                              call gettheta(point_x,point_y,point_z,
     &                             theta_ot)

                              if((abs(theta_ot-bond_thetaot_TM)
     &                             .lt.bond_thetaot_cutoff_TM).AND.
     &                             (abs(theta_pd-bond_thetapd_TM)
     &                             .lt.bond_thetapd_cutoff_TM))then
                                 Prob_Ass=Ass_Rate_TM*time_step
     &                                *AssCalibration ! ready to change
                                 prob=rand3(r3)
                                 if((prob.lt.Prob_Ass).AND.
     &                                (prob.gt.0.000))then
                                    TCR_status_new(k)=1
                                    Plm_MHC_status_new(i,j)=1
                                    Plm_MHC_idx_new(i,j)=k
                                 endif
                              endif
                           endif
                        endif
                     enddo
                  endif
               enddo
c>>   Between CD28 and B7
               do j=1,Plm_B7_num
                  if(Plm_B7_status_new(i,j).eq.0)then
                     do k=1,CD_tot_num
                        if(CD_status_new(k).eq.0)then
c>>>>>>>>>>>>>   calculate distance of two reaction sites between TCR and MHC
                           dist=sqrt((Plm_B7_x_new(i,j,4)-
     &                          CD_x_new(k,4))**2+
     &                          (Plm_B7_y_new(i,j,4)-
     &                          CD_y_new(k,4))**2+
     &                          (Plm_B7_z_new(i,j,4)-
     &                          CD_z_new(k,4))**2)
                           if(dist.lt.bond_dist_cutoff_CB)then
                         
c>>>>>>>>>>>> calculate two angle
c>> perpendicular angle
                              theta_pd=0
                              point_x(1)=CD_x_new(k,1)-CD_x_new(k,4)
                              point_y(1)=CD_y_new(k,1)-CD_y_new(k,4)
                              point_z(1)=CD_z_new(k,1)-CD_z_new(k,4)
                              point_x(2)=0
                              point_y(2)=0
                              point_z(2)=0
                              point_x(3)=Plm_B7_x_new(i,j,1)
     &                             -Plm_B7_x_new(i,j,4)
                              point_y(3)=Plm_B7_y_new(i,j,1)
     &                             -Plm_B7_y_new(i,j,4)
                              point_z(3)=Plm_B7_z_new(i,j,1)
     &                             -Plm_B7_z_new(i,j,4)
                              call gettheta(point_x,point_y,point_z,
     &                             theta_pd)
c>>   orientational angle
                              theta_ot=0
                              point_x(1)=CD_x_new(k,1)-CD_x_new(k,2)
                              point_y(1)=CD_y_new(k,1)-CD_y_new(k,2)
                              point_z(1)=CD_z_new(k,1)-CD_z_new(k,2)
                              point_x(2)=0
                              point_y(2)=0
                              point_z(2)=0
                              point_x(3)=Plm_B7_x_new(i,j,1)
     &                             -Plm_B7_x_new(i,j,2)
                              point_y(3)=Plm_B7_y_new(i,j,1)
     &                             -Plm_B7_y_new(i,j,2)
                              point_z(3)=Plm_B7_z_new(i,j,1)
     &                             -Plm_B7_z_new(i,j,2)
                              call gettheta(point_x,point_y,point_z,
     &                             theta_ot)

                              if((abs(theta_ot-bond_thetaot_CB)
     &                             .lt.bond_thetaot_cutoff_CB).AND.
     &                             (abs(theta_pd-bond_thetapd_CB)
     &                             .lt.bond_thetapd_cutoff_CB))then
                                 Prob_Ass=Ass_Rate_CB*time_step
     &                                *AssCalibration ! ready to change
                                 prob=rand3(r3)
                                 if((prob.lt.Prob_Ass).AND.
     &                                (prob.gt.0.000))then
                                    CD_status_new(k)=1
                                    Plm_B7_status_new(i,j)=1
                                    Plm_B7_idx_new(i,j)=k
                                 endif
                              endif
                           endif
                        endif
                     enddo
                  endif
               enddo
            enddo


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc  2)  complex dissociate 

            do i=1,Plm_tot_num
               do j=1,Plm_MHC_num
                  if(Plm_MHC_status_new(i,j).eq.1)then
                     do k=1,TCR_tot_num
                        if((Plm_MHC_idx_new(i,j).eq.k)
     &                       .AND.(TCR_status_new(k).eq.1))then
                           selected_TCR=k
                        endif
                     enddo
                     Prob_Diss=Diss_Rate_TM*time_step
                     part1=int(rand3(r3)*10)
                     part2=int(rand5(r5)*10)
                     part3=int(rand5(r5)*10)
                     temp=rand3(r3)
                     part4=int(rand3(r3)*10)
                     part5=int(rand4(r4)*10)
                     do j2=1,int(rand4(r4)*10)
                        temp=rand5(r5)
                        temp=rand5(r5)
                     enddo
                     part6=int(rand5(r5)*10)
                     part7=int(rand5(r5)*10)
                     part8=int(rand4(r4)*10)
                     part9=int(rand5(r5)*10)                     
                     prob=real(part1)/10.0+real(part2)/100.0
     &                    +real(part3)/1000.0+real(part4)/10000.0
     &                    +real(part5)/100000.0
     &                    +real(part6)/1000000.0
     &                    +real(part7)/10000000.0
     &                    +real(part8)/100000000.0
     &                    +real(part9)/1000000000.0
                     if((prob.lt.Prob_Diss).AND.(prob.gt.0.000))then
                        TCR_status_new(selected_TCR)=0
                        Plm_MHC_status_new(i,j)=0
                        Plm_MHC_idx_new(i,j)=0
                     endif
                  endif
               enddo
               do j=1,Plm_B7_num
                  if(Plm_B7_status_new(i,j).eq.1)then
                     do k=1,CD_tot_num
                        if((Plm_B7_idx_new(i,j).eq.k)
     &                       .AND.(CD_status_new(k).eq.1))then
                           selected_CD=k
                        endif
                     enddo
                     Prob_Diss=Diss_Rate_CB*time_step
                     part1=int(rand3(r3)*10)
                     part2=int(rand5(r5)*10)
                     part3=int(rand5(r5)*10)
                     temp=rand3(r3)
                     part4=int(rand3(r3)*10)
                     part5=int(rand4(r4)*10)
                     do j2=1,int(rand4(r4)*10)
                        temp=rand5(r5)
                        temp=rand5(r5)
                     enddo
                     part6=int(rand5(r5)*10)
                     part7=int(rand5(r5)*10)
                     part8=int(rand4(r4)*10)
                     part9=int(rand5(r5)*10)                     
                     prob=real(part1)/10.0+real(part2)/100.0
     &                    +real(part3)/1000.0+real(part4)/10000.0
     &                    +real(part5)/100000.0
     &                    +real(part6)/1000000.0
     &                    +real(part7)/10000000.0
     &                    +real(part8)/100000000.0
     &                    +real(part9)/1000000000.0
                     if((prob.lt.Prob_Diss).AND.(prob.gt.0.000))then
                        CD_status_new(selected_CD)=0
                        Plm_B7_status_new(i,j)=0
                        Plm_B7_idx_new(i,j)=0
                     endif
                  endif
               enddo

            enddo

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc         simulation update
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   update the coordinates for all molecules

            do i=1,TCR_tot_num
               do j=1,TCR_res_num
                  TCR_x(i,j)=TCR_x_new(i,j)
                  TCR_y(i,j)=TCR_y_new(i,j)
                  TCR_z(i,j)=TCR_z_new(i,j)
               enddo
               TCR_status(i)=TCR_status_new(i)
            enddo
            do i=1,CD_tot_num
               do j=1,CD_res_num
                  CD_x(i,j)=CD_x_new(i,j)
                  CD_y(i,j)=CD_y_new(i,j)
                  CD_z(i,j)=CD_z_new(i,j)
               enddo
               CD_status(i)=CD_status_new(i)
            enddo

            do i=1,Plm_tot_num
               do j=1,Plm_MHC_num
                  do k=1,Plm_res_num
                     Plm_MHC_x(i,j,k)=Plm_MHC_x_new(i,j,k)
                     Plm_MHC_y(i,j,k)=Plm_MHC_y_new(i,j,k)
                     Plm_MHC_z(i,j,k)=Plm_MHC_z_new(i,j,k)
                     Plm_MHC_xo(i,j,k)=Plm_MHC_xo_new(i,j,k)
                     Plm_MHC_yo(i,j,k)=Plm_MHC_yo_new(i,j,k)
                     Plm_MHC_zo(i,j,k)=Plm_MHC_zo_new(i,j,k)
                  enddo
                  Plm_MHC_status(i,j)=Plm_MHC_status_new(i,j)
                  Plm_MHC_idx(i,j)=Plm_MHC_idx_new(i,j)
               enddo
               do j=1,Plm_B7_num
                  do k=1,Plm_res_num
                     Plm_B7_x(i,j,k)=Plm_B7_x_new(i,j,k)
                     Plm_B7_y(i,j,k)=Plm_B7_y_new(i,j,k)
                     Plm_B7_z(i,j,k)=Plm_B7_z_new(i,j,k)
                     Plm_B7_xo(i,j,k)=Plm_B7_xo_new(i,j,k)
                     Plm_B7_yo(i,j,k)=Plm_B7_yo_new(i,j,k)
                     Plm_B7_zo(i,j,k)=Plm_B7_zo_new(i,j,k)
                  enddo
                  Plm_B7_status(i,j)=Plm_B7_status_new(i,j)
                  Plm_B7_idx(i,j)=Plm_B7_idx_new(i,j)
               enddo
            enddo

            TCR_complex_num=0
            do i=1,TCR_tot_num
               if(TCR_status(i).eq.1)then
                  TCR_complex_num=TCR_complex_num+1
c                  print*,'TCR',i
                  do j=1,Plm_tot_num
                     do k=1,Plm_MHC_num
                        if(Plm_MHC_idx(j,k).eq.i)then
c                           print*,'MHC',j,k
                        endif
                     enddo
                  enddo
               endif
            enddo
            CD_complex_num=0
            do i=1,CD_tot_num
               if(CD_status(i).eq.1)then
                  CD_complex_num=CD_complex_num+1
c                  print*,'CD28',i
               endif
            enddo

            DoubleBind_num=0
            do i=1,Plm_tot_num
               MHC_flag=0
               B7_flag=0
               do j=1,Plm_MHC_num
                  if(Plm_MHC_status(i,j).eq.1)then
                     MHC_flag=1
                  endif
               enddo
               do j=1,Plm_B7_num
                  if(Plm_B7_status(i,j).eq.1)then
                     B7_flag=1
                  endif
               enddo
               if((MHC_flag.eq.1).and.(B7_flag.eq.1))then
                  DoubleBind_num=DoubleBind_num+1
               endif
            enddo

c            print*,mc_time_step,TCR_complex_num,CD_complex_num
c            print*,Plm_MHC_xo_new0(18,1,1),Plm_MHC_x_new0(18,1,1)
c            print*,Plm_MHC_yo_new0(18,1,1),Plm_MHC_y_new0(18,1,1)
ccccccccccccccccccccccccccccccccccccccccccccccccccc
c  output data along the trajectory
ccccccccccccccccccccccccccccccccccccccccccccccccccc           
c            goto 300
c            if((mc_time_step.ge.86000).AND.(mc_time_step.le.87000))then
c            print*,mc_time_step
            if((MOD(mc_time_step,100).eq.0).AND.
     &           (ene_output_flag.eq.1))then
               open (unit=10,file=
     &              'BispRBKMCresult_ene_'//serialnum//'.dat',
     &              status='unknown',access='append')
               write(10,3000) mc_time_step,
     &              TCR_complex_num,CD_complex_num,DoubleBind_num
               close(10)
            endif
 3000       format(I12,1x,I5,1x,I5,1x,I5)
            if((MOD(mc_time_step,10000).eq.0).AND.
     &           (trj_output_flag.eq.1))then
               open (unit=10,file=
     &              'BispRBKMCresult_trj_'//serialnum//'.pdb',
     &              status='unknown',access='append')
               do i=1,TCR_tot_num
                  write(10,2100) 'ATOM  ',i,' CA ','ALA'
     &                 ,i,TCR_x(i,1),TCR_y(i,1),TCR_z(i,1)
                  write(10,2100) 'ATOM  ',i,' CA ','VAL'
     &                 ,i,TCR_x(i,2),TCR_y(i,2),TCR_z(i,2)
                  write(10,2100) 'ATOM  ',i,' CA ','VAL'
     &                 ,i,TCR_x(i,3),TCR_y(i,3),TCR_z(i,3)
                  write(10,2100) 'ATOM  ',i,' CA ','VAL'
     &                 ,i,TCR_x(i,4),TCR_y(i,4),TCR_z(i,4)
               enddo
               write(10,2102) 'TER'    
               do i=1,CD_tot_num
                  write(10,2100) 'ATOM  ',i,' CA ','LYS'
     &                 ,i,CD_x(i,1),CD_y(i,1),CD_z(i,1)
                  write(10,2100) 'ATOM  ',i,' CA ','TYR'
     &                 ,i,CD_x(i,2),CD_y(i,2),CD_z(i,2)
                  write(10,2100) 'ATOM  ',i,' CA ','TYR'
     &                 ,i,CD_x(i,3),CD_y(i,3),CD_z(i,3)
                  write(10,2100) 'ATOM  ',i,' CA ','TYR'
     &                 ,i,CD_x(i,4),CD_y(i,4),CD_z(i,4)
               enddo
               write(10,2102) 'TER'    
               do i=1,Plm_tot_num
                  do j=1,Plm_MHC_num
                     write(10,2100) 'ATOM  ',i,' CA ','LEU'
     &                    ,i,Plm_MHC_x(i,j,1),
     &                    Plm_MHC_y(i,j,1),Plm_MHC_z(i,j,1)
                     write(10,2100) 'ATOM  ',i,' CA ','GLY'
     &                    ,i,Plm_MHC_x(i,j,2),
     &                    Plm_MHC_y(i,j,2),Plm_MHC_z(i,j,2)
                     write(10,2100) 'ATOM  ',i,' CA ','GLY'
     &                    ,i,Plm_MHC_x(i,j,3),
     &                    Plm_MHC_y(i,j,3),Plm_MHC_z(i,j,3)
                     write(10,2100) 'ATOM  ',i,' CA ','GLY'
     &                    ,i,Plm_MHC_x(i,j,4),
     &                    Plm_MHC_y(i,j,4),Plm_MHC_z(i,j,4)
                  enddo
                  write(10,2102) 'TER'  
                  do j=1,Plm_B7_num
                     write(10,2100) 'ATOM  ',i,' CA ','PRO'
     &                    ,i,Plm_B7_x(i,j,1),
     &                    Plm_B7_y(i,j,1),Plm_B7_z(i,j,1)
                     write(10,2100) 'ATOM  ',i,' CA ','PHE'
     &                    ,i,Plm_B7_x(i,j,2),
     &                    Plm_B7_y(i,j,2),Plm_B7_z(i,j,2)
                     write(10,2100) 'ATOM  ',i,' CA ','PHE'
     &                    ,i,Plm_B7_x(i,j,3),
     &                    Plm_B7_y(i,j,3),Plm_B7_z(i,j,3)
                     write(10,2100) 'ATOM  ',i,' CA ','PHE'
     &                    ,i,Plm_B7_x(i,j,4),
     &                    Plm_B7_y(i,j,4),Plm_B7_z(i,j,4)
                  enddo
                  write(10,2102) 'TER'    
               enddo
               write(10,2102) 'TER'    
               write(10,2102) 'END'    
               close(10)
            endif
c            endif
 2100       format(A6,I5,1x,A4,1x,A3,1x,I5,4x,3F8.3)
 2102       format(A3)          
 300        continue
ccccccccccccccccccccccccccccccccccc
cc  end current simulation step
ccccccccccccccccccccccccccccccccccc

         enddo

cccccccccccccccccccccccccccc
cc end current trajectory
cccccccccccccccccccccccccccc

      enddo

cccccccccccccccccccc

      stop
      end

ccccccccccccccccccccccccccccccccccccccccccccccccc

      real  function rand3(r3)
      double precision s,u,v,r3
      s=65536.0
      u=2053.0
      v=13849.0
      m=r3/s
      r3=r3-m*s
      r3=u*r3+v
      m=r3/s
      r3=r3-m*s
      rand3=r3/s
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccc            
ccccccccccccccccccccccccccccccccccccccccccccccccc

      real  function rand4(r4)
      double precision s,u,v,r4
      s=65536.0
      u=2053.0
      v=13849.0
      m=r4/s
      r4=r4-m*s
      r4=u*r4+v
      m=r4/s
      r4=r4-m*s
      rand4=r4/s
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccc                      
ccccccccccccccccccccccccccccccccccccccccccccccccc

      real  function rand5(r5)
      double precision s,u,v,r5
      s=65536.0
      u=2053.0
      v=13849.0
      m=r5/s
      r5=r5-m*s
      r5=u*r5+v
      m=r5/s
      r5=r5-m*s
      rand5=r5/s
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccc                                  
ccccccccccccccccccccccccccccccccccccccccccccccccc            
            
      subroutine gettheta(point_x,point_y,point_z,test_theta)

      implicit none
      real*8 point_x(3),point_y(3),point_z(3),test_theta
ccc
      real*8 lx(2),ly(2),lz(2),lr(2)
      real*8 conv,doth1,doth2,conv2
      integer in

ccccccccccccccccccc
c>>  bond lengt
ccccccccccccccccccc

      do in=1,2
         lx(in)=0
         ly(in)=0
         lz(in)=0
         lr(in)=0
      enddo

      lx(1)=point_x(2)-point_x(1)
      ly(1)=point_y(2)-point_y(1)
      lz(1)=point_z(2)-point_z(1)
      lr(1)=sqrt(lx(1)**2+ly(1)**2+lz(1)**2)
      
      lx(2)=point_x(3)-point_x(2)
      ly(2)=point_y(3)-point_y(2)
      lz(2)=point_z(3)-point_z(2)
      lr(2)=sqrt(lx(2)**2+ly(2)**2+lz(2)**2)
 
cccccccccccccccccccc
c>>  theta value
cccccccccccccccccccc
      
      test_theta=0

      conv=180/3.14159 
      do in=1,1
         doth1=-1*(lx(in+1)*lx(in)+ly(in+1)*ly(in)+lz(in+1)*lz(in))
         doth2=doth1/(lr(in+1)*lr(in))
         if(doth2.gt.1)then
            doth2=1
         endif
         if(doth2.lt.-1)then
            doth2=-1
         endif
         test_theta=acos(doth2)*conv
      enddo

cccccccccccccccc

      return
      end

cccccccccccccccccccccccccccccccccccc
                        
