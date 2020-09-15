      subroutine rot_along_axis(rot_angle,x_axis_bg,y_axis_bg,z_axis_bg,
     &     x_axis_ed,y_axis_ed,z_axis_ed,x_bf_rot,y_bf_rot,z_bf_rot,
     &     x_af_rot,y_af_rot,z_af_rot)
      implicit none
      real*8 rot_angle
      real*8 x_axis_bg,y_axis_bg,z_axis_bg
      real*8 x_axis_ed,y_axis_ed,z_axis_ed
      real*8 x_bf_rot,y_bf_rot,z_bf_rot
      real*8 x_af_rot,y_af_rot,z_af_rot

      real*8 theta
      real*8 x1,y1,z1
      real*8 x2,y2,z2
      real*8 x3,y3,z3
      real*8 x4,y4,z4
      real*8 x5,y5,z5      
      real*8 p,q,ss,mag_pqr
      real*8 u,v,w
      real*8 g,h
      real*8 f1,f2
      real*8 aa,b,c,d

cccccccccccccccccccccccccccc

      theta=rot_angle

      x1=x_axis_bg
      y1=y_axis_bg
      z1=z_axis_bg
      
      x2=x_axis_ed
      y2=y_axis_ed
      z2=z_axis_ed
      
      x4=x_bf_rot
      y4=y_bf_rot
      z4=z_bf_rot
      
cccccccccccc

      p=x2-x1
      q=y2-y1
      ss=z2-z1
      
      mag_pqr=sqrt(p**2+q**2+ss**2)

      u=q*(z4-z1)-ss*(y4-y1)
      v=ss*(x4-x1)-p*(z4-z1)
      w=p*(y4-y1)-q*(x4-x1)
      
      h=sqrt(u**2+v**2+w**2)/mag_pqr

      g=sqrt((x4-x1)**2+(y4-y1)**2+(z4-z1)**2)

      f1=sqrt(g**2-h**2)
      
      f2=sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)

      x3=f1*(x2-x1)/f2+x1
      y3=f1*(y2-y1)/f2+y1
      z3=f1*(z2-z1)/f2+z1

ccccccccccccccccccccccc
c coordinate transfer
ccccccccccccccccccccccc

      x1=x1-x3
      y1=y1-y3
      z1=z1-z3
      
      x2=x2-x3
      y2=y2-y3
      z2=z2-z3
      
      x4=x4-x3
      y4=y4-y3
      z4=z4-z3

cccccccccccccccccc

      d=h**2
              
      c=h**2/f2

      b=(h**2)*COS(theta*3.1415926/180)

      aa=SIN(theta*3.1415926/180)

cccccccccccccccccccc

      x5=(b*x4+aa*c*(z4*(y1-y2)-y4*(z1-z2)))/d+x3
      y5=(b*y4+aa*c*(x4*(z1-z2)-z4*(x1-x2)))/d+y3
      z5=(b*z4+aa*c*(y4*(x1-x2)-x4*(y1-y2)))/d+z3

cccccccccccccccccccccc

      x_af_rot=x5    
      y_af_rot=y5    
      z_af_rot=z5      

cccccccccccccccccccccccccccccccccccc

      return
      end
