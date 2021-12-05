*deck,userCZM     parallel       optimize  USERDISTRIB               gal
      subroutine userCZM (matId, elemId, kMatIntPt, ldstep,isubst,
     &                    keycut, ncomp,nProp, nstatev,
     &                    Time, dTime, Temp, dTemp,
     &                    coords, prop, Strain, dStrain, 
     &                    stress, dsdePl, sedEl, sedPl, statev,
     &                    var1, var2, var3, var4, var5)
c*************************************************************************
c     *** primary function ***
c
c           user cohesive zone model example
c
c           Commands
c             TB,CZM,mat,NTEMP,NPTS,user 
c                TBTEMP if mat. constants are temperature dependent
c                TBDATA define material constants
c*************************************************************************
c     input arguments
c     ===============
c      matId     (int,sc,in)              material #
c      elemId    (int,sc,in)              element #
c      kMatIntPt (int,sc,in)              material integration point #
c      ldstep    (int,sc,in)              load step number
c      isubst    (int,sc,in)              substep number
c      ncomp     (int,sc,in)              number of stress, strain components
c      nProp     (int,sc,in)              Number of material ocnstants
c      nstatev   (int,sc,in)              Number of state variables
c
c      Temp      (dp ,sc,in)              temperature at beginning of time increment
c      dTemp     (dp ,sc,in)              temperature increment 
c      Time      (dp ,sc,in)              time at beginning of increment (t)
c      dTime     (dp ,sc,in)              time increment (dt)
c
c      prop     (dp,ar(nprop),i)          Material constants defined by TB command 
c      Strain   (dp,ar(ncomp),i)          Interface separation at beginning of time increment
c      dStrain  (dp,ar(ncomp),i)          Interface separation increment
c      coords   (dp,ar(3),i)              current coordinates
c
c     output arguments              
c     ======================             
c      stress   (dp,ar(nTesn),io)         Traction stress
c      sedEl    (dp,sc,io)                elastic work
c      sedPl    (dp,sc,io)                plastic work
c      keycut   (int,sc,io)               loading bisect/cut control
c                                         0 - no bisect/cut
c                                         1 - bisect/cut 
c                                         (factor will be determined by ANSYS solution control)
c      dsdePl   (dp,ar(ncomp,ncomp),io)   consistent tangent jacobian matrix
c
c     input output arguments              
c     ======================             
c      statev   (dp,ar(nstatev,io)        user defined solution state variables
c
c     misc.
c     ======================             
c      var1, var2, var3, var4, var5       currently not used
c
c     local variables
c     ======================             
c
c      debugflag (in,sc, l)                debugflag to print debug information
c*****************************************************************************************************************************************************************************
c                                                 USER'S INDICATIONS
c MATERIAL CONSTANTS (imput data):
c In TBDATA there are 4+(2*n)+(2*t1)+(2*t2) material constants: 
c  -the first 4 props are: prop(1)=t_a, prop(2)=n-1, prop(3)=t1-1, prop(4)=t2-1
c  -after there are separations (d) and tractions (sig,tau) for each direction ; i.e  t_a,n,t1,t2,d_1,....,d_n,sig_1,....,sig_n,d_1,...d_t1,tau_1,...tau_t1,...etc.
c
c  prop(2) = number points of normal traction-separation law (modeI) (excluding first point in origin of axes)
c  prop(3) = number points of first tangential t-s law (modeII) (excluding first point in origin of axes)
c  prop(4) = number points of second tangential t-s law (modeIII) (excluding first point in origin of axes)
c
c  t_a        (int,sc) = "type analysis" determines the kind of interpolation of point's tractions-separation laws 
c             if t_a =1  Linear interpolation, otherwise Polynomial interpolation.The polynomial interpolation is calculated with Newton's method whith Divided Differences. 
c
c MATERIAL STATE VARIABLES (imput data):
c In TBDATA there are 4+(ncomp+1) material state variables:
c  -the first 4 state variable are internal Ansys state variable (not specified): kchstat=statev(1), kestat0=statev(2), kestat=statev(3), stiffi=statev(4).
c  -after there are user subroutine state variable:
c     if ncomp=2  state(5)=G_last_old(1), statev(6)=G_last_old(2), statev(7)=d_old.
c     if ncomp=3  state(5)=G_last_old(1), statev(6)=G_last_old(2), statev(7)=G_last_old(3), statev(8)=d_old.
c     Please note; for undamaged CZM material: G_last_old(i)=G0(i) for i=,ncomp, d_old=0.
c
c LOCAL USER SUBROUTINE VARIABLES
c  n          (int,sc) = prop(2)+1 number points of normal traction-separation law (modeI) (including first point in origin of axes)
c  t1         (int,sc) = prop(3)+1 number points of first tangential t-s law (modeII) (including first point in origin of axes)
c  t2         (int,sc) = prop(4)+1 number points of second tangential t-s law (modeIII) (including first point in origin of axes)
c 
c  v          (dp,ar(i,j,k)) = array that contains points coordinates of cohesive laws, it has 3 dimensional:
c             where: i=3, it is number of rows: i=1 there are separations of t-s laws; i=2 there are tractions of t-s laws; i=3 there are areas below segments of the t-s laws.
c                    j=m_p+1, it is number of columns refer to max number of points of t-s laws. 
c                    k=ncomp, it is number of sheets: k=1 for normal t-s law; k=2 for first tangential t-s law; k=3 for second tangential t-s law.
c  
c  v_mix      (dp,ar(3,m_p+1,3)) = array that contains points coordinates of cohesive laws, modified to calculate the cohesive law in a mixed-mode. 
c             In particular in the case of linear interpolation, the points are chosen at regular intervals, while in the case of polynomial interpolation the points are chosen 
c             with intervals given by the Čebyšëv-Gauss-Lobatto nodes, used to minimize the Runge phenomenon. The structure of the v_mix array is similar to that of v, except that
c             for all modes we have the same number of points equal to the number of the cohesive law with maximum number of interpolating points. 
c
c  czm_mix    (dp,ar(3,m_p+1)) = array that contain point's coordinates of cohesive law in mixed-mode.
c
c  T_inter    (dp,ar,(i,k)) = array that contains cohesive tractions, and local derivatives of cohesive laws when damage process accours.
c             where: i=2, it is number of rows: i=1 there are cohesive tractions, i=2 there are local derivatives of cohesive laws.  
c                    k=ncomp, it is number of columns: k=1 for normal t-s law; k=2 for first tangential t-s law; k=3 for second tangential t-s law. 
c
c  T_inter_mix(dp,ar,(i)) = array that contains mixed cohesive traction, and local derivative of mixed cohesive law.
c             where: i=2, i=1 there is mixed cohesive traction, i=2 there is local derivative of mixed cohesive laws.  
c
c  K0         (dp,ar(ncomp)) = array that contains initial penalty stiffness. 
c  d          (dp,ar(ncomp)) = array that contains damage parameters each modes, used for specific case of single mode actived (pure mode).
c  G0         (dp,ar(ncomp)) = array that contains total initial area below t-s laws i.e., critical strain energy release rate: GIc,GIIc,GIIIc.
c  G_last     (dp,ar(ncomp)) = array that contains remaining area under t-s laws curve after damage.
c  G          (dp,ar(3))     = array that contains dissipated strain energy during damage.   
c  G_last_old (dp,ar(ncomp)) = array that contains min historical remaining area below t-s laws (state variables).
c
c  delta_mix  (dp) = equivalent separation in mixed-mode
c  d_mix      (dp) = damage parameter mixed-mode
c  d_old      (dp) = max historical damage parameter (state variable).
c  K0_mix     (dp) = initial penalty stiffness mixed-mode 
c  GT         (dp) = G(1)+G(2)+G(3) 
c 
c                                                USER PARAMETER SETUP
c  m_p        (int,sc) ="max_points" is max number of interpolating points of t-s laws (excluding first point in origin of axes),thus m_p>=max(n,t1,t2)
c  i_a        (int,sc) ="integral_accuracy"; greater is this parameter best will be approssimation of G0 and G_last, calculated with "trapezium method",during polynomial interpolation. 
c  ETA        (dp)= is the esponent of energetic fracture criteria: Benzeggagh and Kenane criterion for 2D mixed_Mode problem, and Reeder criterion for 3D mixed_Mode problem
c  s_t        (logical)= "same t-s law" is a binary parameter used to distinguish the cases the cohesive laws of mode 2 and 3 are same or different: if s_t = true so are same, otherwise no.
c
c Subroutine's author: student Francesco Starace, Master's degree Mechanical Engineering for Design and Production, Federico II University of Naples, academic year 2020-2021 
c Author's contacts: fran.starace@protonmail.com       
c             
c*****************************************************************************************************************************************************************************
c
#include "impcom.inc"
c
      interface
        function S(V,delta) 
        double precision, dimension(:,:,:), intent(in) :: V
        double precision, dimension(3), intent(in) :: delta
        double precision, dimension(2,3) :: S
        endfunction
      endinterface
      interface
        function S_mix(czm_mix,delta_mix) 
        double precision, dimension(:,:), intent(in) :: czm_mix
        double precision, intent(in) :: delta_mix
        double precision, dimension(2) :: S_mix
        endfunction
      endinterface
      INTEGER          matId,elemId,kMatIntPt,ldstep,isubst,keycut,     
     &                 ncomp,nProp,nstatev,t_a,n,t1,t2,i,j,k,l,z,nco
c
c --- COMMON
c          
      COMMON /INDICES/ t_a,n,t1,t2,nco 
c
      DOUBLE PRECISION Time,dTime,Temp,dTemp,sedEl,sedPl,coords(3), 
     &                 prop(nProp), Strain(ncomp), dStrain(ncomp), 
     &                 stress(ncomp), dsdePl(ncomp,ncomp),
     &                 statev(nstatev), var1, var2, var3, var4, var5  
c
c --- PARAMETERS
c
      INTEGER          m_p,i_a
      DOUBLE PRECISION PI,ETA
      LOGICAL          s_t
      PARAMETER       (m_p=2,       !"max_points"
     &                 i_a=50,      !"integral_accuracy"
     &                 PI=3.14159d0,!PI
     &                 ETA=1.5d0,   !esponent of energetic fracture criteria B-K
     &                 s_t=.true.)  !"same t-s law" 
c
      DOUBLE PRECISION delta_mix,K0_mix,d_mix,d_old,d_mix_der,GT,
     &                 y_1,y_2,y_3,y_4,w_0,w_1,w_2,R
      DOUBLE PRECISION, DIMENSION(ncomp) :: k0,d,d_der,G0,G_last,  
     &                                      G_last_old,X,Y
      DOUBLE PRECISION, DIMENSION(3) :: delta,G,puremode,
     &                                  x_1,x_2,x_3,x_4
      DOUBLE PRECISION v(3,m_p+1,ncomp),v_mix(3,m_p+1,3),
     &                 czm_mix(3,m_p+1),T_inter(2,3),T_inter_mix(2) 
c
      INTEGER          kchstat, debugflag, kestat, kestat0, kwr     
      DOUBLE PRECISION stiffi, dsigCap, kdsigCap 
c
c --- debug includes
#include "locknm.inc"
      INTEGER          wrinqr, iott
      EXTERNAL         pplock, ppunlock, wrinqr
c
c*************************************************************************
c
!     INITIALIZATION

      debugflag = 1
      keycut   = 0
      kwr = 0
      sedPl    = 0.0
      stiffi   = 0.0
      sedEl    = 0.0
      sedPl    = 0.0
      keycut   = 0
      kchstat = nint(statev(1))
      kestat0 = nint(statev(2))
      kestat  = kestat0
      
      dsdePl=0.0
      nco=ncomp
      delta_mix=0.0
      k0_mix=0.0
      d_mix=0.0
      d_mix_der=0.0
      GT=0.0
      y_1=0.0
      y_2=0.0
      y_3=0.0
      y_4=0.0
      w_0=0.0
      w_1=0.0
      w_2=0.0
      R=0.0
      k0=0.0
      d=0.0
      d_der=0.0
      G0=0.0
      G_last=0.0
      X=0.0
      Y=0.0
      delta=0.0
      G=0.0
      puremode=0.0
      x_1=0.0
      x_2=0.0
      x_3=0.0
      x_4=0.0
      v=0.0
      v_mix=0.0
      czm_mix=0.0
      T_inter=0.0
      T_inter_mix=0.0   
!     STATE VARIABLES ASSIGNATION IN ARRAYS: d_old(k),G_last_old(k) 
      do i=1,(ncomp+1)
        if(i.le.ncomp)then
          G_last_old(i)=statev(i+4)
        else
          d_old=statev(i+4)
        endif
      enddo
!     DATA ALLOCATION IN ARRAYS: delta(k),v(i,j,k),G0(k),K0(k),v_mix(i,j,k)
      do k=1,ncomp
        if(k.eq.1)then
          delta(k)=0.5*(Strain(k)+abs(Strain(k)))
        else
          delta(k)=abs(Strain(k))
        endif
      enddo
      t_a=nint(prop(1)) 
      n=nint(prop(2))+1
      t1=nint(prop(3))+1
      if(ncomp.eq.3)then
        if(s_t.eq..true.)then
          t2=nint(prop(3))+1
        else  
          t2=nint(prop(4))+1
        endif
      else
        t2=0
      endif
      do k=1,ncomp
        z= 0.5*(k-1)*(k-2)*t2+(k-1)*(3-k)*t1+0.5*(2-k)*(3-k)*n   
        do i=1,2
          do j=2,n
            v(i,j,1)= prop(4+(i-1)*(n-1)+j-1)
          enddo
          do j=2,t1
            v(i,j,2)= prop(4+2*(n-1)+(i-1)*(t1-1)+j-1) 
          enddo
          if(ncomp.eq.3)then
            do j=2,t2 
              if(s_t.eq..true.)then !HP mode 2 ande mode 3 same czm law
                v(i,j,3)=v(i,j,2)
              else                  !HP mode 2 ande mode 3 different czm law           
                v(i,j,3)= prop(4+2*(n-1)+2*(t1-1)+(i-1)*(t2-1)+j-1)
              endif
            enddo
          else
          endif
        enddo 
        IF(t_a.eq.1)then 
          do j=2,z
          v(3,j,k)=0.5*(v(1,j,k)-v(1,j-1,k))*(v(2,j,k)+v(2,j-1,k))
          G0(k)=G0(k)+v(3,j,k)
          enddo
        ELSE
          do j=2,z 
            if(j.eq.2) then
              v(3,j,k)=0.5*v(2,j,k)*v(1,j,k)
            else
              do l=1,i_a-1    
              X(k)=v(1,j-1,k)+((v(1,j,k)-v(1,j-1,k))/i_a)*l
              T_inter=S(V,X)
              v(3,j,k)=v(3,j,k)+(2*T_inter(1,k))
              enddo
              v(3,j,k)=v(3,j,k)+v(2,j-1,k)+v(2,j,k)
              v(3,j,k)=v(3,j,k)*(v(1,j,k)-v(1,j-1,k))/(2*i_a)
            endif  
            G0(k)=G0(k)+v(3,j,k)
          enddo
        ENDIF
        K0(k)=v(2,2,k)/v(1,2,k) 
! v_mix(i,j,k)        
        do j=2,max(n,t1,t2)
          if(j.eq.2)then
            v_mix(1,j,k)=v(1,j,k)
            v_mix(2,j,k)=v(2,j,k)
          else
            IF(t_a.eq.1)then 
              v_mix(1,j,k)=v_mix(1,2,k)+((v(1,z,k)-v(1,2,k))/
     &                     (max(n,t1,t2)-2))*(j-2) !equally spaced nodes
            ELSE
              v_mix(1,j,k)=(0.5*(v(1,z,k)+v(1,1,k)))+(0.5*(v(1,z,k)-
     &                     v(1,1,k))*(-cos(((j*PI)/(max(n,t1,t2)))))) !nodes of Chebyshev Lobatto
            ENDIF
          endif
          Y(k)=v_mix(1,j,k)
          T_inter=S(V,Y)
          v_mix(2,j,k)=T_inter(1,k)
        enddo
        do j=2,max(n,t1,t2) 
          if(j.eq.2) then
            v_mix(3,j,k)=0.5*v_mix(2,j,k)*v_mix(1,j,k)
          else
            do l=1,i_a-1    
            X(K)=v_mix(1,j-1,k)+((v_mix(1,j,k)-v_mix(1,j-1,k))/i_a)*l
            T_inter=S(V,X)
            v_mix(3,j,k)=v_mix(3,j,k)+(2*T_inter(1,K))
            enddo
            v_mix(3,j,k)=v_mix(3,j,k)+v_mix(2,j-1,k)+v_mix(2,j,k)
            v_mix(3,j,k)=v_mix(3,j,k)*(v_mix(1,j,k)-
     &                   v_mix(1,j-1,k))/(2*i_a)
          endif  
        enddo
      enddo
!     CALCULATION OF ARRAYS VALUES: G_last(k),G_last_old(k),G(K)
      T_inter=S(V,delta)
      do k=1,ncomp
        z= 0.5*(k-1)*(k-2)*t2+(k-1)*(3-k)*t1+0.5*(2-k)*(3-k)*n    
        if(delta(k).le.v(1,2,k))then
          G_last(k)=G0(k)
        elseif(delta(k).ge.v(1,z,k))then
          G_last(k)=0
        else
          IF(t_a.eq.1)then 
            do j=2,z-1
              if(delta(k).ge.v(1,j,k).and.delta(k).lt.v(1,j+1,k))then 
                G_last(k)=0.5*(T_inter(1,k)*delta(k)+(v(1,j+1,k)-
     &                    delta(k))*(T_inter(1,k)+v(2,j+1,k)))  
              endif
              if(j.le.z-2) G_last(k)=G_last(k)+v(3,j+2,k)  
            enddo
          ELSE 
            do j=2,z-1
              if(delta(k).ge.v(1,j,k).and.delta(k).le.v(1,j+1,k))then 
                G_last(k)=0.5*T_inter(1,k)*delta(k)                 
                do l=1,i_a-1    
                  X(k)=delta(k)+((v(1,j+1,k)-delta(k))/i_a)*l
                  T_inter=S(V,X)
                  v(3,1,k)=v(3,1,k)+2*T_inter(1,K)
                enddo
                T_inter=S(V,delta)
                v(3,1,k)=v(3,1,k)+T_inter(1,k)+v(2,j+1,k)
                v(3,1,k)=v(3,1,k)*(v(1,j+1,k)-delta(k))/(2*i_a) 
                G_last(k)=G_last(k)+v(3,1,k)
              endif
              if(j.le.z-2) G_last(k)=G_last(k)+v(3,j+2,k)
            enddo
          ENDIF  
        endif
        G_last(k)=min(G_last_old(k),G_last(k)) 
        T_inter=S(V,delta)
        if(delta(k).le.v(1,2,k))then
          G(k)=G0(k)-G_last(k)+(0.5*(T_inter(1,k)*delta(k)))
        else 
          G(k)=G0(k)-G_last(k)
        endif
      enddo
!     CALCULATION OF: GT,puremode(k)
      do k=1,3
      GT=GT+G(k)
        if(G(k).eq.0)then 
          puremode(k)=0
        else
          puremode(k)=1
        endif
      enddo
!     CALCULATION VALUES OF: delta_mix,czm_mix(i,j),k0_mix,d_mix
      delta_mix=((delta(1)**2)+(delta(2)**2)+(delta(3)**2))**0.5
      do j=2,max(n,t1,t2)     
        czm_mix(3,j)=v_mix(3,j,1)+(((v_mix(3,j,2)-v_mix(3,j,1))*
     &    (G(2)/GT))+((v_mix(3,j,3)-v_mix(3,j,1))*(G(3)/GT)))*
     &    (((G(2)+G(3))/GT)**(ETA-1))
        czm_mix(2,j)=((v_mix(2,j,1)**2)+((((v_mix(2,j,2)**2)-
     &    (v_mix(2,j,1)**2))*(G(2)/GT))+(((v_mix(2,j,3)**2)-
     & (v_mix(2,j,1)**2))*(G(3)/GT)))*(((G(2)+G(3))/GT)**(ETA-1)))**0.5
        IF(t_a.eq.1)then 
          czm_mix(1,j)=czm_mix(1,j-1)+(2*czm_mix(3,j)/
     &                 (czm_mix(2,j)+czm_mix(2,j-1)))  
        ELSE
          if(j.eq.2)then
            czm_mix(1,j)=czm_mix(1,j-1)+(2*czm_mix(3,j)/
     &                   (czm_mix(2,j)+czm_mix(2,j-1))) 
          else  
            do k=1,3
              x_1(k)=v_mix(1,j-1,k)
              x_2(k)=v_mix(1,j-1,k)+(((v_mix(1,j,k)-v_mix(1,j-1,k))/3))
              x_3(k)=v_mix(1,j-1,k)+(((v_mix(1,j,k)-v_mix(1,j-1,k))/3)*
     &               2)
              x_4(k)=v_mix(1,j,k)
            enddo 
            T_inter=S(V,x_1)
            y_1=(((T_inter(1,1)**2)+((((T_inter(1,2)**2)-
     &(T_inter(1,1)**2))*(G(2)/GT))+(((T_inter(1,3)**2)-(T_inter(1,1)**
     & 2))*(G(3)/GT)))*(((G(2)+G(3))/GT)**(ETA-1)))**0.5)
            T_inter=S(V,x_2) 
            y_2=(((T_inter(1,1)**2)+((((T_inter(1,2)**2)-
     &(T_inter(1,1)**2))*(G(2)/GT))+(((T_inter(1,3)**2)-(T_inter(1,1)**
     & 2))*(G(3)/GT)))*(((G(2)+G(3))/GT)**(ETA-1)))**0.5)
            T_inter=S(V,x_3) 
            y_3=(((T_inter(1,1)**2)+((((T_inter(1,2)**2)-
     &(T_inter(1,1)**2))*(G(2)/GT))+(((T_inter(1,3)**2)-(T_inter(1,1)**
     & 2))*(G(3)/GT)))*(((G(2)+G(3))/GT)**(ETA-1)))**0.5)
            T_inter=S(V,x_4) 
            y_4=(((T_inter(1,1)**2)+((((T_inter(1,2)**2)-
     &(T_inter(1,1)**2))*(G(2)/GT))+(((T_inter(1,3)**2)-(T_inter(1,1)**
     & 2))*(G(3)/GT)))*(((G(2)+G(3))/GT)**(ETA-1)))**0.5)    
            czm_mix(1,j)=czm_mix(1,j-1)+(czm_mix(3,j)*(8/(y_1+(3*y_2)+
     &                   (3*y_3)+y_4)))
          endif 
        ENDIF    
      enddo
      k0_mix=czm_mix(2,2)/czm_mix(1,2)
      if(delta_mix.le.czm_mix(1,2).or.delta_mix.eq.0)then
        d_mix=0
      elseif(delta_mix.ge.czm_mix(1,max(n,t1,t2)))then 
        d_mix=1
      else
        if((puremode(1)+puremode(2)+puremode(3)).eq.1)then
          do k=1,3
            if(delta(K).eq.0)then
              d(k)=0
              d_der(k)=0
            else
              T_inter=S(V,delta)
              d(k)=1-(T_inter(1,k))/(delta(k)*k0(k))
              d_der(k)=(T_inter(1,k)-T_inter(2,k)*delta(K))/
     &                 (k0(k)*delta(K)**3)
            endif
          enddo
          d_mix=d(1)+d(2)+d(3)
          d_mix_der=d_der(1)+d_der(2)+d_der(3)
        else
        T_inter_mix=S_mix(czm_mix,delta_mix)
        d_mix=1-(T_inter_mix(1))/(delta_mix*k0_mix)
        d_mix_der=(T_inter_mix(1)-T_inter_mix(2)*delta_mix)/
     &            (k0_mix*delta_mix**3)
        endif
      endif
      if(d_mix.gt.1) d_mix=1
!     DETERMINATION OF LOADING CONDITION IN MIXED MODE (DAMAGE EVOLUTION) 
      w_0=((czm_mix(1,max(n,t1,t2)))/100)
      if(delta_mix.le.czm_mix(1,2))then
        w_2=(delta_mix+w_0)*k0_mix
        w_1=(delta_mix-w_0)*k0_mix
      else  
        T_inter_mix=S_mix(czm_mix,(delta_mix+w_0))
        w_2=T_inter_mix(1) 
        T_inter_mix=S_mix(czm_mix,(delta_mix-w_0))
        w_1=T_inter_mix(1) 
      endif  
      R=(w_1*(delta_mix+w_0)-w_2*(delta_mix-w_0))/(w_1-w_2+(k0_mix*
     &  (1-d_old)*2*w_0))
      if(delta_mix.le.R) d_mix=d_old !unloding condition
!     CALCULATION OF ARRAY VALUES: dsdePl(k,k)
      if(delta_mix.le.R)then        !unloading condition
        do i=1,ncomp
          do j=1,ncomp
            if(i.eq.j)then
              if(i.eq.1)then
                if(Strain(i).le.0)then
                  dsdePl(i,j)=k0(i)
                else
                  dsdePl(i,j)=k0(i)*(1-d_old)
                endif
              else  
                dsdePl(i,j)=k0(i)*(1-d_old)
              endif
            else 
              dsdePl(i,j)=0
            endif
          enddo
        enddo
      else                          !loading condition
        do i=1,ncomp
          do j=1,ncomp
            if(i.eq.j)then
              if(i.eq.1)then
                if(Strain(1).le.0)then
                  dsdePl(i,j)=k0(i)
                else
                  dsdePl(i,j)=k0(i)*(1-d_mix*(delta(i)/Strain(i))-
     &                        d_mix_der*Strain(i)*Strain(j))
                endif
              else  
                dsdePl(i,j)=k0(i)*(1-d_mix-d_mix_der*Strain(i)*
     &                      Strain(j))
              endif
            else
              if(Strain(1).lt.0)then
                dsdePl(i,j)=0
              else
                dsdePl(i,j)=k0(i)*(-d_mix_der*Strain(i)*Strain(j))
              endif
            endif
          enddo
        enddo
      endif    
! DETERMINATION OF STATE VARIABLES EACH SUBSTEP: d_old(k),G_last_old(k)
      d_old=d_mix
      do k=1,ncomp    
        G_last_old(k)=G_last(k)
      enddo 
      do i=1,(ncomp+1)
        if(i.le.ncomp)then
          statev(i+4)=G_last_old(i)
        else
          statev(i+4)=d_old
        endif
      enddo
! CALCULATION OF STRESS VALUES
      do k=1,ncomp
        if(k.eq.1)then
          Stress(k)=(1-d_mix)*k0(k)*Strain(k)-d_mix*k0(k)*(0.5*
     &              (-Strain(k)+abs(Strain(k))))
        else
          Stress(k)=(1-d_mix)*k0(k)*Strain(k)
        endif
      enddo   
c
  900 continue
      statev(2) = kchstat
      statev(3) = kestat 
      statev(4) = stiffi 
      sedEl = 0.d0
      sedPl = 0.d0     
c
  990 continue
  
      if (debugflag .gt. 0) then
         kwr = 1
      end if
c
      if(kwr .gt. 0) then
         call pplock(LOCKOT)
         iott = wrinqr(2)
         write(iott,1000) 'userCZM debug :',
     &   ' elemId=',elemId,
     &   ' kMatIntPt=',kMatIntPt,
     &   ' ldstep=',ldstep,
     &   ' isubst=',isubst,
     &   ' time=',Time,
     &   ' dtime=',dTime,
     &   ' Temp=',Temp,
     &   ' dTemp=',dTemp,
     &   ' coord=',coords
         write(iott,'(a,8e12.5/(11x,8e12.5))') 'prop      :',prop
         write(iott,'(a,6e12.5)') 'stress    :',stress
         write(iott,'(a,6e12.5)') 'Strain    :',Strain
         write(iott,'(a,6e12.5)') 'dsdePl(i,i):',
     &    (dsdePl(i,i),i=1,ncomp)
         write(iott,'(a,6e12.5/(10x,6e12.5))')
     &                 'statev    :',(statev(i),i=1,nstatev)
         call ppunlock(LOCKOT)
 1000    format(a/4(a,i5,1x),4(a,i1,1x)/5x,7(a,e12.5,1x),a,3e12.4)
      end if

      return
      end
!*****************************************************************************************************************************************************************************
! THIS FUNCTION S RETURNS AN ARRAY 3x2 THAT CONTAINS VALUES OF: INTERPOLATING FUNCTION AND DERIVATIVE, OF EACH COHESIVE LAWS IN THE 3 DIRECTIONS. 
      pure function S(V,delta)       
      implicit none
      integer i,j,k,l,z,t_a,n,t1,t2,nco
      COMMON /INDICES/ t_a,n,t1,t2,nco
      double precision, dimension(:,:,:), intent(in) :: V
      double precision, dimension(3), intent(in) :: delta
      double precision, dimension(2,3) :: S
      double precision, dimension(max(n,t1,t2),max(n,t1,t2),nco) :: A ! Array that contains Divided Differences, used to calculate coefficient of polynomial.
      double precision, dimension((max(n,t1,t2)-1),
     &                           (max(n,t1,t2)-1),nco) :: A1,A2       ! Arrays used to calculate array C.
      double precision, dimension(max(n,t1,t2),nco) :: B,C            ! Arrays used to calculate powers of the terms of the polynomial function and its derivative. 
      S=0.0
      A=0.0
      A1=0.0
      A2=0.0
      B=0.0
      C=0.0
!************************************************LINEAR INPERPOLATION********************************************************************************
      IF(t_a.eq.1)then 
!     CALCULATION OF TRACTIONS  
        do k=1,nco
          z= 0.5*(k-1)*(k-2)*t2+(k-1)*(3-k)*t1+0.5*(2-k)*(3-k)*n 
          do j=2,z-1
            if(delta(k).ge.v(1,j,k).and.delta(k).lt.v(1,j+1,k))then
              S(1,k)=(v(2,j,k)*(v(1,j+1,k)-delta(k))+v(2,j+1,k)*
     &               (delta(k)-v(1,j,k)))/(v(1,j+1,k)-v(1,j,k))
!     CALCULATION OF DERIVATIVES OF T-S LAWS (dsdePl consistent tangent jacobian matrix)
              S(2,k)=(v(2,j+1,k)-v(2,j,k))/(v(1,j+1,k)-v(1,j,k))
            endif 
          enddo  
        enddo
!********************************************POLYNOMIAL INTERPOLATION********************************************************************************
      ELSE
!     CALCULATION OF TRACTIONS 
        do k=1,nco  
          z= 0.5*(k-1)*(k-2)*t2+(k-1)*(3-k)*t1+0.5*(2-k)*(3-k)*n 
          do i=1,z
            do j=i,z
              if(i.eq.1)then
                A(i,j,k)=v(2,j,k)
                B(i,k)=1
              else
                A(i,j,k)=(A(i-1,j-1,k)-A(i-1,j,k))/
     &                   (v(1,j-i+1,k)-v(1,j,k)) 
                B(i,k)=B(i-1,k)*(delta(k)-v(1,i-1,k))  
              endif
            enddo
            S(1,k)=S(1,k)+A(i,i,k)*B(i,k)
          enddo 
!     CALCULATION OF DERIVATIVES OF T-S LAWS (dsdePl consistent tangent jacobian matrix) 
          do i=1,z-1
            do j=1,z-1
              if(i.eq.j)then
                A1(i,j,k)=1
              else
                A1(i,j,k)=delta(k)-v(1,j,k)
              endif
            enddo          
            do j=1,i
              A2(i,j,k)=1  
              do l=1,i
                A2(i,j,k)=A2(i,j,k)*A1(j,l,k)
              enddo
            enddo
          enddo
          do i=1,z
            do j=1,i-1
              if(i.eq.1)then
                C(i,k)=0
              else  
                C(i,k)=C(i,k)+A2(i-1,j,k) 
              endif
            enddo
            S(2,k)=S(2,k)+A(i,i,k)*C(i,k)
          enddo 
        enddo
      ENDIF
      endfunction S
!*****************************************************************************************************************************************************************************   
! This function S_mix returns an array 1x2 that contains values of: interpolating function and derivative, of cohesive law in mixed-mode.
      pure function S_mix(czm_mix,delta_mix)
      implicit none
      integer i,j,k,l,z,t_a,n,t1,t2,nco
      COMMON /INDICES/ t_a,n,t1,t2,nco
      double precision, dimension(:,:), intent(in) :: czm_mix
      double precision, intent(in) :: delta_mix
      double precision, dimension(2) :: S_mix
      double precision, dimension(max(n,t1,t2),max(n,t1,t2)) :: A  ! Array that contains Divided Differences, used to calculate coefficient of polynomial.
      double precision, dimension((max(n,t1,t2)-1),
     &                           (max(n,t1,t2)-1)) :: A1,A2        ! Arrays used to calculate array C.
      double precision, dimension(max(n,t1,t2)) :: B,C             ! Arrays used to calculate powers of the terms of the polynomial function and its derivative. 
      S_mix=0.0
      A=0.0
      A1=0.0
      A2=0.0
      B=0.0
      C=0.0
!************************************************LINEAR INTERPOLATION********************************************************************************
      IF(t_a.eq.1)then 
!     CALCULATION OF TRACTIONS OF MIXED-MODE      
        do j=2,max(n,t1,t2)-1
          if(delta_mix.ge.czm_mix(1,j).and.delta_mix.lt.
     &       czm_mix(1,j+1))then
            S_mix(1)=(czm_mix(2,j)*(czm_mix(1,j+1)-delta_mix)+
     &               czm_mix(2,j+1)*(delta_mix-czm_mix(1,j)))/
     &               (czm_mix(1,j+1)-czm_mix(1,j))
            S_mix(2)=(czm_mix(2,j+1)-czm_mix(2,j))/
     &               (czm_mix(1,j+1)-czm_mix(1,j))
          endif
        enddo  
!********************************************POLYNOMIAL INTERPOLATION********************************************************************************
      ELSE
        do i=1,max(n,t1,t2)
          do j=i,max(n,t1,t2)
            if(i.eq.1)then
              A(i,j)=czm_mix(2,j)
              B(i)=1
            else
              A(i,j)=(A(i-1,j-1)-A(i-1,j))/
     &               (czm_mix(1,j-i+1)-czm_mix(1,j)) 
              B(i)=B(i-1)*(delta_mix-czm_mix(1,i-1))  
            endif
          enddo
            S_mix(1)=S_mix(1)+A(i,i)*B(i)
        enddo 
!     CALCULATION OF DERIVATIVES OF MIXED-MODE T-S LAWS 
        do i=1,max(n,t1,t2)-1
          do j=1,max(n,t1,t2)-1
            if(i.eq.j)then
              A1(i,j)=1
            else
              A1(i,j)=delta_mix-czm_mix(1,j)
            endif
          enddo          
          do j=1,i
            A2(i,j)=1  
            do l=1,i
              A2(i,j)=A2(i,j)*A1(j,l)
            enddo
          enddo
        enddo
        do i=1,max(n,t1,t2)
          do j=1,i-1
            if(i.eq.1)then
              C(i)=0
            else  
              C(i)=C(i)+A2(i-1,j) 
            endif
          enddo
          S_mix(2)=S_mix(2)+A(i,i)*C(i)
        enddo   
      ENDIF
      endfunction S_mix

      
      