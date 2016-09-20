C---------------------------------------------------
      subroutine update_properties

      include 'SIZE'
      include 'NEKUSE'
      include 'TOTAL'

      integer lxyz,ipoint
      parameter(lxyz=lx1*ly1*lz1)
      real h
      real density,viscosity,specificheat,conductivity
      real volexpcoef,volumefraction,temperature

      do ipoint=1,lxyz*nelv
        h=t(ipoint,1,1,1,ifld_h-1)
        dens(ipoint,1,1,1)=density(h)
        visc(ipoint,1,1,1)=viscosity(h)
        thcap(ipoint,1,1,1)=specificheat(h)
        thcond(ipoint,1,1,1)=conductivity(h)
        beta(ipoint,1,1,1)=volexpcoef(h)
        volfrac(ipoint,1,1,1)=volumefraction(h)
        temper(ipoint,1,1,1)=temperature(h)
      enddo

      return
      END
C---------------------------------------------------
      subroutine point_properties(h,ix,iy,iz,e)
      implicit none

      include 'FLUIDPROP'

      real h
      real density,viscosity,specificheat,conductivity
      real volexpcoef,volumefraction,temperature
      integer ix,iy,iz,e

      dens(ix,iy,iz,e)=density(h)
      visc(ix,iy,iz,e)=viscosity(h)
      thcap(ix,iy,iz,e)=specificheat(h)
      thcond(ix,iy,iz,e)=conductivity(h)
      beta(ix,iy,iz,e)=volexpcoef(h)
      volfrac(ix,iy,iz,e)=volumefraction(h)
      temper(ix,iy,iz,e)=temperature(h)

      return
      END
C---------------------------------------------------
      real function temperature(h)

      include 'SIZE'
      include 'FLUIDPROP'

      real h,cp
      real specificheat

      if(h.le.h_f) then
        cp=specificheat(h)
        temperature=T_sat+(h-h_f)/cp
      else if(h.ge.h_g) then
        cp=specificheat(h)
        temperature=T_sat+(h-h_g)/cp
      else
        temperature=T_sat
      endif

      return
      END
C---------------------------------------------------
      real function volumefraction(h)

      include 'SIZE'
      include 'FLUIDPROP'

      real h

      if(h.le.h_f) then
        volumefraction=0.0
      else if(h.ge.h_g) then
        volumefraction=1.0d0
      else
        volumefraction=rho_f*(h-h_f)/(rho_g*(h_g-h)+rho_f*(h-h_f))
      endif

      return
      END
C---------------------------------------------------
      real function dens_model(h,h0,r0,b0,h1,r1,b1)
      implicit none

      real h,h0,h1,r0,r1,b0,b1
      real K2,K3,hstar

      K2=6.0d0*log(r1/r0)/(h0-h1)-2.0d0*b1-4.0d0*b0
      K3=6.0d0*log(r1/r0)/(h1-h0)+3.0d0*(b1+b0)
      hstar=(h-h0)/(h1-h0)
      dens_model=
     &    r0*exp((h0-h1)*(K3/3.0d0*hstar**3+K2/2.0d0*hstar**2+b0*hstar))

      return
      END
C---------------------------------------------------
      real function beta_model(h,h0,b0,r0,h1,b1,r1) !-(1/rho)(drho/dh)
      implicit none

      real h,h0,h1,r0,r1,b0,b1
      real K2,K3,hstar

      K2=6.0d0*log(r1/r0)/(h0-h1)-2.0d0*b1-4.0d0*b0
      K3=6.0d0*log(r1/r0)/(h1-h0)+3.0d0*(b1+b0)
      hstar=(h-h0)/(h1-h0)
      beta_model=K3*hstar**2+K2*hstar+b0

      return
      END
C---------------------------------------------------
      real function dbeta_model(h,h0,b0,r0,h1,b1,r1) !-(1/rho)(drho/dh)
      implicit none

      real h,h0,h1,r0,r1,b0,b1
      real K2,K3,hstar

      K2=6.0d0*log(r1/r0)/(h0-h1)-2.0d0*b1-4.0d0*b0
      K3=6.0d0*log(r1/r0)/(h1-h0)+3.0d0*(b1+b0)
      hstar=(h-h0)/(h1-h0)
      dbeta_model=(2.0d0*K3*hstar+K2)/(h1-h0)

      return
      END
C---------------------------------------------------
      real function cubicpoly(xhat,x1,f1,f1ph,x2,f2,f2ph)
      implicit none

      real xhat,x1,x2,f1,f2,f1ph,f2ph
      real x,f1p,f2p,A,B,C,D

      x=(xhat-x1)/(x2-x1)
      f1p=f1ph*(x2-x1)
      f2p=f2ph*(x2-x1)
      A=f2p-2.0d0*f2+f1p+2.0d0*f1
      B=3.0d0*f2-f2p-2.0d0*f1p-3.0d0*f1
      C=f1p
      D=f1
      cubicpoly=A*x**3+B*x**2+C*x+D

      return
      end
C---------------------------------------------------
      real function intcubicpoly(xhat,x1,f1,f1ph,x2,f2,f2ph)
      implicit none

      real xhat,x1,x2,f1,f2,f1ph,f2ph
      real x,f1p,f2p,A,B,C,D

      x=(xhat-x1)/(x2-x1)
      f1p=f1ph*(x2-x1)
      f2p=f2ph*(x2-x1)
      A=f2p-2.0d0*f2+f1p+2.0d0*f1
      B=3.0d0*f2-f2p-2.0d0*f1p-3.0d0*f1
      C=f1p
      D=f1
      intcubicpoly=2.5d-1*A*x**4+(B/3.0d0)*x**3+5.0d-1*C*x**2+D*x

      return
      end
C---------------------------------------------------
      real function density(h)

      include 'SIZE'
      include 'FLUIDPROP'

      real h
      real dens_model,beta_model,dbeta_model,intcubicpoly,volumefraction
      real eps,hfm,hfp,hgm,hgp,rhofg,hfg
      real alpha,denom,rhom,dadh,d2adh2,f1,f2,f1p,f2p,rhofm,icp

      hfg=h_g-h_f
      eps=1.0d-4*hfg
      hfm=h_f-eps
      hfp=h_f+eps
      hgm=h_g-eps
      hgp=h_g+eps
      rhofg=rho_g-rho_f

      if(h.le.hfm) then
        density=
     &       dens_model(h,h_l_ref,rho_l_ref,beta_l_ref,h_f,rho_f,beta_f)
      elseif(h.le.hfp) then
        f1=beta_model(hfm,h_l_ref,beta_l_ref,rho_l_ref,h_f,beta_f,rho_f)
        f1p=
     &    dbeta_model(hfm,h_l_ref,beta_l_ref,rho_l_ref,h_f,beta_f,rho_f)
        denom=1.0d0/(rho_g*(h_g-hfp)+rho_f*(hfp-h_f))
        alpha=rho_f*(hfp-h_f)*denom
        rhom=alpha*rho_g+(1.0d0-alpha)*rho_f
        dadh=rho_f*rho_g*hfg*denom**2
        d2adh2=2.0d0*rho_f*rho_g*rhofg*hfg*denom**3
        f2=-rhofg/rhom*dadh
        f2p=(rhofg/rhom*dadh)**2-rhofg/rhom*d2adh2
        rhofm=
     &     dens_model(hfm,h_l_ref,rho_l_ref,beta_l_ref,h_f,rho_f,beta_f)
        icp=intcubicpoly(h,hfm,f1,f1p,hfp,f2,f2p)
        density=rhofm*exp(-(hfp-hfm)*icp)
      elseif(h.le.hgm) then
        alpha=volumefraction(h)
        density=(1.0d0-alpha)*rho_f+alpha*rho_g
      elseif(h.le.hgp) then
        denom=1.0d0/(rho_g*(h_g-hgm)+rho_f*(hgm-h_f))
        alpha=rho_f*(hgm-h_f)*denom
        rhom=alpha*rho_g+(1.0d0-alpha)*rho_f
        dadh=rho_f*rho_g*hfg*denom**2
        d2adh2=2.0d0*rho_f*rho_g*rhofg*hfg*denom**3
        f1=-rhofg/rhom*dadh
        f1p=(rhofg/rhom*dadh)**2-rhofg/rhom*d2adh2
        f2=beta_model(hgp,h_g,beta_g,rho_g,h_v_ref,beta_v_ref,rho_v_ref)
        f2p=
     &    dbeta_model(hgp,h_g,beta_g,rho_g,h_v_ref,beta_v_ref,rho_v_ref)
        icp=intcubicpoly(h,hgm,f1,f1p,hgp,f2,f2p)
        density=rhom*exp(-(hgp-hgm)*icp)
      elseif(h.le.h_v_ref) then
        density=
     &       dens_model(h,h_g,rho_g,beta_g,h_v_ref,rho_v_ref,beta_v_ref)
      else
        density=rho_v_ref
      endif

      return
      END
C---------------------------------------------------
      real function volexpcoef(h) !-(1/rho)(drho/dh)

      include 'SIZE'
      include 'FLUIDPROP'

      real beta_model,dbeta_model,cubicpoly
      real h
      real eps,hfm,hfp,hgm,hgp,rhofg,hfg
      real alpha,denom,rhom,dadh,d2adh2,f1,f2,f1p,f2p,rhofm,icp

      hfg=h_g-h_f
      eps=1.0d-4*hfg
      hfm=h_f-eps
      hfp=h_f+eps
      hgm=h_g-eps
      hgp=h_g+eps
      rhofg=rho_g-rho_f

      if(h.le.(hfm)) then !subcooled liquid
        volexpcoef=
     &       beta_model(h,h_l_ref,beta_l_ref,rho_l_ref,h_f,beta_f,rho_f)
      elseif(h.le.hfp) then !blending region
        f1=beta_model(hfm,h_l_ref,beta_l_ref,rho_l_ref,h_f,beta_f,rho_f)
        f1p=
     &    dbeta_model(hfm,h_l_ref,beta_l_ref,rho_l_ref,h_f,beta_f,rho_f)
        denom=1.0d0/(rho_g*(h_g-hfp)+rho_f*(hfp-h_f))
        alpha=rho_f*(hfp-h_f)*denom
        rhom=alpha*rho_g+(1.0d0-alpha)*rho_f
        dadh=rho_f*rho_g*hfg*denom**2
        d2adh2=2.0d0*rho_f*rho_g*rhofg*hfg*denom**3
        f2=-rhofg/rhom*dadh
        f2p=(rhofg/rhom*dadh)**2-rhofg/rhom*d2adh2
        volexpcoef=cubicpoly(h,hfm,f1,f1p,hfp,f2,f2p)
      elseif(h.le.hgm) then !two-phase mixture
        denom=1.0d0/(rho_g*(h_g-h)+rho_f*(h-h_f))
        alpha=rho_f*(h-h_f)*denom
        rhom=alpha*rho_g+(1.0d0-alpha)*rho_f
        dadh=rho_f*rho_g*hfg*denom**2
        volexpcoef=-rhofg/rhom*dadh
      elseif(h.le.hgp) then !blending region
        denom=1.0d0/(rho_g*(h_g-hgm)+rho_f*(hgm-h_f))
        alpha=rho_f*(hgm-h_f)*denom
        rhom=alpha*rho_g+(1.0d0-alpha)*rho_f
        dadh=rho_f*rho_g*hfg*denom**2
        d2adh2=2.0d0*rho_f*rho_g*rhofg*hfg*denom**3
        f1=-rhofg/rhom*dadh
        f1p=(rhofg/rhom*dadh)**2-rhofg/rhom*d2adh2
        f2=beta_model(hgp,h_g,beta_g,rho_g,h_v_ref,beta_v_ref,rho_v_ref)
        f2p=
     &    dbeta_model(hgp,h_g,beta_g,rho_g,h_v_ref,beta_v_ref,rho_v_ref)
        volexpcoef=cubicpoly(h,hgm,f1,f1p,hgp,f2,f2p)
      elseif(h.le.h_v_ref) then !superheated vapor
        volexpcoef=
     &       beta_model(h,h_g,beta_g,rho_g,h_v_ref,beta_v_ref,rho_v_ref)
      else  !constant density above h_v_ref
        volexpcoef=0.0
      endif

      return
      END
C---------------------------------------------------
      real function viscosity(h)

      include 'SIZE'
      include 'FLUIDPROP'

      real h,alpha,volumefraction

      if(h.le.h_f) then
        viscosity=mu_f
      elseif(h.ge.h_g) then
        viscosity=mu_g
      else
        alpha=volumefraction(h)
        viscosity=(1.0d0-alpha)*mu_f+alpha*mu_g
      endif

      return
      END
C---------------------------------------------------
      real function conductivity(h)

      include 'SIZE'
      include 'FLUIDPROP'

      real h,alpha,volumefraction

      if(h.le.h_f) then
        conductivity=k_f
      elseif(h.ge.h_g) then
        conductivity=k_g
      else
        alpha=volumefraction(h)
        conductivity=(1.0d0-alpha)*k_f+alpha*k_g
      endif

      return
      END
C---------------------------------------------------
      real function specificheat(h)

      include 'SIZE'
      include 'FLUIDPROP'

      real volumefraction,density
      real h,rho,alpha

      if(h.le.h_f) then
        specificheat=cp_f
      elseif(h.ge.h_g) then
        specificheat=cp_g
      else
        alpha=volumefraction(h)
        rho=density(h)
        specificheat=((1.0d0-alpha)*rho_f*cp_f+alpha*rho_g*cp_g)/rho
      endif

      return
      END
C---------------------------------------------------
