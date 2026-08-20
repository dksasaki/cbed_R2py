#-----------------------------------------------------------------------------
# Contains the updated CBED steady-state model function.
#
# cbed_model() extracted verbatim from ss_global_parallel.R (the most
# up-to-date version of the model), so it can be sourced as a standalone
# function library the same way cbed_v1_func.R is, instead of sourcing
# ss_global_parallel.R itself (which is a full driver script: it sets
# hardcoded local paths, runs a parallel foreach over a global grid, writes
# netCDF output, and ends with stop("run complete")).
#
# Differences vs. cbed_v1_func.R's cbed_model():
#   - f.OM1/f.OM2/f.OM3 are no longer computed internally from
#     switch.oxy.dependent.k; the caller must supply them in `pars`
#     (in ss_global_parallel.R these are depth-dependent, see lines 485-487).
#   - k2, k3, k.adj.anoxia rate constants changed.
#   - R.TAlk uses "- 1*R.ODUox" instead of "- 2*R.ODUox".
#   - dODU scales OduDepo by svf.grid/por.grid.
#   - No manual_control_bioturbation / manual_control_bioirrigation /
#     manual_OM_decay_rate branches.
#   - derive() returns flux diagnostics (o2_flux, nh4_flux, no3_flux,
#     tot_denit, anammox, dic_flux, OM_flux, OM_burial) directly instead of
#     the richer rates/transport/grid objects cbed_v1_func.R returns.
#=======================================================================

library(ReacTran)
library(deSolve)
library(marelac)
library(rootSolve)

cbed_model <- function(pars, ...) {

  derive <- function(times, state, pars) {
    with(as.list(pars), {

      #state varibles
      OM1      <- state[(0*pars$N+1):(1*N)]
      OM2    <- state[(1*N+1):(2*N)]
      OM3    <- state[(2*N+1):(3*N)]
      O2     <- state[(3*N+1):(4*N)]
      NH4      <- state[(4*N+1):(5*N)]
      NO3      <- state[(5*N+1):(6*N)]
      ODU      <- state[(6*N+1):(7*N)]  # ODU = oxygen demand unit
      DIC      <- state[(7*N+1):(8*N)]
      TAlk      <- state[(8*N+1):(9*N)]

      #reactions
      # temperature dependence | multiply this with reactions to get temperature dependence. To turn off, set it to 1
      Q10.temp <- Q10^((temp-4)/10)

      # # fix the OM1,2,3 pool ratio
      # if (switch.oxy.dependent.k==1) {
      #   f.OM1 <- 0.70
      #   f.OM2 <- 0.20
      #   f.OM3 <- 0.10
      #   pars$f.OM1 <- f.OM1
      #   pars$f.OM2 <- f.OM2
      #   pars$f.OM3 <- f.OM3
      # }
      # else {
      #
      #   f.OM1 <- 0.70
      #   f.OM2 <- 0.27
      #   f.OM3 <- 0.03
      #   pars$f.OM1 <- f.OM1
      #   pars$f.OM2 <- f.OM2
      #   pars$f.OM3 <- f.OM3
      # }

      #define three pools of OM
      J.OM1 <- f.OM1*J.OM
      J.OM2 <- f.OM2*J.OM
      J.OM3 <- f.OM3*J.OM


      # Change in k with depth Archer et al 2002
      k_scale_z <- 1 # exp(-grid$x.mid/3.5)

      #oxygen dependent OM remin rate switch
      if (switch.oxy.dependent.k==1) {
        #reactivity as function of total POC flux

        k.adj.denit = 0.1        #0.1     # 1 = no effect of denitrification | 0.1 = decay rate decreases by a factor of 10. Generally 1/10
        k.adj.anoxia = 0.001     #0.005
        pars$k.adj.denit <- k.adj.denit
        pars$k.adj.anoxia <- k.adj.anoxia

        k1 <- (1.5*10^-1)*(J.OM)^0.85   # (1.5*10^-1)
        k2 <- (1.5*10^-3)*(J.OM)^0.85   # (2.3*10^-3)
        k3 <- (0.9*10^-4)*(J.OM)^0.85   # (1.3*10^-4)

        #used in cbed paper 1.
        # k1 <- (1.5*10^-1)*(J.OM)^0.85   # (1.5*10^-1)
        # k2 <- (2.3*10^-3)*(J.OM)^0.85   # (2.3*10^-3)
        # k3 <- (1.3*10^-4)*(J.OM)^0.85   # (1.3*10^-4)

        pars$k1 <- k1
        pars$k2 <- k2
        pars$k3 <- k3
      }
      else {
        k.adj.denit = 1      # 1 = no effect of denitrification | 0.1 = decay rate decreases by a factor of 10. Generally 1/10
        k.adj.anoxia = 1
        pars$k.adj.denit <- k.adj.denit
        pars$k.adj.anoxia <- k.adj.anoxia
        #reactivity as function of total POC flux
        k1 <- (1.5*10^-1)*(J.OM)^0.85 *k_scale_z  # (1.5*10^-1)
        k2 <- (1.3*10^-4)*(J.OM)^0.85 *k_scale_z  # (1.3*10^-4)
        k3 <- 0*(1.3*10^-5)*(J.OM)^0.85 *k_scale_z  # (1.3*10^-5)
        pars$k1 <- k1
        pars$k2 <- k2
        pars$k3 <- k3
      }


      ## bgc reactions
      #aerobic respiration
      R.O2.1 <- k1*OM1*(O2/(Ks.O2+O2)) *Q10.temp
      R.O2.2 <- k2*OM2*(O2/(Ks.O2+O2)) *Q10.temp
      R.O2.3 <- k3*OM3*(O2/(Ks.O2+O2)) *Q10.temp

      #OM oxidation by nitrate (denitrification)
      R.NO3.1 <- k.adj.denit * k1*OM1*(NO3/(Ks.NO3+NO3))*(Ks.O2/(Ks.O2+O2)) *Q10.temp
      R.NO3.2 <- k.adj.denit * k2*OM2*(NO3/(Ks.NO3+NO3))*(Ks.O2/(Ks.O2+O2)) *Q10.temp
      R.NO3.3 <- k.adj.denit * k3*OM3*(NO3/(Ks.NO3+NO3))*(Ks.O2/(Ks.O2+O2)) *Q10.temp

      #OM oxidation by all other electron acceptors
      R.ODU.1 <- k.adj.anoxia * k1*OM1*(Ks.O2/(Ks.O2+O2))*(Ks.NO3/(Ks.NO3+NO3)) *Q10.temp
      R.ODU.2 <- k.adj.anoxia * k2*OM2*(Ks.O2/(Ks.O2+O2))*(Ks.NO3/(Ks.NO3+NO3)) *Q10.temp
      R.ODU.3 <- k.adj.anoxia * k3*OM3*(Ks.O2/(Ks.O2+O2))*(Ks.NO3/(Ks.NO3+NO3)) *Q10.temp

      #total OM remineralization
      R.DIC.1 <- R.O2.1+R.NO3.1+R.ODU.1
      R.DIC.2 <- R.O2.2+R.NO3.2+R.ODU.2
      R.DIC.3 <- R.O2.3+R.NO3.3+R.ODU.3

      #nitrification
      R.NOx <- k.NOx*NH4*O2 *Q10.temp

      #anammox | NH4 + NO3 = N2
      R.ana <- k.ana*NH4*NO3 *Q10.temp

      #ODU reoxidation
      R.ODUox <- k.ODUox*ODU*O2 *Q10.temp

      #ODU burial as solid (pyrite burial) #Soetaert et al 1996
      OduDepo    <- (R.ODU.1+R.ODU.2+R.ODU.3)*min(1, 0.233*(w)^0.336 ) *switch.odu.burial

      #calculate TAlk generation and consumption
      R.TAlk <- rNC*svf.grid$mid/por.grid$mid*(R.O2.1+R.O2.2+R.O2.3) +
        (0.8+rNC)*svf.grid$mid/por.grid$mid*(R.NO3.1+R.NO3.2+R.NO3.3) +
        (1+rNC)*svf.grid$mid/por.grid$mid*(R.ODU.1+R.ODU.2+R.ODU.3) -
        2*R.NOx - 1*R.ODUox



      #transport
      trans.OM1 <- tran.1D(C=OM1, flux.up = J.OM1, D=Db.grid, v=v.grid, dx=grid,VF=svf.grid)
      trans.OM2 <- tran.1D(C=OM2, flux.up = J.OM2, D=Db.grid, v=v.grid, dx=grid,VF=svf.grid)
      trans.OM3 <- tran.1D(C=OM3, flux.up = J.OM3, D=Db.grid, v=v.grid, dx=grid,VF=svf.grid)

      trans.O2 <- tran.1D(C=O2, C.up=O2.w, D=Do2.grid, v=u.grid, dx = grid, VF=por.grid)
      trans.NH4 <- tran.1D(C=NH4, C.up=NH4.w, D=D.NH4, v=w.NH4, dx = grid)
      trans.NO3 <- tran.1D(C=NO3, C.up=NO3.w, D=Dno3.grid, v=u.grid, dx = grid, VF=por.grid)
      trans.ODU <- tran.1D(C=ODU, C.up=ODU.w, D=Dodu.grid, v=u.grid, dx = grid, VF=por.grid)
      trans.DIC <- tran.1D(C=DIC, C.up=DIC.w, D=Ddic.grid, v=u.grid, dx = grid, VF=por.grid)
      trans.TAlk <- tran.1D(C=TAlk, C.up=TAlk.w, D=Ddic.grid, v=u.grid, dx = grid, VF=por.grid)

      #ODEs
      dOM1 <- trans.OM1$dC - R.O2.1 - R.NO3.1 - R.ODU.1
      dOM2 <- trans.OM2$dC - R.O2.2 - R.NO3.2 - R.ODU.2
      dOM3 <- trans.OM3$dC - R.O2.3 - R.NO3.3 - R.ODU.3

      dO2 <- trans.O2$dC - svf.grid$mid/por.grid$mid*(R.O2.1+R.O2.2+R.O2.3) - 2*R.NOx - R.ODUox + alpha$mid*(O2.w-O2)


      dNH4 <- 1/phi.N$mid*trans.NH4$dC + svf.grid$mid/phi.N$mid*(rNC.1*R.O2.1+rNC.2*R.O2.2 + rNC.3*R.O2.3 +
                                                                   rNC.1*R.NO3.1+rNC.2*R.NO3.2+rNC.3*R.NO3.3+
                                                                   rNC.1*R.ODU.1+rNC.2*R.ODU.2+rNC.3*R.ODU.3) -
        por.grid$mid/phi.N$mid*(R.NOx+R.ana) + alpha$mid*(NH4.w-NH4)


      dNO3 <- trans.NO3$dC - svf.grid$mid/por.grid$mid*0.8*(R.NO3.1+R.NO3.2+R.NO3.3) + R.NOx - R.ana + alpha$mid*(NO3.w-NO3)

      dODU <- trans.ODU$dC + svf.grid$mid/por.grid$mid*(R.ODU.1+R.ODU.2+R.ODU.3) - R.ODUox - svf.grid$mid/por.grid$mid*OduDepo + alpha$mid*(ODU.w-ODU)

      dDIC <- trans.DIC$dC + svf.grid$mid/por.grid$mid*(R.DIC.1+R.DIC.2+R.DIC.3) + alpha$mid*(DIC.w-DIC)
      dTAlk <- trans.TAlk$dC + R.TAlk + alpha$mid*(TAlk.w-TAlk)

      return(list(c(dOM1=dOM1,
                    dOM2=dOM2,
                    dOM3=dOM3,
                    dO2=dO2,
                    dNH4=dNH4,
                    dNO3=dNO3,
                    dODU=dODU,
                    dDIC=dDIC,
                    dTAlk=dTAlk),
                  o2_flux = trans.O2$flux.up + sum(por.grid$mid*grid$dx*alpha$mid*(pars$O2.w-O2)),
                  nh4_flux = trans.NH4$flux.up + sum(por.grid$mid*grid$dx*alpha$mid*(pars$NH4.w-NH4)),
                  no3_flux = trans.NO3$flux.up + sum(por.grid$mid*grid$dx*alpha$mid*(pars$NO3.w-NO3)),
                  tot_denit = 0.8*sum(svf.grid$mid*grid$dx*(R.NO3.1+R.NO3.2+R.NO3.3)) + 2*sum(por.grid$mid*grid$dx*R.ana),
                  anammox = 2*sum(por.grid$mid*grid$dx*R.ana),
                  dic_flux = trans.DIC$flux.up + sum(por.grid$mid*grid$dx*alpha$mid*(pars$DIC.w-DIC)),
                  OM_flux = trans.OM1$flux.up+trans.OM2$flux.up+trans.OM3$flux.up,
                  OM_burial = trans.OM1$flux.down+trans.OM2$flux.down+trans.OM3$flux.down
      ))


    })

  }

  #grid <- setup.grid.1D(x.up = 0,x.down = 50, N = pars$N)
  #grid <- setup.grid.1D(x.up=0,x.down=pars$L, N=pars$N)
  grid <- setup.grid.1D(x.up=0,x.down=pars$L, N=pars$N,dx.1=0.1)

  #Set Db as Archar 2002
  Db.0 <- 0.0232*(pars$J.OM)^0.85
  Db.grid<- setup.prop.1D(xy = cbind(x=grid$x.mid,
                                     y=(Db.0*exp(-(grid$x.mid/8)^2)*((pars$O2.w*1000)/(pars$O2.w*1000+20)))
  ),grid=grid)


  # Set up porosity and solid volume fraction.
  por.grid <- setup.prop.1D(p.exp,grid=grid,y.0=pars$por.0,
                            y.inf=pars$por.inf,x.att=pars$por.att)
  svf.grid <- setup.prop.1D(xy=cbind(grid$x.mid,1-por.grid$mid),grid=grid)
  # Set up advection terms for solids (v) and solutes (u)
  u.grid <- setup.compaction.1D(v.0=pars$w,por.0=pars$por.0,por.inf=pars$por.inf,
                                por.grid=por.grid)$u
  v.grid <- setup.compaction.1D(v.0=pars$w,por.0=pars$por.0,por.inf=pars$por.inf,
                                por.grid=por.grid)$v

  # Set up diffusion terms
  D <- diffcoeff(S=pars$S,t=pars$temp,species=c("O2","NO3","NH4","Fe","SO4","H2S","CH4","Mn","HCO3"))*100^2*365*86400
  irr.enh <- setup.prop.1D(func = p.exp, grid=grid, y.0=pars$irr.enh.0, y.inf=1,x.L=pars$irr.enh.L)

  # Correct for tortuosity
  Do2.grid  <- setup.prop.1D(xy=cbind(grid$x.mid,
                                      Db.grid$mid+(D$O2*irr.enh$mid)/(1-2*log(por.grid$mid))),
                             grid=grid)
  Dno3.grid <- setup.prop.1D(xy=cbind(grid$x.mid,
                                      Db.grid$mid+(D$NO3*irr.enh$mid)/(1-2*log(por.grid$mid))),
                             grid=grid)
  Dnh4.grid <- setup.prop.1D(xy=cbind(grid$x.mid,
                                      Db.grid$mid+(D$NH4*irr.enh$mid)/(1-2*log(por.grid$mid))),
                             grid=grid)
  Dodu.grid <- setup.prop.1D(xy=cbind(grid$x.mid,
                                      Db.grid$mid+(D$H2S*irr.enh$mid)/(1-2*log(por.grid$mid))),
                             grid=grid)
  Ddic.grid  <- setup.prop.1D(xy=cbind(grid$x.mid,
                                       Db.grid$mid+(D$HCO3*irr.enh$mid)/(1-2*log(por.grid$mid))),
                              grid=grid)
  # Set up parameters for ammonium.
  phi.N <- setup.prop.1D(xy=cbind(grid$x.mid, por.grid$mid+svf.grid$mid*pars$K.equilNH4*pars$ps),grid=grid)
  D.NH4 <- setup.prop.1D(xy=cbind(grid$x.mid, por.grid$mid*Dnh4.grid$mid + svf.grid$mid*pars$K.equilNH4*pars$ps*Db.grid$mid),
                         grid=grid)
  w.NH4 <- setup.prop.1D(xy=cbind(grid$x.mid, por.grid$mid*u.grid$mid + svf.grid$mid*pars$K.equilNH4*pars$ps*v.grid$mid),
                         grid=grid)


  #Bioirrigation as Archer 2002. with 'if loop' to handle complete anoxia
  if (pars$O2.w > 0) {
    I0 <- 11*(((atan((5*pars$J.OM-400)/400))/pi)+0.5) - 0.9 + 20*((pars$O2.w*1000)/(pars$O2.w*1000+10))*exp(-pars$O2.w*1000/10)*(pars$J.OM/(pars$J.OM+30))
  } else {
    I0 <- 0
  }
  alpha <- setup.prop.1D(xy=cbind(grid$x.mid,
                                  I0*exp(-(grid$x.mid/1.8)^2)),
                         grid = grid)


  #initial condition
  yini.ss <- c(rep(rep(0,grid$N),9))  #here 9 is the number of state variables


  # Solve the ODEs (ordinary differential equations)
  out <- steady.1D(y=yini.ss,func=derive,parms=pars,nspec=9,
                   names=c("OM1","OM2","OM3","O2","NH4","NO3", "ODU","DIC","TAlk"),
                   method="runsteady")
  # out <- steady.1D(y=yini.ss,func=derive,parms=pars,nspec=9,
  #                  names=c("OM1","OM2","OM3","O2","NH4","NO3", "ODU","DIC","TAlk"),
  #                  method="stode")

  return(out)

}
