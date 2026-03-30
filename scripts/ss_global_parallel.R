library(ncdf4)
library(oce)
library(DescTools)
library(parallel)
library(doParallel)
library(ReacTran)
library(marelac)
library(deSolve)
library(tictoc)

tic()

setwd("~/Documents/")

#porosity
load("~/Documents/porosity_global.RData")

#btm boundary conditions
nc.path      <- "~/Documents/Fei_global_SR/ocean_cobalt_btm.2013-2017.ann.nc"
nc.path_btm1 <- "~/Documents/Fei_global_SR/ocean_cobalt_tracers_month_z.2013-2017.ann_btm.nc"
nc.path_btm2 <- "~/Documents/Fei_global_SR/ocean_cobalt_omip_tracers_month_z.2013-2017.ann_btm.nc"
nc.path_btm3 <- "~/Documents/Fei_global_SR/ocean_annual_z.2013-2017.ann_btm.nc"

nc      <- nc_open(nc.path)
nc.btm1 <- nc_open(nc.path_btm1)
nc.btm2 <- nc_open(nc.path_btm2)
nc.btm3 <- nc_open(nc.path_btm3)

xh <- nc$dim$xh$vals  #lon
yh <- nc$dim$yh$vals  #lat

# ocean grid 
ocean_grid <- nc_open('~/Documents/Fei_global_SR/ocean_cobalt_btm.static.nc')
#ocean_area <- ncvar_get(ocean_grid, 'areacello')

ocean_depth <- ncvar_get(ocean_grid, 'deptho')
btm_o2 <- ncvar_get(nc, "btm_o2")
btm_temp <- ncvar_get(nc, "btm_temp")

btm_no3 <- ncvar_get(nc.btm1, "no3")
btm_nh4 <- ncvar_get(nc.btm1, "nh4")
btm_dic <- ncvar_get(nc.btm2, "dissic")
btm_talk <- ncvar_get(nc.btm2, "talk")
btm_salinity <- ncvar_get(nc.btm3, "so")


fntot <- (ncvar_get(nc, "fndet_btm") +
            ncvar_get(nc, "fndi_btm") +
            ncvar_get(nc, "fnsm_btm") +
            ncvar_get(nc, "fnmd_btm") +
            ncvar_get(nc, "fnlg_btm"))

# to calculate sedimentation rate. Porosity adjustment will be made in pars loop
w_tem <- ( (ncvar_get(nc, "fcadet_arag_btm")*100/2.71 + 
          ncvar_get(nc, "fcadet_calc_btm")*100/2.94) +
         (ncvar_get(nc, "fsidet_btm" ) + 
            ncvar_get(nc,"fsimd_btm") + 
            ncvar_get(nc, "fsilg_btm"))*60/2.65 +
         (ncvar_get(nc, "ffedet_btm") +
            ncvar_get(nc,"ffedi_btm")+
            ncvar_get(nc, "ffesm_btm")+
            ncvar_get(nc, "ffemd_btm") +
            ncvar_get(nc,"ffelg_btm"))*160/5.24 +
         (ncvar_get(nc, "fpdet_btm") +
            ncvar_get(nc, "fpdi_btm") + 
            ncvar_get(nc, "fpsm_btm") + 
            ncvar_get(nc, "fpmd_btm") +
            ncvar_get(nc,"fplg_btm" ))*120/2.3 + 
         ncvar_get(nc, "flithdet_btm")/2.65 + 
         fntot*6.625*22.4/0.9 )/1e4*60*60*24*365    #/(1-0.8)   # cm/yr


# the steady state model. (copy pasted)
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
      
      # fix the OM1,2,3 pool ratio
      if (switch.oxy.dependent.k==1) {
        f.OM1 <- 0.70
        f.OM2 <- 0.20 
        f.OM3 <- 0.10
        pars$f.OM1 <- f.OM1
        pars$f.OM2 <- f.OM2
        pars$f.OM3 <- f.OM3
      }
      else {
        
        f.OM1 <- 0.70
        f.OM2 <- 0.27 
        f.OM3 <- 0.03
        pars$f.OM1 <- f.OM1
        pars$f.OM2 <- f.OM2
        pars$f.OM3 <- f.OM3
      }
      
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
        k.adj.anoxia = 0.005     #0.005
        pars$k.adj.denit <- k.adj.denit
        pars$k.adj.anoxia <- k.adj.anoxia
        
        k1 <- (1.5*10^-1)*(J.OM)^0.85   # (1.5*10^-1)
        k2 <- (2.3*10^-3)*(J.OM)^0.85   # (2.3*10^-3)
        k3 <- (1.3*10^-4)*(J.OM)^0.85   # (1.3*10^-4)
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

# load pars and run model

#parameter list
pars <- list(
  #J.OM=11.05*365/10, #11.05 8.7
  
  k.NOx=2e5,    #2e5   #1e4  # units: [mmoles^-1 L a^-1]
  k.ana=1e5,       #1e5 
  k.ODUox=1e5,     #1e6 5e6 for Fe2+  
  
  Ks.O2 = 0.008,      #0.008   # units: [umol cm^-3] #RADI: 0.010
  Ks.NO3 = 0.001,     #0.001   # units: [umol cm^-3]  #RADI: 0.030
  
  K.equilNH4=1.6,       # units: [cm^3_pw g]
  Db.max=48, #9 12         # units: [cm^2 a^-1]
  xb.Db=2,
  xb=1,                 # units: [cm]
  xL=12, #12
  alpha.max=10, #10 133       #0.3*365   #bio-irrigation scaling coefficient [a^-1]
  xL.alpha = 20,
  alpha.zi=3,    # bioirrigation depth, where alpha becomes 0. Bohlen etal 2011
  irr.enh.0=1,   #irrigation enhncement factor at surface | max(1, 15.9*depth^-0.43) , depth in m
  irr.enh.L=10,  # depth till with irr.enh is constant. Then decreases to 1. 
  
  O2.w=0.276, #  #0.276     # units: [umol cm^-3]
  NO3.w=0.008,#  #0.020        # units: [umol cm^3]
  NH4.w=0.0006,  # units: [umol cm^3]
  ODU.w = 0,
  DIC.w=2.16,
  TAlk.w=0,
  rNC.1=0.15,
  rNC.2=0.15,
  rNC.3=0.15,
  rNC=0.15,
  
  por.0=0.8,           # porosity at surface units: [cm^3 pw cm^-3 bulk sed]]
  por.inf=0.8,          # porosity at depth units: [cm^3 pm cm^-3 bulk sed]]
  por.att=2,            # porosity attenuation units: [cm]
  
  depth=74,
  #temp = 6.87,
  Q10=1.88,
  S= 34,
  #w = 0.0359,          #sedimentation rate [cm/yr]
  L =20,
  N = 20,                #number of grid cell (20) 
  ps = 2.5,               #sediment denisity g/cm3
  switch.odu.burial=1,    #turn ON OFF (1 0) ODU burial mimicking pyrite burial
  #k.adj.denit = 1,      # 1 = no effect of denitrification | 0.1 = decay rate decreases by a factor of 10. Generally 1/10
  #k.adj.anoxia = 1,        # 1 = no effect of anoxia | 0.1 = decay rate decreases by a factor of 10. Generally 1/1000
  
  switch.oxy.dependent.k = 1,
  
  run_with_climatological_mean = 0, 
  run_with_climatology = 0,
  run_temporal = 1
)


grid <- setup.grid.1D(x.up=0,x.down=pars$L, N=pars$N,dx.1=0.1)
## create new nc file 

#dims
lat <- ncdim_def("Lat", "degree", yh)
lon <- ncdim_def("Lon", "degree", xh)
sed_depth <- ncdim_def("sed_depth", "cm", grid$x.mid)
diag_2D <- ncdim_def("diag_2D", "as denfined", c(1))

#vars
varOM1 <- ncvar_def("OM1","mM", list(lon,lat,sed_depth), 
                    longname = "OM1 conc in sediment", missval = NA)
varOM2 <- ncvar_def("OM2","mM", list(lon,lat,sed_depth), 
                    longname = "OM2 conc in sediment", missval = NA)
varOM3 <- ncvar_def("OM3","mM", list(lon,lat,sed_depth), 
                    longname = "OM3 conc in sediment",missval = NA)
varO2  <- ncvar_def("O2","mM", list(lon,lat,sed_depth), 
                    longname = "O2 conc in sediment",missval = NA)
varNH4 <- ncvar_def("NH4","mM", list(lon,lat,sed_depth), 
                    longname = "NH4 conc in sediment", missval = NA)
varNO3 <- ncvar_def("NO3","mM", list(lon,lat,sed_depth), 
                    longname = "NO3 conc in sediment", missval = NA)
varODU <- ncvar_def("ODU","mM", list(lon,lat,sed_depth), 
                    longname = "ODU conc in sediment", missval = NA)
varDIC <- ncvar_def("DIC","mM", list(lon,lat,sed_depth), 
                    longname = "DIC conc in sediment", missval = NA)
varTAlk <- ncvar_def("TAlk","mM", list(lon,lat,sed_depth), 
                     longname = "TAlk conc in sediment", missval = NA)
varO2_flux <- ncvar_def("O2_flux","umol cm-2 y-1", list(lon,lat,diag_2D), 
                     longname = "benthic oxygen flux", missval = NA)
varNH4_flux <- ncvar_def("NH4_flux","umol cm-2 y-1", list(lon,lat,diag_2D), 
                        longname = "benthic ammonium flux", missval = NA)
varNO3_flux <- ncvar_def("NO3_flux","umol cm-2 y-1", list(lon,lat,diag_2D), 
                        longname = "benthic nitrate flux", missval = NA)
varDIC_flux <- ncvar_def("DIC_flux","umol cm-2 y-1", list(lon,lat,diag_2D), 
                         longname = "benthic DIC flux", missval = NA)
varOM_flux <- ncvar_def("OM_flux","umol C cm-2 y-1", list(lon,lat,diag_2D), 
                         longname = "total OM deposition flux", missval = NA)
varOM_burial <- ncvar_def("OM_burial","umol C cm-2 y-1", list(lon,lat,diag_2D), 
                         longname = "OM burial", missval = NA)
vartot_denit <- ncvar_def("tot_denit","umol N cm-2 y-1", list(lon,lat,diag_2D), 
                        longname = "benthic total denitrification", missval = NA)
varAnammox <- ncvar_def("anammox","umol N cm-2 y-1", list(lon,lat,diag_2D), 
                          longname = "sediment anammox", missval = NA)



# ncss <- nc_create("ncss5.nc", list(varOM1, varOM2, varOM3,
#                                   varO2, varNH4, varNO3,
#                                   varODU, varDIC, varTAlk,
#                                   varO2_flux, varNH4_flux, varNO3_flux,
#                                   varDIC_flux, varOM_flux, varOM_burial,
#                                   vartot_denit, varAnammox))

#ncss <- nc_open("ncss1.nc", write=TRUE)
#nc_close(ncss)


# Set up parallel processing
num_cores <- detectCores() - 2  # Use all but one core
registerDoParallel(cores=num_cores)

# Create a grid of indices to process
ij_grid <- expand.grid(i=seq(1, length(xh)), j=seq(1, length(yh)))
#ij_grid <- expand.grid(i=seq(1, length(xh)), j=seq(1001, 1080))

# Parallel loop
results <- foreach(idx=1:nrow(ij_grid), .packages=c("ncdf4", "oce")) %dopar% {
  i <- ij_grid$i[idx]
  j <- ij_grid$j[idx]
  
  if (!is.na(btm_o2[i,j])) {
    #print(i);print(j)
    ## Open NetCDF file in write mode within each task
    ncss <- nc_open("ncss5.nc", write=TRUE)
    
    if (is.na((ncvar_get(ncss, 'anammox'))[i,j])) {
      
      # Copy parameters
      local_pars <- pars
      local_pars$J.OM <- max(0.0, fntot[i,j] *6.625*1000*3600*24 *365/10)
      local_pars$O2.w <- max(0.0, 1e3*btm_o2[i,j])
      local_pars$w <- max(0.0, w_tem[i,j]/(1-porosity[i,j]))
      local_pars$w <- max(0.0, w_tem[i,j]/(1-pars$por.0))
      local_pars$temp <- btm_temp[i,j]
      local_pars$por.0 <- porosity[i,j]
      local_pars$por.inf <- porosity[i,j]
      local_pars$NO3.w <- max(0.0, 1e3*btm_no3[i,j])
      local_pars$NH4.w <- max(0.0, 1e3*btm_nh4[i,j])
      local_pars$DIC.w <- max(0.0, btm_dic[i,j])
      local_pars$TAlk.w <- max(0.0, btm_talk[i,j])
      local_pars$S <- btm_salinity[i,j]
      local_pars$depth <- ocean_depth[i,j]
      
      #print(i);print(j)
      
      out.ss <- try(cbed_model(local_pars))
      if(inherits(out.ss,"try-error")){
        cat(i,j,"\n")
        #print(i);print(j)
        failed_ij <- data.frame(i,j)
        write.table(failed_ij, "~/Documents/test_ij_5.csv", append = TRUE, 
                    sep = ",", col.names = FALSE, row.names = T)
        n=0
      }
      while(inherits(out.ss, "try-error")) {
        v1 <- max(0.0, fntot[i,j] *6.625*1000*3600*24 *365/10)
        v2 <- ifelse(max(0.0, 1e3*btm_o2[i,j])<0.0001, 0.0, max(0.0, 1e3*btm_o2[i,j]))
        v3 <- max(0.0, w_tem[i,j]/(1-pars$por.0))
        v4 <- btm_temp[i,j]
        v5 <- max(0.0, 1e3*btm_no3[i,j])
        
        set.seed(n)
        tmp <- runif(5)
        local_pars$J.OM <- v1*(1+.001*tmp[1])
        local_pars$O2.w <- v2*(1+.002*tmp[2])
        local_pars$w    <- v3*(1+.001*tmp[3])
        local_pars$temp <- v4*(1+.003*tmp[4])
        local_pars$NO3.w <- v5*(1+.001*tmp[5])
        
        out.ss <- try(cbed_model(local_pars))
        n <- n + 1
      }
      
      # Open NetCDF file in write mode within each task
      #ncss <- nc_open("ncss3.nc", write=TRUE)
      
      # Write results
      ncvar_put(ncss, varOM1, out.ss$y[,"OM1"], start=c(i,j,1), count=c(1,1,20))
      ncvar_put(ncss, varOM2, out.ss$y[,"OM2"], start=c(i,j,1), count=c(1,1,20))
      ncvar_put(ncss, varOM3, out.ss$y[,"OM3"], start=c(i,j,1), count=c(1,1,20))
      ncvar_put(ncss, varO2, out.ss$y[,"O2"], start=c(i,j,1), count=c(1,1,20))
      ncvar_put(ncss, varNH4, out.ss$y[,"NH4"], start=c(i,j,1), count=c(1,1,20))
      ncvar_put(ncss, varNO3, out.ss$y[,"NO3"], start=c(i,j,1), count=c(1,1,20))
      ncvar_put(ncss, varODU, out.ss$y[,"ODU"], start=c(i,j,1), count=c(1,1,20))
      ncvar_put(ncss, varDIC, out.ss$y[,"DIC"], start=c(i,j,1), count=c(1,1,20))
      ncvar_put(ncss, varTAlk, out.ss$y[,"TAlk"], start=c(i,j,1), count=c(1,1,20))
      ncvar_put(ncss, varO2_flux, out.ss$o2_flux, start=c(i,j,1), count=c(1,1,1))
      ncvar_put(ncss, varNH4_flux, out.ss$nh4_flux, start=c(i,j,1), count=c(1,1,1))
      ncvar_put(ncss, varNO3_flux, out.ss$no3_flux, start=c(i,j,1), count=c(1,1,1))
      ncvar_put(ncss, varDIC_flux, out.ss$dic_flux, start=c(i,j,1), count=c(1,1,1))
      ncvar_put(ncss, varOM_flux, out.ss$OM_flux, start=c(i,j,1), count=c(1,1,1))
      ncvar_put(ncss, varOM_burial, out.ss$OM_burial, start=c(i,j,1), count=c(1,1,1))
      ncvar_put(ncss, vartot_denit, out.ss$tot_denit, start=c(i,j,1), count=c(1,1,1))
      ncvar_put(ncss, varAnammox, out.ss$anammox, start=c(i,j,1), count=c(1,1,1))
  }
    # Close NetCDF file
    nc_close(ncss)
  }
  NULL  # Return NULL to minimize memory usage
}

# Clean up parallel backend
stopImplicitCluster()

toc()

# Close the NetCDF file
nc_close(ncss)

stop("run complete")


#======================================================

ncss <- nc_open("ncss1.nc",write = T)

for (j in seq(1, length(yh))) {
  for (i in seq(1, length(xh))) {
    if (!is.na(btm_o2[i,j])) {
      
      pars$J.OM <- max(0.0, fntot[i,j] *6.625*1000*3600*24 *365/10)
      pars$O2.w <- max(0.0, 1e3*btm_o2[i,j])
      pars$w <- max(0.0, w[i,j])
      pars$temp <- btm_temp[i,j]
      
      out.ss <- try(cbed_model(pars))
      if(inherits(out.ss,"try-error")){
        print(cat(i,j,"\n"))
        n=0
      }
      while(inherits(out.ss,"try-error")){
        v1=max(0.0, fntot[i,j] *6.625*1000*3600*24 *365/10)
        v2=max(0.0, 1e3*btm_o2[i,j])
        v3=max(0.0, w[i,j])
        v4=btm_temp[i,j]
        
        set.seed(n)
        tmp=runif(4)
        pars$J.OM=v1*(1+.00001*tmp[1])
        pars$O2.w=v2*(1+.00001*tmp[2])
        pars$w=v3*(1+.00001*tmp[3])
        pars$temp=v4*(1+.00001*tmp[4])
        
        out.ss <- try(cbed_model(pars))
        n = n+1
      }
      
      ncvar_put(ncss, varOM1, out.ss$y[,"OM1"], start = c(i,j,1),
                count=c(1,1,20))
      ncvar_put(ncss, varOM2, out.ss$y[,"OM2"], start = c(i,j,1),
                count=c(1,1,20))
      ncvar_put(ncss, varOM3, out.ss$y[,"OM3"], start = c(i,j,1),
                count=c(1,1,20))
      ncvar_put(ncss, varO2, out.ss$y[,"O2"], start = c(i,j,1),
                count=c(1,1,20))
      ncvar_put(ncss, varNH4, out.ss$y[,"NH4"], start = c(i,j,1),
                count=c(1,1,20))
      ncvar_put(ncss, varNO3, out.ss$y[,"NO3"], start = c(i,j,1),
                count=c(1,1,20))
      ncvar_put(ncss, varODU, out.ss$y[,"ODU"], start = c(i,j,1),
                count=c(1,1,20))
      ncvar_put(ncss, varDIC, out.ss$y[,"DIC"], start = c(i,j,1),
                count=c(1,1,20))
      ncvar_put(ncss, varTAlk, out.ss$y[,"TAlk"], start = c(i,j,1),
                count=c(1,1,20))
    }
  }
}


nc_close(ncss)


# failed_ij=data.frame(i = c(i), j=c(j))
# write.csv(failed_ij, "~/Documents/test_ij_2.csv")
# 
# 
# failed_ij=data.frame(i,j)
# write.table(failed_ij, "~/Documents/test_ij_2.csv", 
#             append = TRUE, 
#             sep = ",", 
#             col.names = FALSE, 
#             row.names = T)


