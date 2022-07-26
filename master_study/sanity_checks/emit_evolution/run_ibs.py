from cpymad.madx import Madx
import yaml
import tree_maker

# Original script from https://github.com/stpapado/calculateIBS_MADX
# Cleanup January 2022

with open("config.yaml", "r") as f:
    config = yaml.safe_load(f)

n_studies = len(config["tflattopconst"])
results = []

for current_study in range(n_studies):

  config["egrx_extra"] *= float(config["dt"][current_study])/3600*1e-6
  config["egry_extra"] *= float(config["dt"][current_study])/3600*1e-6
  
  tree_maker.tag_json.tag_it(config['log_file'], 'started')
  
  mad = Madx()
  
  mad.input('''
  !********************************************************
  ! IBS calculations for the LHC at injection
  ! 
  ! F. Antoniou
  ! June 2014
  ! *******************************************************
  
  ! rename tfs file
  renamefile(filename,save_to):macro={
    system,"mv filename.tfs save_to";
  }
  ''')
  
  mad.input(f'''
  
  ! create beam for protons
  beam,particle=proton,pc={config["En"][current_study]};
  
  ! total duration in seconds
  tmax={config["tflattopconst"][current_study]}*60+1;
  value,tmax;
  
  ''')
  
  mad.input(f'''
  call,file="{config["sequence"][current_study]}";
  
  
  ! twiss
  use,    period=lhcb1; 
  select, flag=twiss, column=name,s,betx,bety,alfx,alfy,dx,dpx,dy,dpy,mux,muy,x,px,y,py;
  !twiss, sequence=lhcb1,chrom,centre, file=twiss.all.b1.data;
  twiss, sequence=lhcb1,chrom,centre;
  
  
  ! constants
  mp=938.2723*1e-3; ! [GeV/c^2]
  rp=1.534698e-18; ![m]
  me=0.51099906*1e-3;
  re=2.8179409e-18; 
  
  ! revolution frequency
  fr=11.2455e3;
  
  ! Constant Cgg  =4pi / 3 * rp/(mp*c^2)^3 for energy loss per turn due to synch rad
  Cgg=4*Pi/3*rp/(mp)^3;
  
  ! Set correct units for intensity & emittance
  Npbtest=1.E+11;
  exn_i=1.E-6;
  eyn_i=1.E-6;
  
  ! bunch length in m from ns
  bl_i={config["blns"][current_study]}*1.E-9*clight/4;
  
  i=0;
  
  ! rf voltage
  V0={config["V0max"][current_study]};
  
  nc=1;
  ! length of accelerator
  RC=table(summ,length);
  
  ! time for 1 revolution
  T0=RC/clight;
  
  h=hrf400;
  show,h;
  
  ! gtransition, gamma, beta, energy
  gt=table(summ,gammatr);
  gamma=beam->gamma;
  betar=beam->beta;
  En=beam->energy;
  E0=En/gamma;
  show,betar;
  
  ! slip factor abs(1/gtrans^2 - 1/gamma^2)
  etap=abs(1/gt^2-1/gamma^2);
  show,etap,En,bl_i,RC,betar,h,T0,V0;
  
  dpp=sqrt(2/pi)*sqrt(nc*(V0-U0)*(sin(bl_i*h*pi*betar/RC))^2)/sqrt(En*h*abs(etap))/betar;
  !el = pi*dpp*4*bl_i/clight*En*1e9; 
  el = 4*pi*dpp*bl_i/clight*En*1e9; ! http://digitalcommons.calpoly.edu/cgi/viewcontent.cgi?article=1441&context=phy_fac
  show,gt,el,dpp,bl_i,V0,En,case;
  !stop;
  
  ! geometric emittance
  exi=exn_i/gamma;
  eyi=eyn_i/gamma;
  
  ! steps duration/time step
  steps=tmax/{config["dt"][current_study]};
  
  tt:=i*{config["dt"][current_study]};
  Npb={config["k"][current_study]}/100*Npbtest;
  Npb0={config["k"][current_study]}/100*Npbtest;
  exi={config["jj"][current_study]}/100*exi;
  eyi={config["jjy"][current_study]}/100*eyi;
  vrf400=V0;
  
  ! initialize beam parameters such as geom emittance, intensity, energy spread, bunch length
  beam,npart=Npb;
  beam,EX=exi;
  beam,EY=eyi;
  beam,sigt=bl_i;
  beam,sige=dpp;
  
  ! synchrotron radiation considered in all dipole magnets
  beam,radiate;
  
  ! plot beta functions & dispersion
  !PLOT, NOTITLE=TRUE, COLOUR=100, HAXIS=S,VAXIS1=BETX,BETY,VAXIS2=DX,DY,interpolate=true;
  !PLOT, NOTITLE=TRUE, COLOUR=100, HAXIS=S,VAXIS1=DX,DY,interpolate=true;
  
  ! when beam, radiate, cmp equilibrium emittances 
  emit,deltap=0;
  show,beam;
  
  !equilibrium emittance
  ex0=beam->ex;
  ey0=beam->ey;
  sp0=beam->sige;
  ss0=beam->sigt;
  En=beam->energy;
  
  show, ex0, ey0, sp0, ss0;
  value, ex0*gamma, ey0*gamma, sp0, ss0;
  
  ! synchrotron radiation integrals 
  twiss,chrom;
  
  I1=table(summ,synch_1);
  I2=table(summ,synch_2);
  I3=table(summ,synch_3);
  I4=table(summ,synch_4);
  !I4=2e-4;
  I5=table(summ,synch_5);
  
  if ({config["flag_SR"]}==1){{
  ! energy loss U0 in one revolution from synchrotron radiation, Uo=1/2pi C E^4 I2 https://indico.cern.ch/event/279729/contributions/1626389/attachments/512375/707123/ElectronDynamicsLRivkin.pdf
  !  U0=Cgg/(2*Pi)*(En)^4*I2;	![GeV/turn]
    U0=Cgg/2/Pi*En^4*I2;
    ! https://arxiv.org/pdf/1507.02213.pdf
    ! damping time tau = 2 * En/Uo *T0 / j where j is damping partition number
    taux=2*En*(RC/clight)/U0;	! Beam size damping time
    tauy=2*En*(RC/clight)/U0;	! Beam size damping time
    taul=2*En*(RC/clight)/U0/2;	! Bunch length damping time
    value, taux, tauy, taul, En, RC, clight, U0, Cgg, I2;
    !exit;
  }}
  else{{
   U0=1.e-15;
   taux=1.e-15;
   tauy=1.e-15;
   taul=1.e-15;
  }}
  
  ! partition numbers
  !Jx=1-I4/I2;
  !Jy=1-I4/I2;
  !Js=2+I4/I2;
  
  ! ez(t) = ez(t0) * exp(-2 *t /taul)
  
  hb=6.582122e-16;
  Cq=55/32/sqrt(3)*hb*clight/(mp*1e9);
  !ex0=Cq*gamma^3/Jx*I5/I2;
  value,I1,I2,I3,I4,I5,U0,En,Cgg,taux/3600,tauy/3600,taul/3600,gamma,ex0*gamma,Cq;
  
  ss0=clight*RC*acos(1-(sp0^2*En*h*pi*betar^2*abs(etap))/(nc*(V0-U0)))/(2*clight*h*pi*betar);
  Qs=sqrt(h*etap/2/pi/betar^2/En*V0*cos(asin(U0/V0)));
  
  show,spp0,ss0,Qs,h,V0,U0;
  !exit;
  beam,radiate=false;
  vrf400=0;
  
  ! exi is the initial emittance
  beam,EX=exi;
  beam,EY=eyi;
  beam,sigt=bl_i;
  beam,sige=dpp;
  beam,npart=Npb;
  show,beam;
  show,exi,eyi,bl_i,dpp;
  !stop;
  exit0=exi;
  eyit0=eyi;
  blit0=bl_i;
  dpit0=dpp;
  
  create,table=ibsscan,column=tt,Npbb,ex0,ey0,sp0,ss0,exin,eyin,dpp,bl_ns,En,Txh,Tyh,Tlh,V0,taux,tauy,taul,L0,tau_Boff;
  
  exIBS:=exi/exit0;
  eyIBS:=eyi/eyit0;
  dppIBS:=dpp/dpit0;
  blIBS:=bl_i/blit0;
  
  exin:=exi*gamma;
  eyin:=eyi*gamma;
  bl_ns:=bl_i/clight*4*1e9;
  
  twiss,centre;
  
  steps=tmax/{config["dt"][current_study]};
  
  Npbb=Npb;
  
  sxst:=sqrt(exi*{config["betastar"]});
  syst:=sqrt(eyi*{config["betastar"]});
  
  Fgeom:=1/(sqrt(1+(bl_i*{config["phi"]}/sqrt(exi*{config["betastar"]})/2)^2));
  Lfact:=fr*{config["nb"]}*(Npbb^2)/(sqrt(2*sxst^2)*sqrt(2*syst^2))/2/pi;
  
  tau_Boff:=Npbb/({config["nIPs"]}*{config["sigmaBOff"]}*L0);
  !value,fr,{config["nb"]},Npbb,sxst,syst,exi,eyi,{config["betastar"]};
  !stop;
  
  
  !ibs;
  !mpla=ibs.tx/3600;
  !value ,sp0,mpla;
  !stop;
  
  while(tt<tmax-1)
    {{        
      L0=Lfact*Fgeom;
      
      if({config["flag_IBS"]}==1)
      {{
        ibs;
        
        ! Beam Size growth times 
        Tx=1/ibs.tx/2;
        Ty=1/ibs.ty/2;
        Tl=1/ibs.tl/2;
  	      	  
        ! Emittance growth times 
        Txh=(1/Tx)/3600/2;	
        Tyh=(1/Ty)/3600/2;
        Tlh=(1/Tl)/3600/2;
      }}    
      else
      {{
        Tx=1e-15;
        Ty=1e-15;
        Tl=1e-15;
        ex0=1e-15;
        ey0=1e-15;
        sp0=1e-15;
        ss0=1e-15;
        
        Txh=(1/Tx)/3600;
        Tyh=(1/Ty)/3600;
        Tlh=(1/Tl)/3600;
      }}
    
      fill,table=ibsscan;
    
      if({config["flag_SR"]}==1)
      {{
        exi=(-ex0+exp(2*{config["dt"][current_study]}*(Tx-1/taux))*(ex0+exi*(-1+Tx*taux)))/(-1+Tx*taux);
        eyi=(-ey0+exp(2*{config["dt"][current_study]}*(Ty-1/tauy))*(ey0+eyi*(-1+Ty*tauy)))/(-1+Ty*tauy);
        dpp=(-sp0+exp({config["dt"][current_study]}*(Tl-1/taul))*(sp0+dpp*(-1+Tl*taul)))/(-1+Tl*taul);
  
      }}    
      else
      {{
        exi=exi*exp(2*{config["dt"][current_study]}*Tx);
        eyi=eyi*exp(2*{config["dt"][current_study]}*Ty);
        dpp=dpp*exp({config["dt"][current_study]}*Tl);
        U0=0;
      }}
      
      bl_i=clight*RC*acos((En*nc*(V0-U0)*betar^2-dpp^2*En^2*h*pi*betar^4*abs(etap))/(En*nc*(V0-U0)*betar^2))/(2*clight*h*pi*betar);
      show, bl_i;
  
      if(bl_i < {config["bl_lev"]}*1e-9*clight/4)
      {{
  	bl_i = {config["bl_lev_max"]}*1e-9*clight/4;
  	dpp=sqrt(2/pi)*sqrt(nc*(V0-U0)*(sin(bl_i*h*pi*betar/RC))^2)/sqrt(En*h*abs(etap))/betar;	
      }}
  
     if({config["flag_extraGrowth"]}==1)
      {{   
  	value, {config["egrx_extra"]},{config["egry_extra"]};
       	exi=exi+{config["egrx_extra"]}/gamma;
       	eyi=eyi+{config["egry_extra"]}/gamma;
      }}
  
  
      beam,EX=exi;
      beam,EY=eyi;
      beam,sigt=bl_i;
      beam,sige=dpp;
      beam,pc=En-E0;
      beam,energy=En;
      beam,gamma=gamma;
      beam,npart=Npbb;
      
      if ({config["flag_BOff"]}==1)
      {{
        Npbb=Npbb/(1+{config["dt"][current_study]}/tau_Boff);
      }}
      show,beam;
  
      value,tt,tau_Boff,Npbb;
      i=i+1;
  }}; // end of timescan
  value,tt;
  el2 = 3.14*dpp*4*bl_i/3e8*En*1e9;
  value,exIBS,eyIBS,blIBS,exi*gamma,eyi*gamma,dpp,bl_i,el2,Txh,Tyh,Tlh,taux/3600,tauy/3600,taul/3600;
  
  !stop;
  blns={config["blns"][current_study]}*100;
  V0={config["V0max"][current_study]}*1e3;
  
  show, V0, blns;
  
  ''')
  
  
  df = mad.table['ibsscan'].dframe(index=mad.table['ibsscan'].tt)
  df["emitx_init"] = df.exin.iloc[0]*1e6
  df["emity_init"] = df.eyin.iloc[0]*1e6
  df["blns_init"]  = df.bl_ns.iloc[0]
  df["En"]         = df.en.iloc[0]
  results.append(df)
  
  #if config["tfs_or_parquet"] == 'tfs':
  #  mad.input(f'''
  #  write,table=ibsscan,file="ibsscan_lhc_he_lev_opt_inj.tfs";
  #  exec,renamefile(ibsscan_lhc_he_lev_opt_inj,{config["save_to"]}_{current_study}.tfs);
  #  
  #  return;
  #  ''')
  #elif config["tfs_or_parquet"] == 'parquet':
  #    mad.table['ibsscan'].dframe(index=mad.table['ibsscan'].tt).to_parquet(f'{config["save_to"]}')
  #else:
  #    print('Wrong format specified for saving results')
import pandas as pd
results = pd.concat(results, axis=0)
results.to_parquet(f'{config["save_to"]}')
tree_maker.tag_json.tag_it(config['log_file'], 'completed')
