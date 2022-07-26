# Comparison of emittance evolution from interpolated values from lookup tables and MADX long run in folder "sanity_checks/emit_evolution"
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import constants 
from scipy.interpolate import RegularGridInterpolator as rgi

lookup_table = pd.read_parquet('lookup_table_En_6800_Vrf_12.0_NEW.parquet')
save_interpolator_to_pickle = False

clight = constants.speed_of_light
RC     = 26658.8832
U0     = 5.945310143e-06
betar  = 1
nc     = 1
h      = 35640
etap   = 0.0003481167216

def growth_rates(dict_interp,npb,exi,eyi,bl_ns):
    try:
        dict_interp = dict_interp['interp_func']
        Nbkeys = np.array(list(dict_interp.keys()))
        argn = [0]
        taux_h = (dict_interp[Nbkeys[argn[0]]]['txh'])((exi,eyi,bl_ns))*Nbkeys[argn[0]]/npb
        tauy_h = (dict_interp[Nbkeys[argn[0]]]['tyh'])((exi,eyi,bl_ns))*Nbkeys[argn[0]]/npb
        taul_h = (dict_interp[Nbkeys[argn[0]]]['tlh'])((exi,eyi,bl_ns))*Nbkeys[argn[0]]/npb
    except:
        print('Check if out of range of interpolated values')
        taux_h = tauy_h = taul_h = np.nan
    return taux_h, tauy_h, taul_h

def compute_evol(dict_interp, npb,exi,eyi,bl_ns,dt,t_fill): 

    ex0  = dict_interp['equilibrium_properties']['ex0']
    ey0  = dict_interp['equilibrium_properties']['ey0']
    sp0  = dict_interp['equilibrium_properties']['sp0']
    ss0  = dict_interp['equilibrium_properties']['ss0']
    taux = dict_interp['equilibrium_properties']["taux"]
    tauy = dict_interp['equilibrium_properties']["tauy"]
    taul = dict_interp['equilibrium_properties']["taul"]
    En   = dict_interp['equilibrium_properties']["En"]
    V0   = dict_interp['equilibrium_properties']["V0"]

    time = [0]
    ex = np.array(exi)
    ey = np.array(eyi)
    bl = np.array(bl_ns)
    npb = np.array(npb)
    bl_ns = np.array(bl_ns)
    tt = 0
    bli = bl_ns*clight/4/1e9
    dpp=np.sqrt(2/np.pi)*np.sqrt(nc*(V0-U0)*(np.sin(bli*h*np.pi*betar/RC))**2)/np.sqrt(En*h*np.abs(etap))/betar;
    while tt < t_fill+dt:
        #print(tt,exi, eyi, bl_ns)
        taux_h, tauy_h, taul_h = growth_rates(dict_interp,npb,exi,eyi,bl_ns)        
        #print(taux_h, tauy_h, taul_h)
        if tt == 0:
            txh = taux_h; tyh = tauy_h; tlh = taul_h;
        else:
            txh = np.vstack((txh, taux_h))
            tyh = np.vstack((tyh, tauy_h))
            tlh = np.vstack((tlh, taul_h))

        Tx = 1./(taux_h*3600.)/2.
        Ty = 1./(tauy_h*3600.)/2.
        Tl = 1./(taul_h*3600.)/2.
        exi = (-ex0+np.exp(2*dt*(Tx-1/taux))*(ex0+exi*(-1+Tx*taux)))/(-1+Tx*taux)
        eyi = (-ey0+np.exp(2*dt*(Ty-1/tauy))*(ey0+eyi*(-1+Ty*tauy)))/(-1+Ty*tauy)
        dpp = (-sp0+np.exp(dt*(Tl-1/taul))*(sp0+dpp*(-1+Tl*taul)))/(-1+Tl*taul)
        bli = clight*RC*np.arccos((En*nc*(V0-U0)*betar**2-dpp**2*En**2*h*np.pi*betar**4*abs(etap))/(En*nc*(V0-U0)*betar**2))/(2*clight*h*np.pi*betar);
        bl_ns = bli/clight*4*1e9
        tt = tt + dt
        ex = np.vstack((ex,exi))
        ey = np.vstack((ey,eyi))
        bl = np.vstack((bl, bl_ns))
        time.append(tt)
    taux_h, tauy_h, taul_h = growth_rates(dict_interp,npb,exi,eyi,bl_ns)
    txh = np.vstack((txh, taux_h))
    tyh = np.vstack((tyh, tauy_h))
    tlh = np.vstack((tlh, taul_h))
    return time, ex, ey, bl, txh, tyh, tlh


exun = np.unique(lookup_table["exin"])
eyun = np.unique(lookup_table["eyin"])
blun = np.unique(lookup_table["bl_ns"])
Nbun = np.unique(lookup_table["npbb"])

dict_interp = {'interp_func': {},
               'equilibrium_properties': {
               'ex0'  : lookup_table['ex0'][0].iloc[0],
               'ey0'  : lookup_table['ey0'][0].iloc[0],
               'sp0'  : lookup_table['sp0'][0].iloc[0],
               'ss0'  : lookup_table['ss0'][0].iloc[0],
               'taux' : lookup_table["taux"][0].iloc[0],
               'tauy' : lookup_table["tauy"][0].iloc[0],
               'taul' : lookup_table["taul"][0].iloc[0],
               'En'   : lookup_table["en"][0].iloc[0],
               'V0'   : lookup_table["v0"][0].iloc[0]}
              }

for nn,Nbb in enumerate(Nbun):
    dict_interp['interp_func'][Nbb] = {}
    df_masked = lookup_table[lookup_table['npbb'] == Nbb]

    Tx = np.zeros((len(exun), len(eyun), len(blun)))
    Ty = np.zeros((len(exun), len(eyun), len(blun)))
    Tl = np.zeros((len(exun), len(eyun), len(blun)))
    print(f"Loop over intensity {nn+1}/{len(Nbun)}")

    for ii,ex in enumerate(exun):
        for jj,ey in enumerate(eyun): 
            for kk,bl in enumerate(blun):
                mask = (df_masked['exin'] == ex) & (df_masked['bl_ns'] == bl) & (df_masked['eyin'] == ey)
                Tx[ii,jj,kk] = df_masked['txh'][mask][0]
                Ty[ii,jj,kk] = df_masked['tyh'][mask][0]
                Tl[ii,jj,kk] = df_masked['tlh'][mask][0]

    fnx = rgi((exun,eyun,blun),Tx)
    fny = rgi((exun,eyun,blun),Ty)
    fnl = rgi((exun,eyun,blun),Tl)
    
    dict_interp['interp_func'][Nbb]['txh'] = fnx
    dict_interp['interp_func'][Nbb]['tyh'] = fny
    dict_interp['interp_func'][Nbb]['tlh'] = fnl
    
if save_interpolator_to_pickle:
    import pickle
    with open('interpolator_vs_nb_En6800_VRF12.0.pkl', 'wb') as f:
        pickle.dump(dict_interp, f)

# Sanity check
df_test = pd.read_parquet("sanity_checks/emit_evolution/IBS_long_output.parquet")

ex_ref   = [1.5, 1.0]
ey_ref   = [1.3, 1.8]
npb_ref  = [1.08e11, 1.2e11]
blns_ref = [1.1, 0.9]

dt = 60 
t_fill = 5*60*60-2*60.
time, ex, ey, bl, txh, tyh, tlh =  compute_evol(dict_interp, npb_ref,ex_ref,ey_ref,blns_ref,dt,t_fill)

fig, ax = plt.subplots(nrows=3, figsize=(8,10))
for key, group in df_test.groupby("npbb"):
    plt.sca(ax[0])
    plt.plot(group.tt/3600., group.exin*1e6, c='lime', lw=5)
    plt.sca(ax[1])
    plt.plot(group.tt/3600., group.eyin*1e6, c='lime', lw=5)
    plt.sca(ax[2])
    plt.plot(group.tt/3600., group.bl_ns, c='lime', lw=5)
    
plt.sca(ax[0])
plt.plot(np.array(time)/3600., ex, lw=3)
plt.ylabel("emit_x")
plt.sca(ax[1])
plt.plot(np.array(time)/3600., ey, lw=3)
plt.ylabel("emit_y")
plt.sca(ax[2])
plt.plot(np.array(time)/3600., bl, lw=3)
plt.ylabel("bl")
plt.xlabel("Time (h)")
#fig.savefig('pngs/1.png')

fig, ax = plt.subplots(nrows=3, figsize=(8,10))
counter=0
for key, group in df_test.groupby("npbb"):
    plt.sca(ax[0])
    plt.plot(group.tt/3600., (group.exin*1e6-ex[:,counter])/ex[:,counter]*100., c='b',lw=2)
    plt.ylabel("emit_x madx-model/model (%)")
    plt.sca(ax[1])
    plt.plot(group.tt/3600., (group.eyin*1e6-ey[:,counter])/ey[:,counter]*100., c='b', lw=2)
    plt.ylabel("emit_y madx-model/model (%)")
    plt.sca(ax[2])
    plt.plot(group.tt/3600., (group.bl_ns-bl[:,counter])/bl[:,counter]*100., c='b', lw=2)
    plt.ylabel("bl madx-model/model (%)")
    counter+=1
#fig.savefig('pngs/2.png')
