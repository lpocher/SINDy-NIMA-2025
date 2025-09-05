import sys
# put path for upper most directory for test for .gs and current files etc.
sys.path.append('/home/lpocher/PySINDy/test')
from warp import *
from beam import createBeam
from settings import params
from util import *
import numpy as np
from tune_util import naff
import sys
# from mphoto import take_photo


setup()#,prefix='cgmfiles/umerlat_main')
# winon()
#def setup():
#    pass

# --- Set four-character run id, comment lines, user's name.
top.pline2   = "Main model of umer lattice"
top.pline1   = "UMER LATTICE Centroid Oscillations"
top.runmaker = "Liam A. Pocher"

#### ---- plotting setting for nice plots  ----- #######
# make nice plots
r_array = range(0,240,1)
b_array = range(240,0,-1)
g_array = np.zeros(240,dtype='uint8')

# --- Locations to save simulation image data (i.e. photos)
# make sure the folder exists of course
FOLDER_LOC = './data/' # folder to save data
if not os.path.exists(FOLDER_LOC):
   os.makedirs(FOLDER_LOC)
fnameid =  sys.argv[0][:sys.argv[0].find(".py")]
FOLDER_LOC_ENV = FOLDER_LOC # folder to save envelope data
fnameid_ENV = 'env' # file name
runEnvEqn = False # run the envelope solver (i.e. not PIC code, just integrates moment equations)

# set values for envelope equation to run
w3d.zmmax = 1 # setting so mphoto.py can run when called in the runEnvEqn script

# --- lattice settings
turns = 300 # number of turns
addearthfield = 0 # include earth field or not
YSECTION = 0 # include Y section or not (replaced with periodic/perfect cell)
numCells = 36 - YSECTION # number of cells, need 36 for 1+ turns. Can run 0 turns and 1 cell ex.

# correction factors for single particle to warp matching
# --- Magnet settings (Amps)
QF_current = 1.826#1.790
QD_current = 1.826#1.816
YQ_current = 5.0
QR1_current = 5.0

# create default umer beam
beam = createBeam(aperture='pencil')
# +++ envelope and slope [m, rad]
top.a0 = beam.a0
top.b0 = beam.b0
top.ap0 = beam.ap0
top.bp0 = beam.bp0
# nocurrent
#top.a0       = 0.0015704195584892587#beam.a0
#top.b0       = 0.0015049406560583435#beam.b0
#top.ap0      = -0.006680161870713737#beam.ap0
#top.bp0      = 0.006596211762159486#beam.bp0
# 0.6mA
top.a0       = 1.25*0.001736182188049734#beam.a0
top.b0       = 0.0016195274771176325#beam.b0
top.ap0      = -0.007277374924081305#beam.ap0
top.bp0      = 0.006893341383986165#beam.bp0
# 6mA
#top.a0       = 1.25*0.003130036371720845
#top.b0       = 0.0027222311657195447
#top.ap0      = -0.012143906683727324
#top.bp0      = 0.010666868163966354

# 0.6 match w/ dipoles
#0.001704633026099904
#0.0016244854649150207
#-0.00705769416493017
#0.006889042456371745

# +++ centroid and angle [m, rad]
top.xcent_s  = 0.982395*1.0e-3#beam.xcent
top.xpcent_s = 0.0#beam.xpcent
top.ycent_s  = 1.0173*1.0e-3
top.ypcent_s = 0.0

# +++ current [A], unnormalized 4*rms emittance [m-rad]
top.ibeam    = beam.ibeam
top.emitx    = beam.emit
top.emity    = beam.emit
top.emit     = top.emit

# boundaries
piperad = 1.0e-2 # pipe radius
top.prwall = piperad

# +++ Set input parameters describing the 3d simulation.
w3d.nx = 64
w3d.ny = 64
w3d.nz = 16
step_size = 5.0e-3
top.dt = (step_size)/beam.vbeam
top.prwall = piperad
top.nhist = 1
top.lrelativ  = True

# --- Set to finite beam.
w3d.xmmin = -piperad
w3d.xmmax =  piperad
w3d.ymmin = -piperad
w3d.ymmax =  piperad

# --- Load Semi-Gaussian cigar beam.
top.npmax = 10000
w3d.distrbtn = "semigaus"

top.lhxrmsz = true
top.lhyrmsz = true

# load dipole setpoints
d_vals = np.append(np.genfromtxt('currents/dvals_06mA.csv'),14.3)
d_vals = np.append(np.genfromtxt('currents/dvals_06mA.csv'),2.5)
 #-----------------------------------------------------------------------------------------------------------------
# -- Build lattice

# standard magnet lengths
D_angle = 10*np.pi/180 # radians
rBend_l_eff = params['D'][3]*params['D'][1]
Q_l_eff_D = params['QD'][3]*params['QD'][1]
Q_l_eff_F = params['QF'][3]*params['QF'][1]
# Standard cell drift lengths
DC_fodo_1 = .08 - Q_l_eff_D/2
DC_fodo_2 = .08 - Q_l_eff_D/2 - rBend_l_eff/2
DC_fodo_3 = .08 - Q_l_eff_F/2 - rBend_l_eff/2
DC_fodo_4 = .08 - Q_l_eff_F/2
# YSECTION 
PD_angle = 6*np.pi/180 # radians
YQ_l_eff = params['YQ'][3]*params['YQ'][1]
QR_l_eff = params['QR1'][3]*params['QR1'][1]
PDRec_l_eff = params['PD'][3]*params['PD'][1]
DC_YQ = .08 - YQ_l_eff/2;
DC_YQ_to_PD = .08 - PDRec_l_eff/2 - YQ_l_eff/2;
DC_PD_to_QR1 = .08 - PDRec_l_eff/2 - QR_l_eff/2;
DC_QR1 = .08 - QR_l_eff/2;

# standard period/cell lattice
top.diposet = False
zloc = 0

for turn in range(turns):

    print("zloc = ", zloc)

    for i in range(numCells):
        addnewdrft(zs=zloc, ze=zloc+DC_fodo_1)
        zloc+=DC_fodo_1
        addnewquad(zs=zloc, ze=zloc+Q_l_eff_D, db=+QD_current*params['QD'][2]*params['QD'][0])
        addnewhele(zs=zloc,ze=zloc+Q_l_eff_D, nn=[3,4], am=[-0.1,-0.05])
        zloc+=Q_l_eff_D
        addnewdrft(zs=zloc, ze=zloc+DC_fodo_2)
        zloc+=DC_fodo_2
        # ------- special dipole stuff -----------------------
        new_angle = c2angle('D'+str(i+1),d_vals[i],efield=0)
        addnewdipo(zs=zloc, ze=zloc+rBend_l_eff, by=-1*params['D'][2]*params['D'][0]*d_vals[i], 
            ta=np.tan(new_angle/2), tb=np.tan(-1*new_angle/2))
        addnewbend(zs=zloc, ze=zloc+rBend_l_eff, rc=rBend_l_eff/new_angle)
	#addnewdrft(zs=zloc, ze=zloc+rBend_l_eff)
        # ------- --------------------------------------------
        zloc+=rBend_l_eff
        addnewdrft(zs=zloc, ze=zloc+DC_fodo_3)
        zloc+=DC_fodo_3
        addnewquad(zs=zloc, ze=zloc+Q_l_eff_F, db=-QF_current*params['QF'][2]*params['QF'][0])
        zloc+=Q_l_eff_F
        addnewdrft(zs=zloc, ze=zloc+DC_fodo_4)
        zloc+=DC_fodo_4
    if YSECTION:
        addnewdrft(zs=zloc, ze=zloc+DC_YQ)
        zloc+=DC_YQ
        addnewquad(zs=zloc, ze=zloc+YQ_l_eff, db=YQ_current*params['YQ'][2]*params['YQ'][0],ph=0.01745*1.0)
	#addnewhele(zs=zloc,ze=zloc+YQ_l_eff, nn=[3],am=[-1.0])
        zloc+=YQ_l_eff
        addnewdrft(zs=zloc, ze=zloc+DC_YQ_to_PD)
        zloc+=DC_YQ_to_PD
        # -------- special PD dipole stuff
        new_angle = c2angle('PD',d_vals[-1],efield=0)
        addnewdipo(zs=zloc, ze=zloc+PDRec_l_eff, by=-1*params['PD'][2]*params['PD'][0]*d_vals[-1], 
            ta=np.tan(new_angle/2), tb=np.tan(-1*new_angle/2))
        addnewbend(zs=zloc, ze=zloc+PDRec_l_eff, rc=PDRec_l_eff/new_angle)
        # ------------------------------
        zloc+=PDRec_l_eff
        addnewdrft(zs=zloc, ze=zloc+DC_PD_to_QR1)
        zloc+=DC_PD_to_QR1
        addnewquad(zs=zloc, ze=zloc+QR_l_eff, db=-QR1_current*params['QR1'][2]*params['QR1'][0])
        zloc+=QR_l_eff
        addnewdrft(zs=zloc, ze=zloc+DC_QR1)
        zloc+=DC_QR1
"""
# basic plots
plotfreq = 5
@callfromafterstep
def runtimeplots():
    # --- Make plots of the transverse distribution in two zwindows.
    if top.it % plotfreq == 0:
        plsys(3)
        beam.ppxy()
        limits(-0.015, +0.015, -0.015, +0.015)
        plsys(4)
        beam.ppxxp(color='density')
        plsys(5)
        beam.ppxy(color='density')
        limits(-0.015, +0.015, -0.015, +0.015)
        plsys(6)
        beam.ppyyp(color='density')
        fma()
        pptrace(color='density')
"""

fma()
# env 
if runEnvEqn:
    env.zl = 0
    env.zu = zloc
    env.dzenv = .005
    package("env")
    generate()
    step()
    penv()
    fma()
    np.savetxt(FOLDER_LOC_ENV + fnameid_env + '.txt',np.array([env.zenv,env.aenv,env.benv]).T)
    if 0:
        from warp import env_match as em
        em.match1(n=100)
        step()
    stophere


w3d.zmmax = 0.1
package("wxy")
generate()
w3d.zmmax = 0.1
steps = int(zloc/step_size)

palette(r_array, g_array, b_array)
top.ncolor = 10

stime = wtime()
if 0:
	step(steps)
else:
    x,y,xp,yp = np.zeros((turns,top.npmax)),np.zeros((turns,top.npmax)),np.zeros((turns,top.npmax)),np.zeros((turns,top.npmax))
    for i in range(turns):
        for j in range(top.npmax):
            x[i,j] = top.pgroup.xp[j]
            y[i,j] = top.pgroup.yp[j]
            xp[i,j] = top.pgroup.uxp[j]
            yp[i,j] = top.pgroup.uyp[j]
        step(int(zloc/turns/step_size))

etime = wtime()

print("Time Running = ", etime-stime)
print("Per a turn = ", (etime-stime)/turns)

steps_per_turn = int(zloc/turns/step_size)

hpxbar(color='blue')
hpybar(color='red',linetype='dash')
fma()
hpepsx(color='blue')
hpepsy(color='red',linetype='dash')
fma()
hppnum()
fma()

# ---- getting index for saving things
idx = np.arange(0,turns*steps_per_turn+1).astype(int)

# --- save Poincare sections at each turn (x,y,ux,uy) for the beam
np.savetxt(FOLDER_LOC+fnameid+'X.txt',x)
np.savetxt(FOLDER_LOC+fnameid+'Y.txt',y)
np.savetxt(FOLDER_LOC+fnameid+'Xp.txt',xp)
np.savetxt(FOLDER_LOC+fnameid+'Yp.txt',yp)
# --- save centroid values in z, x_c, y_c
np.savetxt(FOLDER_LOC+fnameid+'centroid_x_y.txt',np.array([top.hzbar[0,idx,0], top.hxbar[0,idx,0], top.hybar[0,idx,0]]).T)
np.savetxt(FOLDER_LOC+fnameid+'centroid_xp_yp.txt',np.array([top.hzbar[0,idx,0], top.hxpbar[0,idx,0], top.hypbar[0,idx,0]]).T)
# --- save sigma matrix values sequentially
np.savetxt(FOLDER_LOC + 'sigma_xx.txt',np.array([top.hzbar[0,idx,0],top.hxsqbar[0,idx,0]]).T)
np.savetxt(FOLDER_LOC + 'sigma_xxp.txt',np.array([top.hzbar[0,idx,0],top.hxxpbar[0,idx,0]]).T)
np.savetxt(FOLDER_LOC + 'sigma_xy.txt',np.array([top.hzbar[0,idx,0],top.hxybar[0,idx,0]]).T)
np.savetxt(FOLDER_LOC + 'sigma_xyp.txt',np.array([top.hzbar[0,idx,0],top.hxypbar[0,idx,0]]).T)
np.savetxt(FOLDER_LOC + 'sigma_xpxp.txt',np.array([top.hzbar[0,idx,0],top.hxpsqbar[0,idx,0]]).T)
np.savetxt(FOLDER_LOC + 'sigma_yxp.txt',np.array([top.hzbar[0,idx,0],top.hyxpbar[0,idx,0]]).T)
np.savetxt(FOLDER_LOC + 'sigma_xpyp.txt',np.array([top.hzbar[0,idx,0],top.hxpypbar[0,idx,0]]).T)
np.savetxt(FOLDER_LOC + 'sigma_yy.txt',np.array([top.hzbar[0,idx,0],top.hysqbar[0,idx,0]]).T)
np.savetxt(FOLDER_LOC + 'sigma_yyp.txt',np.array([top.hzbar[0,idx,0],top.hyypbar[0,idx,0]]).T)
np.savetxt(FOLDER_LOC + 'sigma_ypyp.txt',np.array([top.hzbar[0,idx,0],top.hypsqbar[0,idx,0]]).T)
np.savetxt(FOLDER_LOC + 'eps_y.txt',np.array([top.hzbar[0,idx,0],top.hepsx[0,idx,0]]).T)
np.savetxt(FOLDER_LOC + 'eps_x.txt',np.array([top.hzbar[0,idx,0],top.hepsy[0,idx,0]]).T)


# ---- print number of partilces leaving the beam
print("Number of particles leaving beam: " + str(top.npmax - top.nplive))

# ---- Leaving the simulation Python Style
raise SystemExit

