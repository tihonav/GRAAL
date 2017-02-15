import sys
import os
import math
import numpy as np
import time
import cPickle 
import datetime
from itertools import izip
#
# Usage example:
#     python $DMPSWSYS/../Analysis/Stk/python/COSM_APR2015/AlignStk.py /dampe/data1/cosmics_dpnc_april_2015/REC2/*root  --bdchnl=$DMPSWSYS/../Calibration/ParametersSTK_FM/bad_chan.txt --align=$DMPSWSYS/../Calibration/ParametersSTK_FM/stk_alignment_step0.txt
#
#


printhelp = False
for arg in sys.argv:
    if arg =="-h" or arg =="-H" or arg[1:].lower() == "help" or arg[2:].lower()=="help" or arg.lower()=="help" or arg.lower()=="h": 
        printhelp = True 
        sys.argv.remove(arg)
    

from ROOT import *
gSystem.Load("libDmpEvent.so")
gROOT.SetBatch(kTRUE)
#gErrorIgnoreLevel = kFatal

MICRON   = 0.001   # all the distances in the ROOT file are in the units of mm
MICRORAD = 0.000001
N_DOF_LADDER  = 10
MAX_POSSIBLE_CLUSTERS = 6
NLADDERS_TRB = 24
NSENSORS  = 4
NAN       = -1.e11

OUTPATH_ENV_NAME = "DMPOUTPATH"

ANGLES = [0.2, 0.4, 0.6, 0.8,1.2, 1.6,  999.]
ANGLES_FORALGORITHM = [999.]  # [0.4, 999.]
#ANGLES = [0.4, 0.6, 999.]

"""
#@ HIGH-ENERGY SELECTION FOR RESIDUAL ANALYSIS (SKIMMING)
BGO_REC_MIN               = 10000  
BGO_REC_MAX               = None
#BGO_STK_DISTANCE_MATCH   = None
BGO_STK_ANGULAR_MATCH     = 0.2 
BGO_STK_COG_MATCH         = 50.
CLUSTER_ENERGY_CUT        = None
MAX_TRACKS_IN_EVENT       = 1
STK_MAX_CLUSTERS          = 999999  # 1200
TRACK_CLUSTER_ENERGY_MAX  = 99999.  # 160.
TRACK_CHISQCUT            = 99999.  # 48.
"""
"""
#@ HIGH-ENERGY SELECTION FOR RESIDUAL ANALYSIS (PLOTS)
BGO_REC_MIN               = 100000   
BGO_REC_MAX               = None
#BGO_STK_DISTANCE_MATCH   = None
BGO_STK_ANGULAR_MATCH     = 0.2 
BGO_STK_COG_MATCH         = 50.
CLUSTER_ENERGY_CUT        = None
MAX_TRACKS_IN_EVENT       = 1
STK_MAX_CLUSTERS          = 1200
TRACK_CLUSTER_ENERGY_MAX  = 160.
TRACK_CHISQCUT            = 48.
"""

#@ NO CUTS
BGO_REC_MIN               = None
BGO_REC_MAX               = None
#BGO_STK_DISTANCE_MATCH   = None
BGO_STK_ANGULAR_MATCH     = None
BGO_STK_COG_MATCH         = None
CLUSTER_ENERGY_CUT        = None
MAX_TRACKS_IN_EVENT       = None
STK_MAX_CLUSTERS          = None
TRACK_CLUSTER_ENERGY_MAX  = None
TRACK_CHISQCUT            = None

class Alignment:
    
    logs = {
            "ERROR" : 0,
            "INFO"  : 1,
            "DEBUG" : 2
            }

    #SENSOR_MARGINS = [0.0 , 95.0, 190.5,  286.0 , 385.0]
    SENSOR_MARGINS = [0.0 , 95.0, 190.5,  286.0 , 400.0]
    
    RESIDUAL_HISTO_X_NAME = "residuals_x"
    RESIDUAL_HISTO_Y_NAME = "residuals_y"
    RESIDUAL_HISTO_NBINS  = 2000
    RESIDUAL_HISTO_MIN    = - 10
    RESIDUAL_HISTO_MAX    =   10
    
    CHISQ_HISTO_X_NAME = "chisqx"
    CHISQ_HISTO_Y_NAME = "chisqy"
    CHISQ_HISTO_NBINS  = 20000
    CHISQ_HISTO_MIN    = 0
    CHISQ_HISTO_MAX    = 100  
    
    COV_HISTO_XX_NAME = "covxx"
    COV_HISTO_YY_NAME = "covyy"
    COV_HISTO_XY_NAME = "covxy"
    COV_HISTO_NBINS  = 20000
    COV_HISTO_MIN    = 0
    COV_HISTO_MAX    = 0.1

    VERTEX_2D_X_NBINS = 1200
    VERTEX_2D_X_MIN = -600
    VERTEX_2D_X_MAX =  600
    VERTEX_2D_Z_NBINS = 3000
    VERTEX_2D_Z_MIN = -500
    VERTEX_2D_Z_MAX =  1000
    FOCUS_DZ_MIN = 12
    FOCUS_DZ_MAX = 14
    VERTEX_DZ_NBINS = 800000
    VERTEX_DZ_MIN   = -400
    VERTEX_DZ_MAX   =  400

    HISTO_INCLINE_MIN   = -10
    HISTO_INCLINE_MAX   =  10
    HISTO_INCLINE_NBINS =  10
    HISTO_NSTRIPS_MIN   =  0
    HISTO_NSTRIPS_MAX   =  100
    HISTO_NSTRIPS_NBINS =  100
    HISTO_CHARGE_MIN   =  0
    HISTO_CHARGE_MAX   =  1000
    HISTO_CHARGE_NBINS =  2000
     
    
    RESIDUAL_2D_HISTO_RES_NBINS  =  400
    RESIDUAL_2D_HISTO_RES_MIN    =  - 2
    RESIDUAL_2D_HISTO_RES_MAX    =    2
    RESIDUAL_2D_HISTO_COOR_NBINS =  380
    RESIDUAL_2D_HISTO_COOR_MIN   = -380
    RESIDUAL_2D_HISTO_COOR_MAX   =  380
    RESIDUAL_2D_HISTO_COORSECTION_NSTEPS = 95 
    
    RESIDUAL_2D_HISTO_RESX_COORX_NAME = "resx_x"
    RESIDUAL_2D_HISTO_RESX_POSV_COORX_NAME = "resx_x_posvy"
    RESIDUAL_2D_HISTO_RESX_NEGV_COORX_NAME = "resx_x_negvy"
    RESIDUAL_2D_HISTO_RESX_COORY_NAME = "resx_y"
    RESIDUAL_2D_HISTO_RESY_COORX_NAME = "resy_x"
    RESIDUAL_2D_HISTO_RESY_COORY_NAME = "resy_y"
    RESIDUAL_2D_HISTO_RESY_POSV_COORY_NAME = "resy_y_posvx"
    RESIDUAL_2D_HISTO_RESY_NEGV_COORY_NAME = "resy_y_negvx"
    
    VERTICAL_TRACKS_ARG    = "--verttracks"
    BAD_CHANNELS_ARG       = "--bdchnl"
    NAME_VERSION           = "--ver"
    SELECTION_NO_CRACKS    = "--nocracks"
    ALIGNMENT_FILE         = "--align"
    SKIP_EVENTS_ARG        = "--skipevents"
    MAX_EVENTS_ARG         = "--maxevents"
    NTRACKS_FILE           = "--ntracksfile"
    TRACKCHISQMAX_ARG      = "--trackchisqmax"      
    TRACKCHISQMAX_AL_ARG   = "--aligntrackchisqmax"      
    ALIGNSENSORS_ARG       = "--alignsensors" 
    UNBIASED_ARG           = "--unbiased"
    FORBID_MORETHANONE_ARG = "--notmorethanone"
    ISRUNNINGONALIGNED_ARG = "--runningonaligned"
    CUT_UNBIASED_RESID     = "--cutunbiasres"
    CUT_OTHERPROJECTION    = "--cutotherprojection"
    DUMP_DATA_ARG          = "--dumpdata"
    NO_EVT_SELECTION_ARG   = "--noevtselection" 
    RUN_ON_COMPR_DATA_ARG  = "--runoncomprdata"
    NO_INSPECT_ARG         = "--noinspect"
    NO_ALIGNMENT_ARG       = "--noalignment"
    AVAILABLE_OPTIONS      = [
                           VERTICAL_TRACKS_ARG,
                           BAD_CHANNELS_ARG,
                           NAME_VERSION,
                           SELECTION_NO_CRACKS,
                           ALIGNMENT_FILE,
                           SKIP_EVENTS_ARG,
                           MAX_EVENTS_ARG,
                           NTRACKS_FILE,
                           TRACKCHISQMAX_ARG,
                           TRACKCHISQMAX_AL_ARG,
                           ALIGNSENSORS_ARG,
                           UNBIASED_ARG,
                           CUT_UNBIASED_RESID,
                           FORBID_MORETHANONE_ARG,
                           ISRUNNINGONALIGNED_ARG,
                           DUMP_DATA_ARG,
                           NO_EVT_SELECTION_ARG,
                           RUN_ON_COMPR_DATA_ARG,
                           NO_INSPECT_ARG,
                           NO_ALIGNMENT_ARG
                           ]
      
    
    CHISQ_NAME       = "chisq/ndof"
    CHISQ_NBINS      = 10000
    CHISQ_MIN        = 0
    CHISQ_MAX        = 100
    
    
    HISTO_NXCLUSTER_VS_Y_NAME  = "NXCLUSTER_VS_Y"
    HISTO_NYCLUSTER_VS_X_NAME  = "NYCLUSTER_VS_X"
    HISTO_NXCLUSTER_VS_Y_MIN   = -380
    HISTO_NXCLUSTER_VS_Y_MAX   =  380
    HISTO_NXCLUSTER_VS_Y_NBINS =  76000
    
    
    HITO_CENTRAL_CRACK_X = "central_crack_x"
    HITO_CENTRAL_CRACK_Y = "central_crack_y"
    HITO_CENTRAL_CRACK_MIN   = -10
    HITO_CENTRAL_CRACK_MAX   =  10
    HITO_CENTRAL_CRACK_NBINS = 10000
    
    
    MIN_GOOG_XY_HITS            = 6
    TRACK_TAN_NOCUT_MIN         = -1.
    TRACK_TAN_NOCUT_MAX         = 9999.
    TRACK_TAN_VERTCUT_MIN       = -1.   #0.6    #-1.
    TRACK_TAN_VERTCUT_MAX       = 0.4   # 9999. #0.2
    TRACK_CRACK_REGION_NSTRIPS  = 5 # RO stips ==> ~1mm
    TRACK_CRACK_REGION_DISTANCE = 3 # mm
    
    N_RO_STRIPS = 384
    N_STK_LAYERS =  6
   
    
    #MATRIX_M_ALYERS_RANK = 6 # 3 translations + 3 rotation angles  
    
    OUT_DERIVATIVES_FILE = "derivatives.dump"
    OUT_DATA_FILE        = "data.npy"
    ROOTFILE_PREFIX      = "out_ALIGN"
    TMP_ROOT_FILE_NAME   = "tmp.root"
    

    def __get_silicon_sensor_index__(self, coordinate):  # 0 - 3
        for i in xrange(len(self.SENSOR_MARGINS)-1):
            if abs(coordinate) < self.SENSOR_MARGINS[i]: continue
            if abs(coordinate) > self.SENSOR_MARGINS[i+1]: continue
            return i
        self.__print_error__("__get_silicon_sensor_index__: unrecognized coordinate: %f ==> throwing Exception!"%coordinate)
        raise Exception
    
    
    def get_linear_fit(self,x, y, getderivatives = False):
        #@
        #@  y = ax + b
        #@
        
        #self.__measuretime_start__(self.get_linear_fit.__name__)
        
        n   = 0
        Ex  = 0
        Ex2 = 0
        Ey  = 0
        Exy = 0
        for i in xrange(len(x)):
            argument = x[i]
            value    = y[i]
            Ex+= argument
            Ex2+=argument**2
            Ey+=value
            Exy+=argument*value
            n+=1
        
        
        """ 
        a = (n*Exy - Ex* Ey )/(n*Ex2 - Ex**2)
        b = (Ex * Exy - Ey * Ex2) / (Ex**2 - n * Ex2)
        """
        factor1 = n*Exy - Ex* Ey
        factor2 = n*Ex2 - Ex**2
        a = factor1 / factor2
        
        factor3 = Ex * Exy - Ey * Ex2
        factor4 = Ex**2 - n * Ex2
        b = factor3 / factor4
        
        if not getderivatives:
            return a, b
        
        
        # evaluate derivatives
        #da_darg = []
        #da_dval = []
        #db_darg = []
        #db_dval = []
        da_darg = np.array([0.] * len(x))
        da_dval = np.array([0.] * len(x))
        db_darg = np.array([0.] * len(x))
        db_dval = np.array([0.] * len(x))
        for i  in xrange(len(x)):
            #
            part1 = np.multiply(n,y[i]) - Ey
            part2 = np.multiply(n,x[i]) - Ex
            dadarg = part1 - np.divide(
                                       np.multiply(factor1, np.multiply(2,part2)), 
                                       factor2)
            dadarg = np.divide(dadarg, factor2)
            #da_darg.append(dadarg)
            da_darg[i] = dadarg
            #
            dadval = np.divide( part2 , factor2)
            #da_dval.append(dadval)
            da_dval[i] = dadval
            #
            part3 = Exy + np.multiply(y[i], Ex) - np.multiply(2, np.multiply(Ey, x[i])) 
            part4 = np.multiply(2, Ex - np.multiply(n,x[i]))
            dbdarg = part3 - np.divide(
                                       np.multiply(factor3, part4),
                                       factor4)
            dbdarg = np.divide(dbdarg, factor4)
            #db_darg.append(dbdarg)
            db_darg[i] = dbdarg
            #
            dbdval = np.divide(
                               np.multiply(x[i],Ex) - Ex2,
                               factor4
                               )
            #db_dval.append(dbdval)
            db_dval[i] = dbdval
            
        # get residuals
        residuals = []
        for i in xrange(len(x)):
            arg = x[i]
            val = y[i]
            projval = np.multiply(a,arg) + b
            residuals.append(val - projval)
        residuals = np.array(residuals)
         
            
        # chi-sq
        chisq = np.sum([np.multiply(res,res) for res in residuals])
        
        
        #self.__measuretime_stop__(self.get_linear_fit.__name__)
        
        return a, b, da_darg, da_dval, db_darg, db_dval, residuals, chisq
    
    # ancillary derivatives for the linear fit       





    def dx_dpar(self,x,y,indx):
        #@
        #@ Parameters:
        #@    deltax
        #@    deltay
        #@    deltaz
        #@    dthetax
        #@    dthetay
        #@    dthetaz
        #@    alphax
        #@    alphay
        #@    betax
        #@    betay 
        #return [1.,  0.,  0.,  0.,  0.,  -y , 0.,    0. ,  0.,     0.  ]
        if indx==0:
            return 1.
        elif indx==5:   
            return -y
        return 0.

             
    def dy_dpar(self,x,y,indx):
        #return [0.,  1.,  0.,  0.,  0.,   x , 0.,    0. ,  0.,     0.  ]
        if indx==1:
            return 1.
        elif indx==5:   
            return x
        return 0.
    
    def dz_dpar(self,x,y,indx):
        #return [0.,  0.,  1.,  y, -x,    0. , x**2 , y**2, x**3 , y**3 ]
        if          indx==2:
            return 1.
        elif        indx==3:   
            return y
        elif        indx==4:   
            return -x
        elif        indx==6:   
            return x**2
        elif        indx==7:   
            return y**2
        elif        indx==8:   
            return x**3
        elif        indx==9:   
            return y**3
        return 0.
    
    def dval_dpar(self,clusterladders, x, y, ladder, parameter_index, dval_dpar_function):
        #@ val:  arguments of the linear fit
        #@ par:  alignment parameter
        #@ parameter_index: 0 - 5   (3 translations and 3 rotations)

        #self.__measuretime_start__(self.dval_dpar.__name__)
        result = np.array([0.]*len(clusterladders))
        for i in xrange(len(clusterladders)):
            if clusterladders[i] != ladder:
                continue
            #result[i] = dval_dpar_function(x[i], y[i]) [parameter_index]
            result[i] = dval_dpar_function(x[i], y[i], parameter_index)
        #self.__measuretime_stop__(self.dval_dpar.__name__)
        return result
    


    """
    def dval_dpar(self,clusterladders, x, y, ladder, parameter_index, dval_dpar_function):
        
        #@ val:  arguments of the linear fit
        #@ par:  alignment parameter
        #@ parameter_index: 0 - 5   (3 translations and 3 rotations)
        
        #self.__measuretime_start__(self.dval_dpar.__name__)

        #result = np.array([0.]*len(clusterladders))
        result = []
        for x_val, y_val, lad_val in izip(x,y,clusterladders):
            if lad_val!= ladder:
                result.append(0.)
            else:
                #result.append(dval_dpar_function(x_val, y_val)[parameter_index])
                result.append(dval_dpar_function(x_val, y_val,parameter_index))
        result =  np.array(result)
        #self.__measuretime_stop__(self.dval_dpar.__name__)
        return result
    """ 
      
    def dval_dpar_xtrack(self,clusterladders, x, y, ladder, parameter_index):
        return self.dval_dpar(clusterladders, x, y, ladder, parameter_index, self.dx_dpar)
       
    def dval_dpar_ytrack(self,clusterladders, x, y, ladder, parameter_index):
        return self.dval_dpar(clusterladders, x, y, ladder, parameter_index, self.dy_dpar)
        
    def darg_dpar(self,clusterladders, x, y, ladder, parameter_index):
        #@ val:  arguments of the linear fit
        #@ par:  alignment parameter
        #@ parameter_index: 0 - 5   (3 translations and 3 rotations)
        #self.__measuretime_start__(self.darg_dpar.__name__)
        #result = []
        result = np.array([0.]*len(clusterladders))
        for i in xrange(len(clusterladders)):
            if clusterladders[i] != ladder:
                #result.append(0.)
                continue
            
            
            #result.append(
            #              self.dz_dpar(x[i], y[i])[parameter_index]
            #              )
            #result[i] = self.dz_dpar(x[i], y[i])[parameter_index]
            result[i] = self.dz_dpar(x[i], y[i],parameter_index)
        #self.__measuretime_stop__(self.darg_dpar.__name__) 
        return result
    
    def da_dpar(self,da_dval, dval_dpar, da_darg, darg_dpar):
        #@
        #@ a is a parameter of linear fit y = ax + b
        #@
        #self.__measuretime_start__(self.da_dpar.__name__)
        ##result =  np.sum(np.multiply(da_dval, dval_dpar)) + np.sum(np.multiply(da_darg, darg_dpar))
        ##result =  np.sum(np.array(da_dval) *  np.array(dval_dpar)) + np.sum(np.array(da_darg) *np.array( darg_dpar))
        #result =  np.sum(da_dval *  dval_dpar) + np.sum(da_darg * darg_dpar)
        #self.__measuretime_stop__(self.da_dpar.__name__)
        #return result

        #return  np.sum(da_dval *  dval_dpar) + np.sum(da_darg * darg_dpar)
        return  np.sum(np.multiply(da_dval ,dval_dpar)) + np.sum(np.multiply(da_darg, darg_dpar))
    
    def db_dpar(self,db_dval, dval_dpar, db_darg, darg_dpar):
        #@ b is a parameter of linear fit y = ax + b
        return self.da_dpar(db_dval, dval_dpar, db_darg, darg_dpar)
    
        
    
    
    
    
    def Init(self):
        #
        # parse command line arguments
        #
        self.log = self.logs["DEBUG"]
        
        self.badchannels = {}
        self.read_bad_channels()
        
        self.track_tan_cut_min = self.TRACK_TAN_NOCUT_MIN
        self.track_tan_cut_max = self.TRACK_TAN_NOCUT_MAX
        self.do_straight_tracks()

        self.N_STK_LADDERS = 192
        self.__all_sensors__ = False
        self.do_all_sensors()
        
        self.versiontag = ""
        self.di_version_tag()
        
        self.exclude_cracks = False
        self.do_exclude_cracks()
        
        self.alignments = {}
        self.get_initial_alignment()
        
        self.skipevents = 0
        self.get_skip_events()
        
        self.get_max_events()

        self.__unbiased__ = False
        self.get_unbiased()
        
        self.__is_ntracks_file__ = False
        self.normalize_n_tracks = [1.] * self.N_STK_LADDERS
        self.get_n_tracks_file()
            
        self.chi_sq_cut = 999999999.
        self.get_track_chi_sq_cut()


        self.aligntrack_chi_sq_cut = 9999999.
        self.get_aligntrack_chi_sq_cut()

        self.__unbiased_track_residuals_max__ = None
        self.get_cut_unbiased_track_residuals()

        self.__unbiased_track_residuals_otherprojection_max__ = None
        self.get_cut_other_projection()

        self.__forbid_more_than_one__ = False
        self.get_forbid_more_than_one()

        self.__running_on_aligned__ = False
        self.get_is_running_on_aligned()

        self.__dump_data__ = False
        self.get_is_dump_data()

        self.__run_on_compr_data__ = False
        self.get_run_on_copr_data()

        self.__no_evt_selection__ = False
        self.get_no_evt_selection()

        self.__no_inspect__ = False
        self.get_no_inspect()

        self.__no_alignment__ = False
        self.get_no_alignment()
        
        
        args = [arg for arg in sys.argv if arg.strip()[0]=="-"]
        if args:
            self.__print_error__("Unrecognized parameter: " + args[0]+" ==> throwing exception!")
            raise Exception
        
        
        #self.inname = sys.argv[1]
        #self.fin = TFile(self.inname, "READ")
        #self.t = self.fin.Get("CollectionTree")
        
        if self.__run_on_compr_data__:
            fname = sys.argv[1]
            self.loaded_data = np.load(fname)
            nfiles = 1
        else:
            self.t = TChain("CollectionTree")
            nfiles, fname = self.__add_files_to_chain__()  
            self.stkclusters  = TClonesArray("DmpStkSiCluster")
            #self.stkadccounts = TClonesArray("DmpStkLadderAdc")
            self.stktracks    = TClonesArray("DmpStkTrack")
            self.bgorec       = DmpEvtBgoRec()
            self.bgohits      = DmpEvtBgoHits()
            self.t.SetBranchAddress("StkKalmanTracks",     self.stktracks)
            self.t.SetBranchAddress("StkClusterCollection",  self.stkclusters)
            #self.t.SetBranchAddress("DmpStkLadderAdcCollection", self.stkadccounts)
            self.t.SetBranchAddress("DmpEvtBgoRec",self.bgorec)
            self.t.SetBranchAddress("DmpEvtBgoHits", self.bgohits)
            # branches
            self.b_stktracks   = self.t.GetBranch("StkKalmanTracks")
            self.b_bgohits     = self.t.GetBranch("DmpEvtBgoHits")
            self.b_stkclusters = self.t.GetBranch("StkClusterCollection")
            self.b_bgorec      = self.t.GetBranch("DmpEvtBgoRec")
        
        #outputfilename
        if nfiles == 1:
            self.OUT_DERIVATIVES_FILE = self.OUT_DERIVATIVES_FILE.split(".")[0] + "_" + fname.split("/")[-1].split(".")[0] + "." +  self.OUT_DERIVATIVES_FILE.split(".")[1]
            self.OUT_DATA_FILE        = self.OUT_DATA_FILE.split(".")[0]        + "_" + fname.split("/")[-1].split(".")[0] + "." +  self.OUT_DATA_FILE .split(".")[1]
            self.ROOTFILE_PREFIX      = self.ROOTFILE_PREFIX + "_" + fname.split("/")[-1].split(".")[0] + "_"  
        if OUTPATH_ENV_NAME in os.environ:
            self.OUT_DERIVATIVES_FILE = os.environ[OUTPATH_ENV_NAME] + "/" +  self.OUT_DERIVATIVES_FILE
            self.ROOTFILE_PREFIX      = os.environ[OUTPATH_ENV_NAME] + "/" +  self.ROOTFILE_PREFIX  
            self.OUT_DATA_FILE        = os.environ[OUTPATH_ENV_NAME] + "/" +  self.OUT_DATA_FILE
         


        print "\n\n\n\n\n\n\n"
        print "test1 ", self.OUT_DERIVATIVES_FILE
        print "test2 ", self.ROOTFILE_PREFIX
        print "\n\n\n\n\n\n\n"
   
        if self.versiontag:
            self.OUT_DERIVATIVES_FILE = self.OUT_DERIVATIVES_FILE.split(".")[0] + self.versiontag + "." + self.OUT_DERIVATIVES_FILE.split(".")[-1]
            self.OUT_DATA_FILE        = self.OUT_DATA_FILE.split(".")[0]        + self.versiontag + "." + self.OUT_DATA_FILE.split(".")[-1]
             
        
        # derivatives
        self.dchish_dpars = { 
                             "derivatives": [[0. for j in xrange(N_DOF_LADDER)] for i in xrange(self.N_STK_LADDERS)],
                             "ntracks"    : 0,
                             "chisq"      : 0.0
                             }
        self.tracks_for_lad = [[0 for angle1 in ANGLES_FORALGORITHM] for i in xrange(self.N_STK_LADDERS)]
        
        
        # histograms
        self.histo_track_unnormalized_chisq_for_derivatives = TH1F("track_unnormalized_chisq", "track_unnormalized_chisq", 20000, 0, 200)
        self.histo_alignedtrack_chisq = TH1F("alignedtrack_chisq", "alignedtrack_chisq", 20000, 0, 200)
        self.histo_track_kalman_chisq_for_derivatives = TH1F("track_kalman_chisq", "track_kalman_chisq", 20000, 0, 200)
        self.histo_track_angle_x = TH1F("track_tanx", "track_tanx", 1000, 0, 10)
        self.histo_track_angle_y = TH1F("track_tany", "track_tany", 1000, 0, 10)
        self.histo_track_angle   = TH1F("track_tan",  "track_tany", 1000, 0, 10)
        self.histo_track_anglex_chisqx = TH2F("anglex_chisqx","anglex_chisqx", 1000, 0, 10, 1000, 0, 10)
        self.histo_track_angley_chisqy = TH2F("angley_chisqy","angley_chisqy", 1000, 0, 10, 1000, 0, 10)
        self.histo_track_angle_chisq   = TH2F("angle_chisq",  "angle_chisq",   1000, 0, 10, 1000, 0, 10)
        self.histo_nclustersbefore = TH1F("histo_nclustersbefore","histo_nclustersbefore", 10000, 0, 10000)
        self.histo_nclustersafter = TH1F("histo_nclustersafter","histo_nclustersafter", 10000, 0, 10000)
        self.histo_distance    = TH1F("histo_distance", "histo_distance",       10000, 0 , 100)
        self.histo_angdistance = TH1F("histo_angdistance", "histo_angdistance", 10000, 0 , 2)
        self.histo_chargex     = TH1F("chargex", "chargex", 100000,  0, 10000.)
        self.histo_chargey     = TH1F("chargey", "chargey", 100000,  0, 10000.)
        
        if not self.__run_on_compr_data__:
            self.histo_incline_nstrips_2D_x = TH2F("incline_nstrips_2D_x",  "incline_nstrips_2D_x",   	
				self.HISTO_INCLINE_NBINS, 
				self.HISTO_INCLINE_MIN, 
				self.HISTO_INCLINE_MAX, 
				self.HISTO_NSTRIPS_NBINS,
				self.HISTO_NSTRIPS_MIN, 
				self.HISTO_NSTRIPS_MAX)
            self.histo_incline_nstrips_2D_y = TH2F("incline_nstrips_2D_y",  "incline_nstrips_2D_y",   
				self.HISTO_INCLINE_NBINS, 
				self.HISTO_INCLINE_MIN, 
				self.HISTO_INCLINE_MAX, 
				self.HISTO_NSTRIPS_NBINS,
				self.HISTO_NSTRIPS_MIN,
				self.HISTO_NSTRIPS_MAX)
            self.histo_incline_nstrips_2D = TH2F("incline_nstrips_2D",  "incline_nstrips_2D",   
				self.HISTO_INCLINE_NBINS, 
				self.HISTO_INCLINE_MIN, 
				self.HISTO_INCLINE_MAX, 
				self.HISTO_NSTRIPS_NBINS,
				self.HISTO_NSTRIPS_MIN,
				self.HISTO_NSTRIPS_MAX)
            self.histo_incline_charge_2D_x = TH2F("incline_charge_2D_x",  "incline_charge_2D_x",      
				self.HISTO_INCLINE_NBINS, 
				self.HISTO_INCLINE_MIN, 
				self.HISTO_INCLINE_MAX, 
				self.HISTO_CHARGE_NBINS,
				self.HISTO_CHARGE_MIN,
				self.HISTO_CHARGE_MAX)
            self.histo_incline_charge_2D_y = TH2F("incline_charge_2D_y",  "incline_charge_2D_y",      
				self.HISTO_INCLINE_NBINS, 
				self.HISTO_INCLINE_MIN, 
				self.HISTO_INCLINE_MAX, 
				self.HISTO_CHARGE_NBINS, 
				self.HISTO_CHARGE_MIN,
				self.HISTO_CHARGE_MAX)
            self.histo_incline_charge_2D = TH2F("incline_charge_2D",  "incline_charge_2D",      
				self.HISTO_INCLINE_NBINS, 
				self.HISTO_INCLINE_MIN, 
				self.HISTO_INCLINE_MAX, 
				self.HISTO_CHARGE_NBINS,
				self.HISTO_CHARGE_MIN,
				self.HISTO_CHARGE_MAX)
            self.histo_incline_chargevacorrected_2D_x = TH2F("incline_chargevacorrected_2D_x",  "incline_chargevacorrected_2D_x",      
				self.HISTO_INCLINE_NBINS, 
				self.HISTO_INCLINE_MIN, 
				self.HISTO_INCLINE_MAX, 
				self.HISTO_CHARGE_NBINS,
				self.HISTO_CHARGE_MIN,
				self.HISTO_CHARGE_MAX)
            self.histo_incline_chargevacorrected_2D_y = TH2F("incline_chargevacorrected_2D_y",  "incline_chargevacorrected_2D_y",      
				self.HISTO_INCLINE_NBINS, 
				self.HISTO_INCLINE_MIN, 
				self.HISTO_INCLINE_MAX, 
				self.HISTO_CHARGE_NBINS, 
				self.HISTO_CHARGE_MIN, 
				self.HISTO_CHARGE_MAX)
            self.histo_incline_chargevacorrected_2D = TH2F("incline_chargevacorrected_2D",  "incline_chargevacorrected_2D",      
				self.HISTO_INCLINE_NBINS, 
				self.HISTO_INCLINE_MIN, 
				self.HISTO_INCLINE_MAX, 
				self.HISTO_CHARGE_NBINS, 
				self.HISTO_CHARGE_MIN, 
				self.HISTO_CHARGE_MAX)
            self.histo_incline_nstrips_x = []
            self.histo_incline_nstrips_y = []
            self.histo_incline_nstrips   = []
            self.histo_incline_charge_x  = [] 
            self.histo_incline_charge_y  = [] 
            self.histo_incline_charge    = [] 
            self.histo_incline_chargevacorr_x = []
            self.histo_incline_chargevacorr_y = []
            self.histo_incline_chargevacorr   = []
            for angle in ANGLES:
                self.histo_incline_nstrips_x.append(TH1F("nstrips_x_incl_%.1f"%angle, "nstrips_x_incl_%.1f"%angle, self.HISTO_NSTRIPS_NBINS, self.HISTO_NSTRIPS_MIN, self.HISTO_NSTRIPS_MAX ))
                self.histo_incline_nstrips_y.append(TH1F("nstrips_y_incl_%.1f"%angle, "nstrips_y_incl_%.1f"%angle, self.HISTO_NSTRIPS_NBINS, self.HISTO_NSTRIPS_MIN, self.HISTO_NSTRIPS_MAX ))
                self.histo_incline_nstrips.append(TH1F("nstrips_incl_%.1f"%angle, "nstrips_incl_%.1f"%angle, self.HISTO_NSTRIPS_NBINS, self.HISTO_NSTRIPS_MIN, self.HISTO_NSTRIPS_MAX ))
                self.histo_incline_charge_x.append(TH1F("charge_x_incl_%.1f"%angle, "charge_x_incl_%.1f"%angle, self.HISTO_CHARGE_NBINS, self.HISTO_CHARGE_MIN, self.HISTO_CHARGE_MAX))
                self.histo_incline_charge_y.append(TH1F("charge_y_incl_%.1f"%angle, "charge_y_incl_%.1f"%angle, self.HISTO_CHARGE_NBINS, self.HISTO_CHARGE_MIN, self.HISTO_CHARGE_MAX))
                self.histo_incline_charge.append(TH1F("charge_incl_%.1f"%angle, "charge_incl_%.1f"%angle, self.HISTO_CHARGE_NBINS, self.HISTO_CHARGE_MIN, self.HISTO_CHARGE_MAX))
                self.histo_incline_chargevacorr_x.append(TH1F("chargevacorr_x_incl_%.1f"%angle, "chargevacorr_x_incl_%.1f"%angle, self.HISTO_CHARGE_NBINS, self.HISTO_CHARGE_MIN, self.HISTO_CHARGE_MAX)) 
                self.histo_incline_chargevacorr_y.append(TH1F("chargevacorr_y_incl_%.1f"%angle, "chargevacorr_y_incl_%.1f"%angle, self.HISTO_CHARGE_NBINS, self.HISTO_CHARGE_MIN, self.HISTO_CHARGE_MAX))
                self.histo_incline_chargevacorr.append(TH1F("chargevacorr_incl_%.1f"%angle, "chargevacorr_incl_%.1f"%angle, self.HISTO_CHARGE_NBINS, self.HISTO_CHARGE_MIN, self.HISTO_CHARGE_MAX)) 
            
        
        # output ROOT file
        #self.fout = TFile(self.ROOTFILE_PREFIX + self.versiontag +".root", "RECREATE")
        self.fout = TFile(self.TMP_ROOT_FILE_NAME, "RECREATE")

        # output data file
        if self.__dump_data__:
            self.data_dump_all = {"x":[], "z_x":[], "y_x":[], "ladx":[], "y":[], "z_y":[], "x_y":[], "lady":[], "covxx":[], "covyy":[], "vals":[]}

        # some common loops
        self.ALL_LADS     = np.array(xrange(self.N_STK_LADDERS))
        self.ALL_PARS_DOF = np.array(xrange(N_DOF_LADDER))
    
    """    
    def __create_v_unity_matrix__(self):
        self.V = TMatrixD(2,2)
        self.V[0][0] = 1
        self.V[1][1] = 1
        self.V[0][1] = 0
        self.V[1][0] = 0
    """
        
    """
    def __get_v_matrix(self, track, layerfromcalorimeter):
        return self.V
    """
    
    def __measuretime_start__(self , name):
        if name not in self.__time_elapsed__.keys():
            self.__time_elapsed__[name] = 0.
            self.__number_of_calls__[name] = 0.
        #self.__time_elapsed_start__[name] = datetime.datetime.now()
        self.__time_elapsed_start__[name] = time.clock()
    
    def __measuretime_stop__(self , name):
        if name not in self.__time_elapsed_start__:
            print "[__measuretime_stop__] Function time measurment not started!"
            raise Exception
        
        if self.__time_elapsed_start__[name] is None:
            print "[__measuretime_stop__] Function time measurment was stopped but not started again!"
            raise Exception
        
        #elapsed = datetime.datetime.now() - self.__time_elapsed_start__[name]
        #elapsed = elapsed.total_seconds() *1000. + elapsed.total_microseconds()*0.001
        elapsed = time.clock() - self.__time_elapsed_start__[name]
        self.__time_elapsed__[name]+= elapsed
        self.__number_of_calls__[name]+= 1
        self.__time_elapsed_start__[name] = None
        #self.__time_elapsed_function__ = None
        
        
    
    def __init__(self):
        self.histo_vertex_dz = TH1F("vertex_dz", "vertex_dz", self.VERTEX_DZ_NBINS,  self.VERTEX_DZ_MIN, self.VERTEX_DZ_MAX)
        self.histo_vertex_xz = TH2F("vertex_xz", "vertex_xz", self.VERTEX_2D_X_NBINS,  self.VERTEX_2D_X_MIN, self.VERTEX_2D_X_MAX, self.VERTEX_2D_Z_NBINS, self.VERTEX_2D_Z_MIN, self.VERTEX_2D_Z_MAX)
        self.histo_vertex_yz = TH2F("vertex_yz", "vertex_yz", self.VERTEX_2D_X_NBINS,  self.VERTEX_2D_X_MIN, self.VERTEX_2D_X_MAX, self.VERTEX_2D_Z_NBINS, self.VERTEX_2D_Z_MIN, self.VERTEX_2D_Z_MAX)
        self.histo_vertex_focus_dz = TH1F("vertex_focus_dz", "vertex_focus_dz", self.VERTEX_DZ_NBINS,  self.VERTEX_DZ_MIN, self.VERTEX_DZ_MAX)
        self.histo_vertex_focus_xz = TH2F("vertex_focus_xz", "vertex_focus_xz", self.VERTEX_2D_X_NBINS,  self.VERTEX_2D_X_MIN, self.VERTEX_2D_X_MAX, self.VERTEX_2D_Z_NBINS, self.VERTEX_2D_Z_MIN, self.VERTEX_2D_Z_MAX)
        self.histo_vertex_focus_yz = TH2F("vertex_focus_yz", "vertex_focus_yz", self.VERTEX_2D_X_NBINS,  self.VERTEX_2D_X_MIN, self.VERTEX_2D_X_MAX, self.VERTEX_2D_Z_NBINS, self.VERTEX_2D_Z_MIN, self.VERTEX_2D_Z_MAX)
	self.ngoodtracks = TH1F("ngoodtracks","ngoodtracks",100,0,100)
        self.residx = [{} for i in xrange(len(ANGLES) + 1)]
        self.residy = [{} for i in xrange(len(ANGLES) + 1)]
        self.chisqx = [{} for i in xrange(len(ANGLES) + 1)]
        self.chisqy = [{} for i in xrange(len(ANGLES) + 1)]
        self.covxx  = [{} for i in xrange(len(ANGLES) + 1)]
        self.covyy  = [{} for i in xrange(len(ANGLES) + 1)]
        self.covxy  = [{} for i in xrange(len(ANGLES) + 1)]
        self.residx_vs_x = [{} for i in xrange(len(ANGLES) + 1)]
        self.residx_vs_x_posv = [{} for i in xrange(len(ANGLES) + 1)]
        self.residx_vs_x_negv = [{} for i in xrange(len(ANGLES) + 1)]
        self.residx_vs_y = [{} for i in xrange(len(ANGLES) + 1)]
        self.residy_vs_y = [{} for i in xrange(len(ANGLES) + 1)]
        self.residy_vs_y_posv = [{} for i in xrange(len(ANGLES) + 1)]
        self.residy_vs_y_negv = [{} for i in xrange(len(ANGLES) + 1)]
        self.residy_vs_x = [{} for i in xrange(len(ANGLES) + 1)]
        
	self.resid_all_lad_y = [{} for i in xrange(len(ANGLES) + 1)]
	self.resid_all_lad_x = [{} for i in xrange(len(ANGLES) + 1)]

        self.trackchisqx = TH1F(
                               self.CHISQ_NAME, 
                               self.CHISQ_NAME, 
                               self.CHISQ_NBINS, 
                               self.CHISQ_MIN, 
                               self.CHISQ_MAX 
                               )
        self.trackchisqy = TH1F(
                               self.CHISQ_NAME, 
                               self.CHISQ_NAME, 
                               self.CHISQ_NBINS, 
                               self.CHISQ_MIN, 
                               self.CHISQ_MAX 
                               )
        #self.matrix_M_layerx = {}
        #self.matrix_M_layery = {}
        #self.matrix_NU_layerx = {}
        #self.matrix_NU_layery = {}
        
        self.nxcluster_vs_y = [
                               TH1F(
                                    self.HISTO_NXCLUSTER_VS_Y_NAME +"_layer_"+str(i),
                                    self.HISTO_NXCLUSTER_VS_Y_NAME +"_layer_"+str(i),
                                    self.HISTO_NXCLUSTER_VS_Y_NBINS,
                                    self.HISTO_NXCLUSTER_VS_Y_MIN,
                                    self.HISTO_NXCLUSTER_VS_Y_MAX
                                    )
                               for i in xrange(self.N_STK_LAYERS)
                               ]
        self.nycluster_vs_x = [
                               TH1F(
                                    self.HISTO_NYCLUSTER_VS_X_NAME +"_layer_"+str(i),
                                    self.HISTO_NYCLUSTER_VS_X_NAME +"_layer_"+str(i),
                                    self.HISTO_NXCLUSTER_VS_Y_NBINS,
                                    self.HISTO_NXCLUSTER_VS_Y_MIN,
                                    self.HISTO_NXCLUSTER_VS_Y_MAX
                                    )
                               for i in xrange(self.N_STK_LAYERS)
                               ]
        
        self.centralcrack_x =  TH1F(
                                    self.HITO_CENTRAL_CRACK_X,
                                    self.HITO_CENTRAL_CRACK_X,
                                    self.HITO_CENTRAL_CRACK_NBINS,
                                    self.HITO_CENTRAL_CRACK_MIN,
                                    self.HITO_CENTRAL_CRACK_MAX
                                    )
                              
        self.centralcrack_y =  TH1F(
                                    self.HITO_CENTRAL_CRACK_Y,
                                    self.HITO_CENTRAL_CRACK_Y,
                                    self.HITO_CENTRAL_CRACK_NBINS,
                                    self.HITO_CENTRAL_CRACK_MIN,
                                    self.HITO_CENTRAL_CRACK_MAX
                                    )
        
        self.MAX_ENTRIES = None
        
        # time measurement
        self.__time_elapsed__ = {}
        self.__time_elapsed_start__ = {}
        self.__number_of_calls__ = {}
        
        
                               
        #self.fout = TFile(self.inname.split("/")[-1].split(".root")[0] + "_ALIGN.root", "RECREATE")
        
    def __add_files_to_chain__(self):
        nfiles = 0
        fname = None
        for f in sys.argv[1:]:
            if f[0]=="-": continue
            self.__print_info__("Using file: "+f)
            self.t.Add(f)
            nfiles+=1
            fname = f
        return nfiles, fname 
        
        
    def Run(self):
        self.__measuretime_start__("TOTAL")
        
        if self.__run_on_compr_data__:
            #
            #  Running on compressed data
            #
            nentries = len(self.loaded_data["x"])
            if self.MAX_ENTRIES is not None:
                nentries = min(nentries,self.skipevents + self.MAX_ENTRIES)

            entry     = -1
            processed =  0
            #  ----------- event loop COMPRESSED data -----------------
            for [self.data_x, 
                 self.data_z_x, 
                 self.data_y_x, 
                 self.data_ladx, 
                 self.data_covxx, 
                 self.data_y,
                 self.data_z_y, 
                 self.data_x_y, 
                 self.data_lady, 
                 self.data_covyy, 
                 self.data_val] in izip(self.loaded_data["x"],
                                        self.loaded_data["z_x"],
                                        self.loaded_data["y_x"],
                                        self.loaded_data["ladx"],
                                        self.loaded_data["covxx"],
                                        self.loaded_data["y"],
                                        self.loaded_data["z_y"],
                                        self.loaded_data["x_y"],
                                        self.loaded_data["lady"],
                                        self.loaded_data["covyy"],
                                        self.loaded_data["vals"]):
                entry+=1
                if entry < self.skipevents: continue
                sys.stdout.write("\r Processing file:    %10d         %10d / %10d ( %2d%% )"%(entry, entry - self.skipevents,nentries, 100 * (entry-self.skipevents)/nentries))
                sys.stdout.flush()

                self.process_event()
                
                processed+=1
                if processed == self.MAX_ENTRIES: break
            #  ............ end of event loop COMPRESSED data
        

        else:
            #
            # Running on root data
            #
            nentries = self.t.GetEntries()
            if self.MAX_ENTRIES is not None:
                nentries = min(nentries,self.skipevents + self.MAX_ENTRIES)
            
            #  ----------- event loop ROOT  data -----------------
            for i in xrange(self.skipevents, nentries):
            
                sys.stdout.write("\r Processing file:    %10d         %10d / %10d ( %2d%% )"%(i, i - self.skipevents,nentries, 100 * (i-self.skipevents)/nentries))
                sys.stdout.flush()
            
                #self.t.GetEntry(i)
                self.b_stktracks.GetEntry(i) 
                self.b_bgohits.GetEntry(i)
                self.b_stkclusters.GetEntry(i)
                self.b_bgorec.GetEntry(i)

                self.process_event()
            #  ............ end of event loop ROOT data
            
        self.__measuretime_stop__("TOTAL")
            
            
    def __print_debug__(self, message):
        if self.log >= self.logs["DEBUG"]:
            print "DEBUG : ", message
            
    def __print_info__(self, message):
        if self.log >= self.logs["INFO"]:
            print "INFO  : ", message
            
    def __print_error__(self, message):
        if self.log >= self.logs["ERROR"]:
            print "ERROR : ", message
            
    







    def process_event(self):
        # 1. Loop over tracks and find good ones
        #if not self.__no_evt_selection__:      
        if not self.__run_on_compr_data__:
            all_track = []
            all_slope = []
            all_xclusters = []
            all_covxx = []
            all_yclusters = []
            all_covyy = []
            all_a_x_nonalign = []
            all_b_x_nonalign = []
            all_a_y_nonalign = []
            all_b_y_nonalign = []
            all_trackchisq_unnormalized = []
            all_chisqx = []
            all_chisqy = []
            all_clusterx_y = []
            all_clustery_x = []
            all_chisqkalman = []
            all_maxclustercharge = []
            all_minclustercharge = []
            all_bgoangdist = []
            all_bgocogdist = [] 


            # bgo energy cuts
            energy = sum([en for en in self.bgohits.fEnergy])
            if BGO_REC_MAX is not None and energy> BGO_REC_MAX: return
            if BGO_REC_MIN is not None and energy< BGO_REC_MIN: return
            all_bgoenergy = energy

            # stk cluster mu;ltiplicity cut
            self.histo_nclustersbefore.Fill(self.stkclusters.GetLast()+1)
            if STK_MAX_CLUSTERS is not None and self.stkclusters.GetLast()+1 > STK_MAX_CLUSTERS : return
            self.histo_nclustersafter.Fill(self.stkclusters.GetLast()+1)
            all_nclusters = self.stkclusters.GetLast()+1


            # stk multiplicity cut
            if MAX_TRACKS_IN_EVENT is not None and self.stktracks.GetLast() + 1 > MAX_TRACKS_IN_EVENT: return
            all_ntracks = self.stktracks.GetLast() + 1

            # bgo direction
            bgo_x = self.bgorec.GetTrajectoryLocation2D().x()
            bgo_y = self.bgorec.GetTrajectoryLocation2D().y()
            bgo_z = self.bgorec.GetTrajectoryLocation2D().z()
            bgo_incl_x = self.bgorec.GetSlopeXZ()
            bgo_incl_y = self.bgorec.GetSlopeYZ()
            bgo_intercept_x = bgo_x - bgo_incl_x * bgo_z
            bgo_intercept_y = bgo_y - bgo_incl_y * bgo_z

            # bgo COG
            bgo_cog_x = 0.
            bgo_cog_y = 0.
            bgo_cog_zx = 0.
            bgo_cog_zy = 0.
            bgo_cog_ex = 0.
            bgo_cog_ey = 0.
            for i in xrange(self.bgohits.fEnergy.size()):
                e = self.bgohits.fEnergy[i]
		x = self.bgohits.GetHitX(i);
		y = self.bgohits.GetHitY(i);
		z = self.bgohits.GetHitZ(i);
		if self.bgohits.IsHitMeasuringX(i):
                    bgo_cog_x  += e*x
                    bgo_cog_zx += e*z
                    bgo_cog_ex += e
		else:
                    bgo_cog_y  += e*y
                    bgo_cog_zy += e*z
                    bgo_cog_ey += e
            bgo_cog_x  = bgo_cog_x  /  bgo_cog_ex 
            bgo_cog_y  = bgo_cog_y  /  bgo_cog_ey 
            bgo_cog_zx = bgo_cog_zx /  bgo_cog_ex 
            bgo_cog_zy = bgo_cog_zy /  bgo_cog_ey 


            for j in xrange(self.stktracks.GetLast()+1):
                track_tmp = self.stktracks.ConstructedAt(j)
                #self.__fill_histogrms__beforeselection__(track)
                if track_tmp.getNhitXY()<self.MIN_GOOG_XY_HITS: continue
                #slopex_good = True
                #slopey_good = True
                #if math.fabs(track.getTrackParams().getSlopeX()) > self.track_tan_cut_max: slopex_good = False
                #if math.fabs(track.getTrackParams().getSlopeX()) < self.track_tan_cut_min: slopex_good = False
                #if math.fabs(track.getTrackParams().getSlopeY()) > self.track_tan_cut_max: slopey_good = False
                #if math.fabs(track.getTrackParams().getSlopeY()) < self.track_tan_cut_min: slopey_good = False
                #if not slopex_good and not slopey_good: continue
                slope_tmp =  math.sqrt(track_tmp.getTrackParams().getSlopeX()**2 + track_tmp.getTrackParams().getSlopeY()**2)
                if slope_tmp < self.track_tan_cut_min: continue
                if slope_tmp > self.track_tan_cut_max: continue
                if self.exclude_cracks and self.is_crack_track(track_tmp): continue

                # some bgo-track match cuts
                stk_incl_x = track_tmp.getTrackParams().getSlopeX()
                stk_incl_y = track_tmp.getTrackParams().getSlopeY()
                stk_x = track_tmp.getImpactPoint().x()
                stk_y = track_tmp.getImpactPoint().y()
                stk_z = track_tmp.getImpactPoint().z()
                stk_intercept_x = stk_x - stk_incl_x * stk_z
                stk_intercept_y = stk_y - stk_incl_y * stk_z
                distance = (stk_intercept_x - bgo_intercept_x)**2 + (stk_intercept_y - bgo_intercept_y) **2
                angdistance = (stk_incl_x - bgo_incl_x)**2 + (stk_incl_y - bgo_incl_y)**2 
                cogdistance = (stk_x + stk_incl_x * (bgo_cog_zx - stk_z ) - bgo_cog_x)**2 + (stk_y + stk_incl_y * (bgo_cog_zy - stk_z ) - bgo_cog_y)**2 
                angdistance = angdistance**0.5
                cogdistance = cogdistance**0.5
                self.histo_distance.Fill(cogdistance)
                self.histo_angdistance.Fill(angdistance)
                #if BGO_STK_DISTANCE_MATCH is not None and distance    > BGO_STK_DISTANCE_MATCH **2 : continue             
                if BGO_STK_ANGULAR_MATCH  is not None and angdistance > BGO_STK_ANGULAR_MATCH:    continue
                if BGO_STK_COG_MATCH      is not None and cogdistance > BGO_STK_COG_MATCH:        continue
                if TRACK_CHISQCUT         is not None and track_tmp.getChi2() > TRACK_CHISQCUT:      continue

                # some additional track parameters
                chisqkalman_tmp = track_tmp.getChi2()
                inclfactor = (1+stk_incl_x**2 + stk_incl_y**2) ** 0.5
                maxclustercharge_tmp = 0.
                minclustercharge_tmp = 999999.
                
                xclusters_tmp = []
                covxx_tmp     = []
                for k in xrange(track_tmp.GetNPoints()):
                    xcluster = track_tmp.GetClusterX(k, self.stkclusters)
                    if not xcluster: continue
                    if self.is_cluster_bad_channel(xcluster): continue
                    if xcluster.getEnergy() > maxclustercharge_tmp: maxclustercharge_tmp = xcluster.getEnergy()
                    if xcluster.getEnergy() < minclustercharge_tmp: minclustercharge_tmp = xcluster.getEnergy()
                    xclusters_tmp.append(xcluster)
                    covxx_tmp.append(track_tmp.getCovXX(k))
                    
                yclusters_tmp = []
                covyy_tmp     = []
                for k in xrange(track_tmp.GetNPoints()):
                    ycluster = track_tmp.GetClusterY(k, self.stkclusters)
                    if not ycluster: continue
                    if self.is_cluster_bad_channel(ycluster): continue
                    if ycluster.getEnergy() > maxclustercharge_tmp: maxclustercharge_tmp = ycluster.getEnergy()
                    if ycluster.getEnergy() < minclustercharge_tmp: minclustercharge_tmp = ycluster.getEnergy()
                    yclusters_tmp.append(ycluster)
                    covyy_tmp.append(track_tmp.getCovYY(k))
         
                # some cluster cuts
                if TRACK_CLUSTER_ENERGY_MAX is not None and maxclustercharge_tmp > TRACK_CLUSTER_ENERGY_MAX * inclfactor: continue
                if CLUSTER_ENERGY_CUT is not None       and minclustercharge_tmp < CLUSTER_ENERGY_CUT: continue
                    
    
    
    
                #
                #  Filter tracks
                #
                #trackchisq_unnormalized = chisqx + chisqy               
                args_x_nonalign  = [cl.GetZ() for cl in xclusters_tmp]
                vals_x_nonalign  = [cl.GetX() for cl in xclusters_tmp]   
                #
                args_y_nonalign  = [cl.GetZ() for cl in yclusters_tmp]
                vals_y_nonalign  = [cl.GetY() for cl in yclusters_tmp]   
                #
                a_x_nonalign_tmp, b_x_nonalign_tmp = self.get_linear_fit(args_x_nonalign, vals_x_nonalign, getderivatives = False)
                res_x_nonalign = self.get_residuals(args_x_nonalign, vals_x_nonalign , a_x_nonalign_tmp, b_x_nonalign_tmp)         
                #
                a_y_nonalign_tmp, b_y_nonalign_tmp = self.get_linear_fit(args_y_nonalign, vals_y_nonalign, getderivatives = False)
                res_y_nonalign = self.get_residuals(args_y_nonalign, vals_y_nonalign , a_y_nonalign_tmp, b_y_nonalign_tmp)         
                #
                chisqx_tmp= np.sum([np.multiply(res[1],res[1]) for res in res_x_nonalign])
                chisqy_tmp = np.sum([np.multiply(res[1],res[1]) for res in res_y_nonalign])
                trackchisq_unnormalized_tmp = chisqx_tmp + chisqy_tmp
                #trackchisq_unnormalized = sum()
                if trackchisq_unnormalized_tmp>self.chi_sq_cut: #or chisqx>self.chi_sq_cut/2. or chisqy> self.chi_sq_cut/2.: 
                    return
    
                """
                # in case of duplicate - drop the event
                if self.__forbid_more_than_one__ and all_track:
                    self.__print_debug__( "Using not more than one good track! ==> exiting event!")
                    return None , None, None
                """


                clusterx_y_tmp = [a_y_nonalign_tmp*xclusters_tmp[i].GetZ() + b_y_nonalign_tmp  for i in xrange(len(xclusters_tmp))]
                clustery_x_tmp = [a_x_nonalign_tmp*yclusters_tmp[i].GetZ() + b_x_nonalign_tmp  for i in xrange(len(yclusters_tmp))]

                # found good track
                all_track.append(track_tmp)
                all_slope.append(slope_tmp)
                all_xclusters.append(xclusters_tmp)
                all_covxx.append(covxx_tmp)
                all_yclusters.append(yclusters_tmp)
                all_covyy.append(covyy_tmp)
                all_a_x_nonalign.append(a_x_nonalign_tmp)
                all_b_x_nonalign.append(b_x_nonalign_tmp)
                all_a_y_nonalign.append(a_y_nonalign_tmp)
                all_b_y_nonalign.append(b_y_nonalign_tmp)
                all_trackchisq_unnormalized.append(trackchisq_unnormalized_tmp)
                all_chisqx.append(chisqx_tmp)
                all_chisqy.append(chisqy_tmp)
                all_clusterx_y.append(clusterx_y_tmp)
                all_clustery_x.append(clustery_x_tmp)
                all_chisqkalman.append(chisqkalman_tmp)
                all_maxclustercharge.append(maxclustercharge_tmp)
                all_minclustercharge.append(minclustercharge_tmp)
                all_bgoangdist.append(angdistance)
                all_bgocogdist.append(cogdistance)
    

    
            # 2. check good tracks
            self.ngoodtracks.Fill(len(all_track))
            if len(all_track)==2:
                x,zx,y,zy = self.get_tracks_intersection(all_track[0],all_track[1])
                dz = zx-zy
                self.histo_vertex_dz.Fill(dz)
                self.histo_vertex_xz.Fill(x,zx)
                self.histo_vertex_yz.Fill(y,zy)
                if dz>=self.FOCUS_DZ_MIN and dz<=self.FOCUS_DZ_MAX:
                    self.histo_vertex_focus_dz.Fill(dz)
                    self.histo_vertex_focus_xz.Fill(x,zx)
                    self.histo_vertex_focus_yz.Fill(y,zy)
    
            if self.__forbid_more_than_one__ and len(all_track)>1:
                self.__print_debug__( "Using not more than one good track! ==> exiting event!")
                return
    
        # 3. If dump to file option - save info to file and do not process anything
        if self.__dump_data__ and not self.__run_on_compr_data__:
            for j in xrange(len(all_track)):
                y, z_y, x_y  = [NAN]*MAX_POSSIBLE_CLUSTERS, [NAN]*MAX_POSSIBLE_CLUSTERS, [NAN]*MAX_POSSIBLE_CLUSTERS
                x, z_x, y_x  = [NAN]*MAX_POSSIBLE_CLUSTERS, [NAN]*MAX_POSSIBLE_CLUSTERS, [NAN]*MAX_POSSIBLE_CLUSTERS
                covxx, covyy  = [NAN]*MAX_POSSIBLE_CLUSTERS, [NAN]*MAX_POSSIBLE_CLUSTERS
                ladx, lady = [0]*MAX_POSSIBLE_CLUSTERS, [0]*MAX_POSSIBLE_CLUSTERS
                for c in xrange(len(all_xclusters[j])):
                    x[c]     = all_xclusters[j][c].GetX()
                    z_x[c]   = all_xclusters[j][c].GetZ()
                    y_x[c]   = all_clusterx_y[j][c]
                    ladx[c]  = self.__get_ladder_ultimate__(all_xclusters[j][c], all_xclusters[j][c].GetX(), all_clusterx_y[j][c])
                    covxx[c] = all_covxx[j][c] 
                for c in xrange(len(all_yclusters[j])):
                    y[c]     = all_yclusters[j][c].GetY()
                    z_y[c]   = all_yclusters[j][c].GetZ()
                    x_y[c]   = all_clustery_x[j][c]
                    lady[c] = self.__get_ladder_ultimate__(all_yclusters[j][c], all_clustery_x[j][c], all_yclusters[j][c].GetY())
                    covyy[c] = all_covyy[j][c] 
                x, z_x, y_x, y, z_y, x_y  = np.array(x), np.array(z_x), np.array(y_x), np.array(y), np.array(z_y), np.array(x_y)
                ladx, lady = np.array(ladx), np.array(lady)
                covxx, covyy = np.array(covxx), np.array(covyy)
                
                vals = np.array([                   #  ====================== DEFINITION OF COMPRESSED DATA ANCILLARY PARAMETERS =====================
                    all_a_x_nonalign[j],            #      slope x a
                    all_b_x_nonalign[j],            #      slope x b
                    all_a_y_nonalign[j],            #      slope y a
                    all_b_y_nonalign[j],            #      slope y b
                    all_trackchisq_unnormalized[j], #      non normalized Chi-sqaure
                    all_chisqx[j],                  #      squared residuals x
                    all_chisqy[j],                  #      squared residuals y
                    all_ntracks,                    #      N ALL TRACKS IN EVENT
                    all_nclusters,                  #      N ALL CLUSTERS IN EVENT 
                    all_bgoenergy,                  #      BGO ENERGY
                    all_chisqkalman[j],             #      Track: Chi-square of Klaman track
                    all_maxclustercharge[j],        #      Track: max cluster charge 
                    all_minclustercharge[j],        #      Track: min cluster charge 
                    all_bgoangdist[j],              #      Track: ANGULAR MATCH TO BGO 
                    all_bgocogdist[j]               #      Track: DISTANCE MATCH TO BGO
                    ])                              #  ====================== DEFINITION OF COMPRESSED DATA ANCILLARY PARAMETERS =====================

                self.data_dump_all["x"].append(x)
                self.data_dump_all["z_x"].append(z_x)
                self.data_dump_all["y_x"].append(y_x)
                self.data_dump_all["ladx"].append(ladx)
                self.data_dump_all["y"].append(y)
                self.data_dump_all["z_y"].append(z_y)
                self.data_dump_all["x_y"].append(x_y)
                self.data_dump_all["lady"].append(lady)
                self.data_dump_all["covxx"].append(covxx)
                self.data_dump_all["covyy"].append(covyy)
                self.data_dump_all["vals"].append(vals)
                
        
        
        

        # 4a. One good track per event found - inspect it
        if not self.__run_on_compr_data__ and not self.__no_alignment__:
            for j in xrange(len(all_track)):
                track = all_track[j]
                xclusters = all_xclusters[j]
                covxx = all_covxx[j]
                yclusters = all_yclusters[j]
                covyy = all_covyy[j]
                a_x_nonalign = all_a_x_nonalign[j]
                b_x_nonalign = all_b_x_nonalign[j]
                a_y_nonalign = all_a_y_nonalign[j]
                b_y_nonalign = all_b_y_nonalign[j]
                trackchisq_unnormalized = all_trackchisq_unnormalized[j]
                chisqx = all_chisqx[j]
                chisqy = all_chisqy[j]
                clusterx_y = all_clusterx_y[j]
                clustery_x = all_clustery_x[j]
                #clusterx_y = [a_y_nonalign*xclusters[i].GetZ() + b_y_nonalign  for i in xrange(len(xclusters))]
                #clustery_x = [a_x_nonalign*yclusters[i].GetZ() + b_x_nonalign  for i in xrange(len(yclusters))]
    
                angle = (a_y_nonalign**2 + a_x_nonalign**2)**0.5
                self.histo_track_unnormalized_chisq_for_derivatives.Fill(trackchisq_unnormalized) 
                self.histo_track_kalman_chisq_for_derivatives.Fill(track.getChi2())
                self.histo_track_angle_x.Fill(a_x_nonalign)
                self.histo_track_angle_y.Fill(a_y_nonalign)
                self.histo_track_angle.Fill(angle)
                self.histo_track_anglex_chisqx.Fill(a_x_nonalign, chisqx)
                self.histo_track_angley_chisqy.Fill(a_y_nonalign, chisqy)
                self.histo_track_angle_chisq.Fill(angle,trackchisq_unnormalized)


                # fill cluster information
                a_total = math.sqrt(a_x_nonalign**2 + a_y_nonalign**2)
                for cluster in xclusters:                
                    #signal = sum([cluster.GetAdcValue(v,self.stkadccounts) for v in xrange(cluster.getNstrip())])
                    signalvacorr = sum([cluster.GetSignal(v) for v in xrange(cluster.getNstrip())])
                    self.histo_incline_nstrips_2D_x.Fill(a_x_nonalign,cluster.getNstrip())
                    #self.histo_incline_charge_2D_x.Fill(a_x_nonalign,signal)
                    self.histo_incline_chargevacorrected_2D_x.Fill(a_x_nonalign,signalvacorr) 
                    anglei = self.get_matching_angle(a_x_nonalign)
                    if anglei>=0:
                        self.histo_incline_nstrips_x[anglei].Fill(cluster.getNstrip())
                        #self.histo_incline_charge_x[anglei].Fill(signal)
                        self.histo_incline_chargevacorr_x[anglei].Fill(signalvacorr)


                for cluster in yclusters:                
                    #signal = sum([cluster.GetAdcValue(v,self.stkadccounts) for v in xrange(cluster.getNstrip())])
                    signalvacorr = sum([cluster.GetSignal(v) for v in xrange(cluster.getNstrip())])
                    self.histo_incline_nstrips_2D_y.Fill(a_y_nonalign,cluster.getNstrip())
                    #self.histo_incline_charge_2D_y.Fill(a_y_nonalign,signal)
                    self.histo_incline_chargevacorrected_2D_y.Fill(a_y_nonalign,signalvacorr)
                    anglei = self.get_matching_angle(a_y_nonalign)
                    if anglei>=0:
                        self.histo_incline_nstrips_y[anglei].Fill(cluster.getNstrip())
                        #self.histo_incline_charge_y[anglei].Fill(signal)
                        self.histo_incline_chargevacorr_y[anglei].Fill(signalvacorr)

                for cluster in xclusters + yclusters:
                    #signal = sum([cluster.GetAdcValue(v,self.stkadccounts) for v in xrange(cluster.getNstrip())])
                    signalvacorr = sum([cluster.GetSignal(v) for v in xrange(cluster.getNstrip())])
                    self.histo_incline_nstrips_2D.Fill(a_total,cluster.getNstrip())
                    #self.histo_incline_charge_2D.Fill(a_total,signal)
                    self.histo_incline_chargevacorrected_2D.Fill(a_total,signalvacorr) 
                    anglei = self.get_matching_angle(a_total)
                    if anglei>=0:
                        self.histo_incline_nstrips[anglei].Fill(cluster.getNstrip())
                        #self.histo_incline_charge[anglei].Fill(signal)
                        self.histo_incline_chargevacorr[anglei].Fill(signalvacorr)



            
                #@ Get ChiSq derivatives for track and inspect track
                #deriv_for_track, tracks_for_lad, trackchisq = self.__process_tracks_get_chisq_derivatives__(track, xclusters, yclusters,  covxx, covyy, clusterx_y, clustery_x)
                deriv_for_track, tracks_for_lad, trackchisq = self.__process_tracks_get_chisq_derivatives_RECFILE__(track, xclusters, yclusters,  covxx, covyy, clusterx_y, clustery_x)
                if deriv_for_track is not None and trackchisq < self.aligntrack_chi_sq_cut:
           
                    self.residuals_cut_passed = True
    
                    # inspect track 
                    if not self.__no_inspect__:
                        self.__process_tracks_inspect__(track, xclusters, yclusters,  covxx, covyy,  clusterx_y, clustery_x)
                        self.histo_alignedtrack_chisq.Fill(trackchisq)

                    # append derivatives
                    if self.residuals_cut_passed:
                        self.dchish_dpars["derivatives"] = np.add(self.dchish_dpars["derivatives"], deriv_for_track)
                        self.dchish_dpars["ntracks"] += 1
                        self.dchish_dpars["chisq"]   += trackchisq
                        self.tracks_for_lad = np.add(self.tracks_for_lad, tracks_for_lad)
            
            
        # 4b. One good track per event found - inspect it
        if self.__run_on_compr_data__ and not self.__no_alignment__:
            deriv_for_track, tracks_for_lad, trackchisq =  self.__process_tracks_get_chisq_derivatives_COMPRFILE__()


            # some track sellection condition
            track_inclx     = self.data_val[0]
            track_incly     = self.data_val[2]
            event_ntracktot = self.data_val[7]
            event_nclusters = self.data_val[8]
            event_bgoenergy = self.data_val[9]
            track_chisqkalm = self.data_val[10]
            track_maxclchrg = self.data_val[11]
            track_minclchrg = self.data_val[12]
            track_bgoangdst = self.data_val[13]
            track_bgocogdst = self.data_val[14]
            inclfactor = (1 + track_inclx**2 + track_incly**2) ** 0.5

            # -------
            trackcond = True
            # -------
            if MAX_TRACKS_IN_EVENT is not None and event_ntracktot > MAX_TRACKS_IN_EVENT: trackcond = False
            # -------
            if BGO_REC_MAX is not None and event_bgoenergy> BGO_REC_MAX: trackcond = False
            if BGO_REC_MIN is not None and event_bgoenergy< BGO_REC_MIN: trackcond = False
            # -------
            self.histo_nclustersbefore.Fill(event_nclusters)
            if STK_MAX_CLUSTERS is not None and event_nclusters > STK_MAX_CLUSTERS : trackcond = False
            self.histo_nclustersafter.Fill(event_nclusters)
            # -------
            if BGO_STK_ANGULAR_MATCH  is not None and track_bgoangdst > BGO_STK_ANGULAR_MATCH: trackcond = False
            self.histo_angdistance.Fill(track_bgoangdst)
            # -------
            if BGO_STK_COG_MATCH is not None and track_bgocogdst > BGO_STK_COG_MATCH: trackcond = False
            self.histo_distance.Fill(track_bgocogdst)
            # -------
            if TRACK_CHISQCUT is not None and track_chisqkalm > TRACK_CHISQCUT: trackcond = False
            self.histo_track_kalman_chisq_for_derivatives.Fill(track_chisqkalm)
            # -------
            if TRACK_CLUSTER_ENERGY_MAX is not None and track_maxclchrg > TRACK_CLUSTER_ENERGY_MAX * inclfactor: trackcond = False
            if CLUSTER_ENERGY_CUT is not None and track_minclchrg < CLUSTER_ENERGY_CUT: trackcond = False
            # -------

            
            
            if trackcond and deriv_for_track is not None and trackchisq < self.aligntrack_chi_sq_cut:
                

                self.residuals_cut_passed = True

                # inspect track 
                if not self.__no_inspect__:
                    #self.__process_tracks_inspect__(track, xclusters, yclusters,  covxx, covyy,  clusterx_y, clustery_x)

                    goodx = [val for val in self.data_x if not self.__is_nan__(val)] 
                    goody = [val for val in self.data_y if not self.__is_nan__(val)] 
                    if len(goodx)>=self.MIN_GOOG_XY_HITS and len(goody)>=self.MIN_GOOG_XY_HITS:
                        self.__process_tracks_inspect__(None, goodx, goody, self.data_covxx, self.data_covyy, None, None)
                        self.histo_alignedtrack_chisq.Fill(trackchisq)

                # append derivatives
                if self.residuals_cut_passed:
                    self.dchish_dpars["derivatives"] = np.add(self.dchish_dpars["derivatives"], deriv_for_track)
                    self.dchish_dpars["ntracks"] += 1
                    self.dchish_dpars["chisq"]   += trackchisq
                    self.tracks_for_lad = np.add(self.tracks_for_lad, tracks_for_lad)

            
            
                             
            
    def __get_ladder_ultimate__(self, cluster, x=None, y=None):
        if not self.__all_sensors__:
            return cluster.getLadderHardware()
        if cluster.isX():
            return cluster.getLadderHardware()*NSENSORS + self.__get_silicon_sensor_index__(y) 
        else:
            return cluster.getLadderHardware()*NSENSORS + self.__get_silicon_sensor_index__(x) 


    def get_tracks_intersection(self, track1, track2):
        track1_x = track1.getImpactPoint().x()
        track1_y = track1.getImpactPoint().y()
        track1_z = track1.getImpactPoint().z()
        track1_incl_x = track1.getTrackParams().getSlopeX()
        track1_incl_y = track1.getTrackParams().getSlopeY()
        track2_x = track2.getImpactPoint().x()
        track2_y = track2.getImpactPoint().y()
        track2_z = track2.getImpactPoint().z()
        track2_incl_x = track2.getTrackParams().getSlopeX()
        track2_incl_y = track2.getTrackParams().getSlopeY()
        #@
        ax1 = track1_incl_x
        bx1 = track1_x - track1_incl_x * track1_z
        ay1 = track1_incl_y
        by1 = track1_y - track1_incl_y * track1_z
        #@
        ax2 = track2_incl_x
        bx2 = track2_x - track2_incl_x * track2_z
        ay2 = track2_incl_y
        by2 = track2_y - track2_incl_y * track2_z
        #@
        z_cross_x = (bx2 - bx1) * 1.0 / (ax1 - ax2) 
        z_cross_y = (by2 - by1) * 1.0 / (ay1 - ay2) 
        #@
        cross_x = ax1 * z_cross_x + bx1
        cross_y = ay1 * z_cross_y + by1
        #@
        return cross_x,z_cross_x,cross_y,z_cross_y
      
    def __is_nan__(self,val):
        return val < NAN*0.1
          
    def __process_tracks_get_chisq_derivatives_COMPRFILE__(self):
        self.__measuretime_start__("test0")
        argsx, valsx, ladx, argsy, valsy, lady = [], [], [], [], [], []
        #normalizex, normalizey = [], []
        for x_tmp, z_x_tmp, y_x_tmp, ladx_tmp, covxx_tmp in izip(
                 self.data_x,
                 self.data_z_x,
                 self.data_y_x,
                 self.data_ladx,
                 self.data_covxx):
            if self.__is_nan__(x_tmp): continue

            argsx.append(self.get_aligned_z( {"lad":ladx_tmp, "z":z_x_tmp} , x_tmp, y_x_tmp, None))
            valsx.append(self.get_aligned_x( {"lad":ladx_tmp, "x":x_tmp},           y_x_tmp, None))
            ladx.append(ladx_tmp)  
            #normalizex.append(1.0) 

        for y_tmp, z_y_tmp, x_y_tmp, lady_tmp, covyy_tmp in izip(
                 self.data_y,
                 self.data_z_y,
                 self.data_x_y,
                 self.data_lady,
                 self.data_covyy):
            if self.__is_nan__(y_tmp): continue
            argsy.append(self.get_aligned_z( {"lad":lady_tmp, "z":z_y_tmp} , x_y_tmp, y_tmp, None))
            valsy.append(self.get_aligned_y( {"lad":lady_tmp, "y":y_tmp}   , x_y_tmp,        None))
            lady.append(lady_tmp)  
            #normalizey.append(1.0) 
        """
        print 
        print    "argsx: ",argsx
        print    "valsx: ",valsx
        print    "ladx:  ",ladx 
        print    "argsy: ",argsy
        print    "valsy: ",valsy
        print    "lady:  ",lady
        print    "normalizex: ",normalizex 
        print    "normalizey: ",normalizey
        print 
        raise SystemExit
        """
        argsx, valsx, ladx, argsy, valsy, lady = np.array(argsx), np.array(valsx), np.array(ladx), np.array(argsy), np.array(valsy), np.array(lady)

        if len(argsx)<self.MIN_GOOG_XY_HITS: return None , None, None
        if len(argsy)<self.MIN_GOOG_XY_HITS: return None , None, None
        self.__measuretime_stop__("test0")
        
        if not self.__is_ntracks_file__:
            normalizex = None
            normalizey = None
        else:
            normalizex = [ self.normalize_n_tracks[i]  for i in ladx]
            normalizey = [ self.normalize_n_tracks[i]  for i in lady]
            


        return self.__process_tracks_get_chisq_derivatives__(argsx,valsx,ladx,argsy,valsy,lady,normalizex,normalizey)

 
        """
        justtest= [self.data_x, 
                 self.data_z_x, 
                 self.data_ladx, 
                 self.data_covxx, 
                 self.data_y,
                 self.data_z_y, 
                 self.data_lady, 
                 self.data_covyy, 
                 self.data_val] 
         """
     
       
    def __process_tracks_get_chisq_derivatives_RECFILE__(self,track, xclusters, yclusters,  covxx, covyy,  clusterx_y, clustery_x):
        self.__measuretime_start__("test0")
        if len(xclusters)<self.MIN_GOOG_XY_HITS: return None , None, None
        if len(yclusters)<self.MIN_GOOG_XY_HITS: return None , None, None
        
        
        argsx  = [self.get_aligned_z(xclusters[i] ,xclusters[i].GetX(), clusterx_y[i], track)  for i in xrange(len(xclusters))]
        valsx  = [self.get_aligned_x(xclusters[i] ,                     clusterx_y[i], track)  for i in xrange(len(xclusters))]
        #ladx   = [c.getLadderHardware() for c in xclusters] 
        ladx = [self.__get_ladder_ultimate__(xclusters[i], xclusters[i].GetX(), clusterx_y[i]) for i in xrange(len(xclusters))]
        
        #argsy = [self.get_aligned_z(c) for c in yclusters] 
        #valsy = [self.get_aligned_y(c) for c in yclusters]
        
        argsy = [self.get_aligned_z(yclusters[i] ,clustery_x[i], yclusters[i].GetY(),track)  for i in xrange(len(yclusters))]
        valsy = [self.get_aligned_y(yclusters[i] ,clustery_x[i]                     ,track)  for i in xrange(len(yclusters))]
        #lady   = [c.getLadderHardware() for c in yclusters]
        lady = [self.__get_ladder_ultimate__(yclusters[i], clustery_x[i], yclusters[i].GetY()) for i in xrange(len(yclusters))]
        
        #normalizex = [1.0 / np.power(self.normalize_n_tracks[cl.getLadderHardware()],2)  for cl in xclusters]
        if self.__is_ntracks_file__:
            #normalizex = [1.0 / np.power(self.normalize_n_tracks[  self.__get_ladder_ultimate__(xclusters[i], xclusters[i].GetX() , clusterx_y[i])         ],2)  for i in xrange(len(xclusters))]
            #normalizey = [1.0 / np.power(self.normalize_n_tracks[  self.__get_ladder_ultimate__(yclusters[i], clustery_x[i]       , yclusters[i].GetY())   ],2)  for i in xrange(len(yclusters))]
            normalizex = [self.normalize_n_tracks[  self.__get_ladder_ultimate__(xclusters[i], xclusters[i].GetX() , clusterx_y[i])         ]  for i in xrange(len(xclusters))]
            normalizey = [self.normalize_n_tracks[  self.__get_ladder_ultimate__(yclusters[i], clustery_x[i]       , yclusters[i].GetY())   ]  for i in xrange(len(yclusters))]
        else:
            normalizex = None
            normalizey = None
        #normalizey = [1.0 / np.power(self.normalize_n_tracks[cl.getLadderHardware()],2)  for cl in yclusters]

        argsx, valsx, ladx, argsy, valsy, lady = np.array(argsx), np.array(valsx), np.array(ladx), np.array(argsy), np.array(valsy), np.array(lady)

        return self.__process_tracks_get_chisq_derivatives__(argsx,valsx,ladx,argsy,valsy,lady,normalizex,normalizey)
        self.__measuretime_stop__("test0")




        
    def __process_tracks_get_chisq_derivatives__(self,argsx,valsx,ladx,argsy,valsy,lady,normalizex,normalizey):
        
        ax, bx, dax_darg, dax_dval, dbx_darg, dbx_dval, residualsx, chisqx = self.get_linear_fit(argsx, valsx, getderivatives = True)
        ay, by, day_darg, day_dval, dby_darg, dby_dval, residualsy, chisqy = self.get_linear_fit(argsy, valsy, getderivatives = True)

        ix = self.get_matching_algorithmangle(ax)
        iy = self.get_matching_algorithmangle(ay)

        dchish_dpars = [[0. for j in xrange(N_DOF_LADDER)] for i in xrange(self.N_STK_LADDERS)]
        tracks_for_lad = [[0 for angle1 in ANGLES_FORALGORITHM] for i in xrange(self.N_STK_LADDERS)]
        self.__measuretime_start__("test1")
        
        alllads = []
        for i in ladx:
            if i not in alllads: alllads.append(i)
        for i in lady:
            if i not in alllads: alllads.append(i)

        """
        # DEBUG information
        for i in xrange(len(normalizex)):
            thenorm = normalizex[i]
            if thenorm[ix][iy]: continue
            thelad  = ladx[i] 
            print "WARNING: no track normaliztion for inclination for ladder ", thelad, " ix, iy = ",ix,iy, "  angles = ", ax, ay

        # DEBUG information
        for i in xrange(len(normalizey)):
            thenorm = normalizey[i]
            if thenorm[ix][iy]: continue
            thelad  = lady[i] 
            print "WARNING: no track normaliztion for inclination for ladder ", thelad, " ix, iy = ",ix,iy, "  angles = ", ax, ay
        """



        #for lad in xrange(self.N_STK_LADDERS):
        #for lad in self.ALL_LADS:
        for lad in alllads:
            #if lad not in ladx and lad not in lady: continue

            if lad in ladx:
                tracks_for_lad[lad][ix] = 1
            elif lad in lady:
                tracks_for_lad[lad][iy] = 1
                
            #for par_index in xrange(N_DOF_LADDER):
            for par_index in self.ALL_PARS_DOF:
                




                dvalx_dpar =  self.dval_dpar_xtrack(ladx, valsx, valsy, lad, par_index)
                dargx_dpar =  self.darg_dpar(ladx, valsx, valsy, lad, par_index)
                dax_dpar = self.da_dpar(dax_dval, dvalx_dpar, dax_darg, dargx_dpar)
                dbx_dpar = self.db_dpar(dbx_dval, dvalx_dpar, dbx_darg, dargx_dpar)
                
                
                dresx_dpar =  np.subtract(
                                          dvalx_dpar,
                                          np.multiply(dax_dpar, argsx)
                                          )
                dresx_dpar = np.subtract(
                                         dresx_dpar,
                                         np.multiply(ax,dargx_dpar)
                                         )
                dresx_dpar = np.subtract(
                                         dresx_dpar,
                                         dbx_dpar
                                         )
                
                #dresx_dpar =  [dvalx_dpar_val - dax_dpar * argsx_val for dvalx_dpar_val,argsx_val in izip(dvalx_dpar, argsx)]
                #dresx_dpar =  [dresx_dpar_val - ax*dargx_dpar_val for dresx_dpar_val,dargx_dpar_val in izip(dresx_dpar,dargx_dpar)]
                #dresx_dpar =  [dresx_dpar_val - dbx_dpar for dresx_dpar_val in dresx_dpar]
                
                """ 
                print 
                print
                print "dvalx_dpar=", dvalx_dpar  
                print "dax_dpar=", dax_dpar  
                print "argsx=", argsx  
                print 
                print
                print 
                print
                print "dresx_dpar=", dresx_dpar  
                print "ax=", ax  
                print "dargx_dpar=", dargx_dpar  
                print 
                print
                print
                print "dresx_dpar=", dresx_dpar
                print "dbx_dpar=", dbx_dpar  
                print 
                print
                raise SystemExit
                """
                

                
                
                
                dvaly_dpar =  self.dval_dpar_ytrack(lady, valsx, valsy, lad, par_index)
                dargy_dpar =  self.darg_dpar(lady, valsx, valsy, lad, par_index)
                day_dpar = self.da_dpar(day_dval, dvaly_dpar, day_darg, dargy_dpar)
                dby_dpar = self.db_dpar(dby_dval, dvaly_dpar, dby_darg, dargy_dpar)
                
                
                dresy_dpar =  np.subtract(
                                          dvaly_dpar,
                                          np.multiply(day_dpar, argsy)
                                          )
                dresy_dpar = np.subtract(
                                         dresy_dpar,
                                         np.multiply(ay,dargy_dpar)
                                         )
                dresy_dpar = np.subtract(
                                         dresy_dpar,
                                         dby_dpar
                                         )
                
                
                
                #dresy_dpar =  [dvaly_dpar_val - day_dpar* argsy_val for dvaly_dpar_val,argsy_val in izip(dvaly_dpar,argsy)]
                #dresy_dpar =  [dresy_dpar_val - ay*dargy_dpar_val for dresy_dpar_val,dargy_dpar_val in izip(dresy_dpar,dargy_dpar)]
                #dresy_dpar =  [dresy_dpar_val - dby_dpar for dresy_dpar_val in dresy_dpar]
               
                
                
                
                #dchisc_dpar =  sum(
                #                  [residualsx[i] * dresx_dpar[i] * normalizex[i] for i in xrange(len(residualsx))]
                #                  ) +  sum(
                #                  [residualsy[i] * dresy_dpar[i] * normalizey[i] for i in xrange(len(residualsy))]
                #                  )
                
                if normalizex is not None and normalizey is not None:
                    dchisc_dpar = np.sum(
                                         np.multiply(np.multiply(residualsx,dresx_dpar),[self.__wrap_normaliztions_with_angles__(norm,ix) for norm in normalizex])
                              ) + np.sum(
                                         np.multiply(np.multiply(residualsy,dresy_dpar),[self.__wrap_normaliztions_with_angles__(norm,iy) for norm in normalizey])
                              )
                else:
                    
                    dchisc_dpar = np.sum(np.multiply(residualsx,dresx_dpar)) + np.sum(np.multiply(residualsy,dresy_dpar) )
                    #dchisc_dpar = sum([a*b for a,b in izip(residualsx,dresx_dpar)]) + sum([a*b for a,b in izip(residualsy,dresy_dpar)])
                   
                                          
                    
                #dchisc_dpar = sum(np.multiply(residualsx, dresx_dpar)) + sum(np.multiply(residualsy, dresy_dpar))
               
                dchish_dpars[lad][par_index] = dchisc_dpar  
                
                
                
        self.__measuretime_stop__("test1")
        
        if normalizex is not None and normalizey is not None:
            trackchisq = sum(
                             [residualsx[i] * residualsx[i] * self.__wrap_normaliztions_with_angles__(normalizex[i],ix) for i in xrange(len(residualsx))]
                            ) + sum(
                             [residualsy[i] * residualsy[i] * self.__wrap_normaliztions_with_angles__(normalizey[i],iy) for i in xrange(len(residualsy))]
                            )
        else:
            trackchisq = np.sum(np.multiply(residualsx,residualsx)) + np.sum(np.multiply(residualsy,residualsy))


          
        
                
        return dchish_dpars, tracks_for_lad, trackchisq


    def __wrap_normaliztions_with_angles__(self,norm,i):
        if norm[i]:
            return 1. / norm[i]
        else:
            return 1.
            
            
    #
    # Functions to be used before any track selection \ correction
    #
    def __fill_histogrms__beforeselection__(self, track):
        """
        This is not very elegant part ... 
        """
        allz = [-46, -79, -111, -144, -176, -210]
        dz   = 5
        dr   = 1
        COMPAREPOINT = -100
        
        
        if track.getNhitXY()<3: return
        if track.GetNPoints()<self.N_STK_LAYERS:
            return
        
        
        for i in xrange(track.GetNPoints()):
            #clusterx = track.GetClusterX(i, self.stkclusters)
            #clustery = track.GetClusterY(i, self.stkclusters)
            x = track.getHitX(i)
            y = track.getHitX(i)
            for m in xrange(self.stkclusters.GetLast()+1):
                cluster = self.stkclusters.ConstructedAt(m)
                if math.fabs(cluster.GetZ() - allz[i]) > dz: continue
                if cluster.isX() and math.fabs(cluster.GetX()-x)< dr:
                    self.nxcluster_vs_y[i].Fill(y)
                if cluster.isY() and math.fabs(cluster.GetY()-y)< dr:
                    self.nycluster_vs_x[i].Fill(x)
            
            
            
        # inspect crack in the middle of tray
        track_x_positivey = []
        track_x_negativey = []
        track_y_positivex = []
        track_y_negativex = []    
        for i in xrange(track.GetNPoints()):
            clusterx = track.GetClusterX(i,self.stkclusters)
            clustery = track.GetClusterY(i,self.stkclusters)
            if not clusterx or not clustery: continue
            if clusterx.GetX()>0:
                track_y_positivex.append(clustery)
            else:
                track_y_negativex.append(clustery)
            
            if clustery.GetY()>0:
                track_x_positivey.append(clusterx)
            else:
                track_x_negativey.append(clusterx)
                
              
                
        if len(track_x_positivey)>1 and len(track_x_negativey)>1:
            a1, b1 = self.get_linear_fit([self.get_aligned_z(cl,None,None,track) for cl in track_x_positivey], 
                                    [self.get_aligned_x(cl,None,track) for cl in track_x_positivey] )
            a2, b2 = self.get_linear_fit([self.get_aligned_z(cl,None,None,track) for cl in track_x_negativey], 
                                    [self.get_aligned_x(cl,None,track) for cl in track_x_negativey] )
            val1 = a1*COMPAREPOINT + b1
            val2 = a2*COMPAREPOINT + b2 
            self.centralcrack_x.Fill(val1-val2)
            
        if len(track_y_positivex)>1 and len(track_y_negativex)>1:
            a1, b1 = self.get_linear_fit([self.get_aligned_z(cl,None,None,track) for cl in track_y_positivex], 
                                    [self.get_aligned_y(cl,None,track) for cl in track_y_positivex] )
            a2, b2 = self.get_linear_fit([self.get_aligned_z(cl,None,None,track) for cl in track_y_negativex], 
                                    [self.get_aligned_y(cl,None,track) for cl in track_y_negativex] )
            val1 = a1*COMPAREPOINT + b1
            val2 = a2*COMPAREPOINT + b2
            self.centralcrack_y.Fill(val1-val2)
    #
    #    
    def __get_track_cluster__(self, track, index, direction):
        if direction.lower()=="x":
            return track.GetClusterX(index, self.stkclusters)
        elif direction.lower()=="y":
            return track.GetClusterY(index, self.stkclusters)
        else:
            raise Exception("[__get_track_cluster__] Unrecognized situation")
    #
    #   
    def __get_cluster_coordinate__(self, cluster, direction):
        if direction.lower()=="x":
            return cluster.GetX()
        elif direction.lower()=="y":
            return cluster.GetY()
        else:
            raise Exception("[__get_cluster_coordinate__] Unrecognized situation")
    #
    #
    def __is_crack_track__(self, track, direction):
        for k in xrange(track.GetNPoints()-1):
            cluster = self.__get_track_cluster__(track, k, direction)
            if not cluster: continue
            for m in xrange(k, track.GetNPoints()):
                next_cluster = self.__get_track_cluster__(track, m, direction)
                if not next_cluster: continue
                if math.fabs(self.__get_cluster_coordinate__(cluster, direction) - 
                             self.next_cluster(cluster, direction)) > self.TRACK_CRACK_REGION_DISTANCE: continue
                condition1 =  cluster.getIndex1() >= self.TRACK_CRACK_REGION_NSTRIPS 
                condition2 =  cluster.getIndex1() <  self.N_RO_STRIPS - self.TRACK_CRACK_REGION_DISTANCE + 1 -  cluster.getNstrip()
                if condition1 and condition2: continue
                
                
                
                
                
                
        
        
    def __process_tracks_align_plane__(self, track,  xclusters, yclusters,  covxx, covyy, planefromcalorimeter):
        #
        #  planefromcalorimeter = 0,1,2, ..., 5
        #
        if len(xclusters)>=self.MIN_GOOG_XY_HITS: return 
        if len(yclusters)>=self.MIN_GOOG_XY_HITS: return
        

    def __process_tracks_inspect__(self, track, xclusters, yclusters, covxx, covyy,  clusterx_y, clustery_x):
        if self.__unbiased__:
            self.__process_tracks_inspect_unbiased__(track, xclusters, yclusters, covxx, covyy,  clusterx_y, clustery_x)
        else:
            self.__process_tracks_inspect_biased__(track, xclusters, yclusters, covxx, covyy,  clusterx_y, clustery_x)

    def get_matching_angle(self,a):
        incl = abs(a)
        for i in xrange(len(ANGLES)):
            if incl > ANGLES[i]: continue
            return i
        return -1

    def get_matching_algorithmangle(self,a):
        incl = abs(a)
        for i in xrange(len(ANGLES_FORALGORITHM)):
            if incl > ANGLES_FORALGORITHM[i]: continue
            return i
        return -1

    
    def get_matching_histos_y(self, ay):
        tng = abs(ay)
        for i in xrange(len(ANGLES)+1):
            if i==0:
                if tng > ANGLES[0]: continue
            elif i<len(ANGLES):
                if tng<ANGLES[i-1] or tng>ANGLES[i]: continue
            yield self.residy[i], self.residy_vs_x[i], self.residy_vs_y[i], self.residy_vs_y_posv[i], self.residy_vs_y_negv[i], self.chisqy[i], self.covyy[i], self.resid_all_lad_y[i]
            
        
    def get_matching_histos_x(self, ax):
        tng = abs(ax)
        for i in xrange(len(ANGLES)+1):
            if i==0:
                if tng > ANGLES[0]: continue
            elif i<len(ANGLES):
                if tng<ANGLES[i-1] or tng>ANGLES[i]: continue
            yield self.residx[i], self.residx_vs_x[i], self.residx_vs_x_posv[i], self.residx_vs_x_negv[i], self.residx_vs_y[i], self.chisqx[i], self.covxx[i], self.resid_all_lad_x[i]



    def __process_tracks_inspect_unbiased__(self, track, xclusters, yclusters, covxx, covyy,  clusterx_y, clustery_x):
        if len(xclusters)==MAX_POSSIBLE_CLUSTERS and len(yclusters)==MAX_POSSIBLE_CLUSTERS:


            if self.__run_on_compr_data__:
                argsx  = [self.get_aligned_z({"lad":self.data_ladx[i], "z": self.data_z_x[i]}, self.data_x[i], self.data_y_x[i], None)  for i in xrange(len(self.data_x))]
                valsx  = [self.get_aligned_x({"lad":self.data_ladx[i], "x": self.data_x[i]},                   self.data_y_x[i], None)  for i in xrange(len(self.data_x))]
                clzx   = self.data_z_x
                ladidx = self.data_ladx
        
                argsy  = [self.get_aligned_z({"lad":self.data_lady[i], "z": self.data_z_y[i]}, self.data_x_y[i], self.data_y[i], None)  for i in xrange(len(self.data_y))]
                valsy  = [self.get_aligned_y({"lad":self.data_lady[i], "y": self.data_y[i]},   self.data_x_y[i],                 None)  for i in xrange(len(self.data_y))]
                clzy   = self.data_z_y
                ladidy = self.data_lady 
                
            else:
                argsx  = [self.get_aligned_z(xclusters[i], xclusters[i].GetX(), clusterx_y[i],track)  for i in xrange(len(xclusters))]
                valsx  = [self.get_aligned_x(xclusters[i],                      clusterx_y[i],track)  for i in xrange(len(xclusters))]
                clzx   = [c.GetZ() for c in xclusters]
                ladidx = [c.getLadderHardware() for c in xclusters] 
        
                argsy = [self.get_aligned_z(yclusters[i], clustery_x[i], yclusters[i].GetY(),track)  for i in xrange(len(yclusters))]
                valsy = [self.get_aligned_y(yclusters[i], clustery_x[i]                     ,track)  for i in xrange(len(yclusters))]
                clzy  = [c.GetZ() for c in yclusters]
                ladidy = [c.getLadderHardware() for c in yclusters] 

            ax, bx = self.get_linear_fit(argsx, valsx)
            ay, by = self.get_linear_fit(argsy, valsy)

            if not self.__run_on_compr_data__:
                self.trackchisqx.Fill(track.getChi2())

            #resid = self.get_residuals(argsx, valsx, ax, bx)
            #for i in xrange(len(resid)):
            for i in xrange(len(argsx)):
                #arg = resid[i][0]
                #dif = resid[i][1]
                clz = clzx[i]
                ladid = ladidx[i] 
                #histotag = "_z_%d_mm"%clz
                

                # add histograms
                for histo_residx, histo_residx_vs_x, histo_residx_vs_x_posv, histo_residx_vs_x_negv, histo_residx_vs_y, histo_chisqx, histo_covxx, histo_lad_res_x in self.get_matching_histos_x(ax):
                    if clz not in histo_residx: 
                        self.create_control_histos_x(clz, histo_residx, histo_residx_vs_x, histo_residx_vs_x_posv, histo_residx_vs_x_negv, histo_residx_vs_y, histo_chisqx, histo_covxx)
                    if ladid not in histo_lad_res_x:
                        self.create_control_histos_ladders_x(histo_lad_res_x, ladid)
 
                    
                    
                # Fill histograms
                arg = argsx[i]
                val = valsx[i]
                arguments = argsx[:i]+argsx[i+1:]
                values    =  valsx[:i]+valsx[i+1:]
                a, b = self.get_linear_fit(arguments, values)
                projectedx = arg * a + b
                dif = val  - projectedx

                
                # Apply track traightness cut
                if self.__unbiased_track_residuals_max__:
                    resx_unbiased = self.get_residuals(arguments, values , a, b)
                    if [res for res in resx_unbiased if abs(res[1])> self.__unbiased_track_residuals_max__] or abs(dif)>self.__unbiased_track_residuals_max__: 
                        self.residuals_cut_passed = False
                        continue

                # cut in other projection
                if self.__unbiased_track_residuals_otherprojection_max__:
                    resy_all = self.get_residuals(argsy, valsy , ay, by)
                    if [res for res in resy_all if abs(res[1])> self.__unbiased_track_residuals_otherprojection_max__]: 
                        self.residuals_cut_passed = False
                        continue
                   



                #positivey = True if xclusters[i].getLadderHardware()/NLADDERS_TRB in [2,3] else False
                positivey = True if ladidx[i]/ (NSENSORS if self.ALIGNSENSORS_ARG else 1) /NLADDERS_TRB in [2,3] else False
                
                for histo_residx, histo_residx_vs_x, histo_residx_vs_x_posv, histo_residx_vs_x_negv, histo_residx_vs_y, histo_chisqx, histo_covxx, histo_lad_res_x in self.get_matching_histos_x(ax):
                    histo_residx[clz].Fill(dif)
                    histo_lad_res_x[ladid].Fill(dif)
                    histo_residx_vs_x[clz].Fill(dif, ax*arg + bx)
                    if positivey: histo_residx_vs_x_posv[clz].Fill(dif, ax*arg + bx)
                    else:         histo_residx_vs_x_negv[clz].Fill(dif, ax*arg + bx)
                    histo_residx_vs_y[clz].Fill(dif, ay*arg + by)
                    histo_chisqx[clz].Fill(dif * dif / covxx[i])
                    histo_covxx[clz].Fill(covxx[i]) 
                
            
 
        #if len(yclusters)>=self.MIN_GOOG_XY_HITS:
            if not self.__run_on_compr_data__:
                self.trackchisqy.Fill(track.getChi2())

            #resid = self.get_residuals(argsy, valsy, ay, by)
            #for i in xrange(len(resid)):
            for i in xrange(len(argsy)):
                #arg = resid[i][0]
                #dif = resid[i][1]
                clz = clzy[i]
                ladid = ladidy[i] 
                #histotag = "_z_%d_mm"%clz
                
                #add histograms
                for histo_residy, histo_residy_vs_x, histo_residy_vs_y, histo_residy_vs_y_posv, histo_residy_vs_y_negv, histo_chisqy, histo_covyy, histo_lad_res_y in self.get_matching_histos_y(ay):
                    if clz not in histo_residy: 
                        self.create_control_histos_y(clz, histo_residy, histo_residy_vs_x, histo_residy_vs_y, histo_residy_vs_y_posv, histo_residy_vs_y_negv, histo_chisqy, histo_covyy)
                    if ladid not in histo_lad_res_y:
                        self.create_control_histos_ladders_y(histo_lad_res_y, ladid)
                
                arg = argsy[i]
                val = valsy[i]
                arguments = argsy[:i]+argsy[i+1:]
                values    =  valsy[:i]+valsy[i+1:]
                a, b = self.get_linear_fit(arguments, values)
                projectedy = arg * a + b
                dif = val  - projectedy

                # Apply track traightness cut
                if self.__unbiased_track_residuals_max__:
                    resy_unbiased = self.get_residuals(arguments, values , a, b)
                    if [res for res in resy_unbiased if abs(res[1])> self.__unbiased_track_residuals_max__] or abs(dif)>self.__unbiased_track_residuals_max__: 
                        self.residuals_cut_passed = False
                        continue

                # cut in other projection
                if self.__unbiased_track_residuals_otherprojection_max__:
                    resx_all = self.get_residuals(argsx, valsx , ax, bx)
                    if [res for res in resx_all if abs(res[1])> self.__unbiased_track_residuals_otherprojection_max__]: 
                        self.residuals_cut_passed = False
                        continue

                #positivex = True if yclusters[i].getLadderHardware()/NLADDERS_TRB in [0,1] else False
                positivex = True if ladidy[i]/(NSENSORS if self.ALIGNSENSORS_ARG else 1)/NLADDERS_TRB in [0,1] else False
                for histo_residy, histo_residy_vs_x, histo_residy_vs_y, histo_residy_vs_y_posv, histo_residy_vs_y_negv, histo_chisqy, histo_covyy, histo_lad_res_y in self.get_matching_histos_y(ay):
                    histo_residy[clz].Fill(dif)
                    histo_lad_res_y[ladid].Fill(dif)
                    histo_residy_vs_x[clz].Fill(dif, ax*arg + bx)
                    histo_residy_vs_y[clz].Fill(dif, ay*arg + by)
                    if positivex: histo_residy_vs_y_posv[clz].Fill(dif, ay*arg + by)
                    else:         histo_residy_vs_y_negv[clz].Fill(dif, ay*arg + by)
                    histo_chisqy[clz].Fill(dif*dif / covyy[i])
                    histo_covyy[clz].Fill( covyy[i])

        

    def __process_tracks_inspect_biased__(self, track, xclusters, yclusters, covxx, covyy,  clusterx_y, clustery_x):
        if self.__run_on_compr_data__:
            print "Running 'biased' plots with compressed data is not implemented (OBSOLETE)!"
            print "     Please consider revising the code in order not to use 'biased' plots! ==> throwing Exception!"
            raise Exception()


        if len(xclusters)==MAX_POSSIBLE_CLUSTERS and len(yclusters)==MAX_POSSIBLE_CLUSTERS:



            argsx  = [self.get_aligned_z(xclusters[i], xclusters[i].GetX(), clusterx_y[i],track)  for i in xrange(len(xclusters))]
            valsx  = [self.get_aligned_x(xclusters[i],                      clusterx_y[i],track)  for i in xrange(len(xclusters))]
            clzx   = [c.GetZ() for c in xclusters]
            ax, bx = self.get_linear_fit(argsx, valsx)
        
            argsy = [self.get_aligned_z(yclusters[i], clustery_x[i], yclusters[i].GetY(),track)  for i in xrange(len(yclusters))]
            valsy = [self.get_aligned_y(yclusters[i], clustery_x[i]                     ,track)  for i in xrange(len(yclusters))]
            clzy  = [c.GetZ() for c in yclusters]
            ay, by = self.get_linear_fit(argsy, valsy)

            self.trackchisqx.Fill(track.getChi2())
            resid = self.get_residuals(argsx, valsx, ax, bx)
            for i in xrange(len(resid)):
                arg = resid[i][0]
                dif = resid[i][1]
                clz = clzx[i]
                #histotag = "_z_%d_mm"%clz
                
                # add histograms
                for histo_residx, histo_residx_vs_x, histo_residx_vs_x_posv, histo_residx_vs_x_negv, histo_residx_vs_y, histo_chisqx, histo_covxx, tmp in self.get_matching_histos_x(ax):
                    if clz not in histo_residx: 
                        self.create_control_histos_x(clz, histo_residx, histo_residx_vs_x, histo_residx_vs_x_posv, histo_residx_vs_x_negv, histo_residx_vs_y, histo_chisqx, histo_covxx)
                    
                    
                # Fill histograms
                positivey = True if xclusters[i].getLadderHardware()/NLADDERS_TRB in [2,3] else False
                for histo_residx, histo_residx_vs_x, histo_residx_vs_x_posv, histo_residx_vs_x_negv, histo_residx_vs_y, histo_chisqx, histo_covxx, tmp in self.get_matching_histos_x(ax):
                    histo_residx[clz].Fill(dif)
                    histo_residx_vs_x[clz].Fill(dif, ax*arg + bx)
                    if positivey: histo_residx_vs_x_posv[clz].Fill(dif, ax*arg + bx)
                    else:         histo_residx_vs_x_negv[clz].Fill(dif, ax*arg + bx)
                    histo_residx_vs_y[clz].Fill(dif, ay*arg + by)
                    histo_chisqx[clz].Fill(dif * dif / covxx[i])
                    histo_covxx[clz].Fill(covxx[i]) 
                
            
 
            self.trackchisqy.Fill(track.getChi2())
            resid = self.get_residuals(argsy, valsy, ay, by)
            for i in xrange(len(resid)):
                arg = resid[i][0]
                dif = resid[i][1]
                clz = clzy[i]
                #histotag = "_z_%d_mm"%clz
                
                #add histograms
                for histo_residy, histo_residy_vs_x, histo_residy_vs_y, histo_residy_vs_y_posv, histo_residy_vs_y_negv, histo_chisqy, histo_covyy, tmp in self.get_matching_histos_y(ay):
                    if clz not in histo_residy: 
                        self.create_control_histos_y(clz, histo_residy, histo_residy_vs_x, histo_residy_vs_y, histo_residy_vs_y_posv, histo_residy_vs_y_negv, histo_chisqy, histo_covyy)
                
                positivex = True if yclusters[i].getLadderHardware()/NLADDERS_TRB in [0,1] else False
                for histo_residy, histo_residy_vs_x, histo_residy_vs_y, histo_residy_vs_y_posv, histo_residy_vs_y_negv, histo_chisqy, histo_covyy, tmp in self.get_matching_histos_y(ay):
                    histo_residy[clz].Fill(dif)
                    histo_residy_vs_x[clz].Fill(dif, ax*arg + bx)
                    histo_residy_vs_y[clz].Fill(dif, ay*arg + by)
                    if positivex: histo_residy_vs_y_posv[clz].Fill(dif, ay*arg + by)
                    else:         histo_residy_vs_y_negv[clz].Fill(dif, ay*arg + by)
                    histo_chisqy[clz].Fill(dif*dif / covyy[i])
                    histo_covyy[clz].Fill( covyy[i])

         
    def create_control_histos_x(self, clz, residx, residx_vs_x, residx_vs_x_posv, residx_vs_x_negv, residx_vs_y, chisqx, covxx):
        histotag = "_z_%d_mm"%clz
        residx[clz] = TH1F(self.RESIDUAL_HISTO_X_NAME+histotag, 
                                                                   self.RESIDUAL_HISTO_X_NAME + histotag, 
                                                                   self.RESIDUAL_HISTO_NBINS, 
                                                                   self.RESIDUAL_HISTO_MIN, 
                                                                   self.RESIDUAL_HISTO_MAX   
                                                                   )
        residx_vs_x[clz] = TH2F(self.RESIDUAL_2D_HISTO_RESX_COORX_NAME + histotag,
                                                  self.RESIDUAL_2D_HISTO_RESX_COORX_NAME + histotag,
                                                  self.RESIDUAL_2D_HISTO_RES_NBINS,
                                                  self.RESIDUAL_2D_HISTO_RES_MIN,
                                                  self.RESIDUAL_2D_HISTO_RES_MAX,
                                                  self.RESIDUAL_2D_HISTO_COOR_NBINS,
                                                  self.RESIDUAL_2D_HISTO_COOR_MIN,
                                                  self.RESIDUAL_2D_HISTO_COOR_MAX
                                                  )
        residx_vs_x_posv[clz] = TH2F(
                                                  self.RESIDUAL_2D_HISTO_RESX_POSV_COORX_NAME + histotag,
                                                  self.RESIDUAL_2D_HISTO_RESX_POSV_COORX_NAME + histotag,
                                                  self.RESIDUAL_2D_HISTO_RES_NBINS,
                                                  self.RESIDUAL_2D_HISTO_RES_MIN,
                                                  self.RESIDUAL_2D_HISTO_RES_MAX,
                                                  self.RESIDUAL_2D_HISTO_COOR_NBINS,
                                                  self.RESIDUAL_2D_HISTO_COOR_MIN,
                                                  self.RESIDUAL_2D_HISTO_COOR_MAX
                                                  )
        residx_vs_x_negv[clz] = TH2F(
                                                  self.RESIDUAL_2D_HISTO_RESX_NEGV_COORX_NAME + histotag,
                                                  self.RESIDUAL_2D_HISTO_RESX_NEGV_COORX_NAME + histotag,
                                                  self.RESIDUAL_2D_HISTO_RES_NBINS,
                                                  self.RESIDUAL_2D_HISTO_RES_MIN,
                                                  self.RESIDUAL_2D_HISTO_RES_MAX,
                                                  self.RESIDUAL_2D_HISTO_COOR_NBINS,
                                                  self.RESIDUAL_2D_HISTO_COOR_MIN,
                                                  self.RESIDUAL_2D_HISTO_COOR_MAX
                                                  )
        residx_vs_y[clz] = TH2F(self.RESIDUAL_2D_HISTO_RESX_COORY_NAME + histotag,
                                                  self.RESIDUAL_2D_HISTO_RESX_COORY_NAME + histotag,
                                                  self.RESIDUAL_2D_HISTO_RES_NBINS,
                                                  self.RESIDUAL_2D_HISTO_RES_MIN,
                                                  self.RESIDUAL_2D_HISTO_RES_MAX,
                                                  self.RESIDUAL_2D_HISTO_COOR_NBINS,
                                                  self.RESIDUAL_2D_HISTO_COOR_MIN,
                                                  self.RESIDUAL_2D_HISTO_COOR_MAX
                                                  )
                    
        chisqx[clz]  = TH1F(self.CHISQ_HISTO_X_NAME+ histotag,
                                             self.CHISQ_HISTO_X_NAME+ histotag,
                                             self.CHISQ_HISTO_NBINS,
                                             self.CHISQ_HISTO_MIN,
                                             self.CHISQ_HISTO_MAX
                                             )  
                    
        covxx[clz]  = TH1F(self.COV_HISTO_XX_NAME+ histotag,
                                            self.COV_HISTO_XX_NAME+ histotag,
                                            self.COV_HISTO_NBINS,
                                            self.COV_HISTO_MIN,
                                            self.COV_HISTO_MAX
                                            )


    def create_control_histos_y(self, clz, residy, residy_vs_x, residy_vs_y, residy_vs_y_posv, residy_vs_y_negv, chisqy, covyy):
        histotag = "_z_%d_mm"%clz
        residy[clz] = TH1F(self.RESIDUAL_HISTO_Y_NAME + histotag, 
                                                                   self.RESIDUAL_HISTO_Y_NAME + histotag, 
                                                                   self.RESIDUAL_HISTO_NBINS, 
                                                                   self.RESIDUAL_HISTO_MIN, 
                                                                   self.RESIDUAL_HISTO_MAX   
                                                                   )
        residy_vs_x[clz] = TH2F(self.RESIDUAL_2D_HISTO_RESY_COORX_NAME + histotag,
                                                  self.RESIDUAL_2D_HISTO_RESY_COORX_NAME + histotag,
                                                  self.RESIDUAL_2D_HISTO_RES_NBINS,
                                                  self.RESIDUAL_2D_HISTO_RES_MIN,
                                                  self.RESIDUAL_2D_HISTO_RES_MAX,
                                                  self.RESIDUAL_2D_HISTO_COOR_NBINS,
                                                  self.RESIDUAL_2D_HISTO_COOR_MIN,
                                                  self.RESIDUAL_2D_HISTO_COOR_MAX
                                                  )
        residy_vs_y[clz] = TH2F(self.RESIDUAL_2D_HISTO_RESY_COORY_NAME + histotag,
                                                  self.RESIDUAL_2D_HISTO_RESY_COORY_NAME + histotag,
                                                  self.RESIDUAL_2D_HISTO_RES_NBINS,
                                                  self.RESIDUAL_2D_HISTO_RES_MIN,
                                                  self.RESIDUAL_2D_HISTO_RES_MAX,
                                                  self.RESIDUAL_2D_HISTO_COOR_NBINS,
                                                  self.RESIDUAL_2D_HISTO_COOR_MIN,
                                                  self.RESIDUAL_2D_HISTO_COOR_MAX
                                                  )
        residy_vs_y_posv[clz] = TH2F(
                                                  self.RESIDUAL_2D_HISTO_RESY_POSV_COORY_NAME + histotag,
                                                  self.RESIDUAL_2D_HISTO_RESY_POSV_COORY_NAME + histotag,
                                                  self.RESIDUAL_2D_HISTO_RES_NBINS,
                                                  self.RESIDUAL_2D_HISTO_RES_MIN,
                                                  self.RESIDUAL_2D_HISTO_RES_MAX,
                                                  self.RESIDUAL_2D_HISTO_COOR_NBINS,
                                                  self.RESIDUAL_2D_HISTO_COOR_MIN,
                                                  self.RESIDUAL_2D_HISTO_COOR_MAX
                                                  )
        residy_vs_y_negv[clz] = TH2F(
                                                  self.RESIDUAL_2D_HISTO_RESY_NEGV_COORY_NAME + histotag,
                                                  self.RESIDUAL_2D_HISTO_RESY_NEGV_COORY_NAME + histotag,
                                                  self.RESIDUAL_2D_HISTO_RES_NBINS,
                                                  self.RESIDUAL_2D_HISTO_RES_MIN,
                                                  self.RESIDUAL_2D_HISTO_RES_MAX,
                                                  self.RESIDUAL_2D_HISTO_COOR_NBINS,
                                                  self.RESIDUAL_2D_HISTO_COOR_MIN,
                                                  self.RESIDUAL_2D_HISTO_COOR_MAX
                                                  )
        chisqy[clz]  = TH1F(self.CHISQ_HISTO_Y_NAME+ histotag,
                                             self.CHISQ_HISTO_Y_NAME+ histotag,
                                             self.CHISQ_HISTO_NBINS,
                                             self.CHISQ_HISTO_MIN,
                                             self.CHISQ_HISTO_MAX
                                             )
                    
        covyy[clz]  = TH1F(self.COV_HISTO_YY_NAME+ histotag,
                                            self.COV_HISTO_YY_NAME+ histotag,
                                            self.COV_HISTO_NBINS,
                                            self.COV_HISTO_MIN,
                                            self.COV_HISTO_MAX
                                            )




    def create_control_histos_ladders_y(self, resid_lad_y, lad):
        histoname = "residue_lady_%d"%lad
        resid_lad_y[lad] = TH1F(histoname,
	                        histoname,
                                self.RESIDUAL_HISTO_NBINS, 
                                self.RESIDUAL_HISTO_MIN, 
                                self.RESIDUAL_HISTO_MAX   
                                )

    def create_control_histos_ladders_x(self, resid_lad_x, lad):
        histoname = "residue_ladx_%d"%lad
        resid_lad_x[lad] = TH1F(histoname,
	                        histoname,
                                self.RESIDUAL_HISTO_NBINS, 
                                self.RESIDUAL_HISTO_MIN, 
                                self.RESIDUAL_HISTO_MAX   
                                )

	
            
    def get_residuals(self, args, vals , a, b):
        result = []
        for i in xrange(len(args)):
            arg = args[i]
            val = vals[i]
            projval = a*arg + b
            diff    = val - projval
            result.append([arg, diff])
        return result
    
    
    def do_straight_tracks(self):
        arg = [arg for arg in sys.argv if self.VERTICAL_TRACKS_ARG in arg]
        if not arg: return
        #value = arg.split("=")[-1]
        #if value.lower() == "true":
        self.track_tan_cut_min =  self.TRACK_TAN_VERTCUT_MIN
        self.track_tan_cut_max =  self.TRACK_TAN_VERTCUT_MAX
        self.__print_info__("Using the cut on tracks inclanation > "+str(self.track_tan_cut_min))
        self.__print_info__("Using the cut on tracks inclanation < "+str(self.track_tan_cut_max))
        #elif value.lower() == "false":
        #    pass
        #else:
        #    self.__print_error__("Unknown value for the argument " + arg + " ==> exiting!")
        #    raise SystemExit
        sys.argv.remove(arg[0])

    def get_unbiased(self):
        arg = [arg for arg in sys.argv if self.UNBIASED_ARG in arg]
        if not arg: return
        self.__unbiased__ = True
        sys.argv.remove(arg[0])


    def get_forbid_more_than_one(self):
        arg = [arg for arg in sys.argv if self.FORBID_MORETHANONE_ARG in arg]
        if not arg: return
        self.__forbid_more_than_one__ = True
        sys.argv.remove(arg[0])


    def get_is_running_on_aligned(self):
        arg = [arg for arg in sys.argv if self.ISRUNNINGONALIGNED_ARG in arg]
        if not arg: return
        self.__running_on_aligned__ = True
        sys.argv.remove(arg[0])

    def get_is_dump_data(self):
        arg = [arg for arg in sys.argv if self.DUMP_DATA_ARG in arg]
        if not arg: return
        self.__dump_data__ = True
        sys.argv.remove(arg[0])


    def get_run_on_copr_data(self):
        arg = [arg for arg in sys.argv if self.RUN_ON_COMPR_DATA_ARG in arg]
        if not arg: return
        self.__run_on_compr_data__ = True
        sys.argv.remove(arg[0])

    def get_no_evt_selection(self):
        arg = [arg for arg in sys.argv if self.NO_EVT_SELECTION_ARG in arg]
        if not arg: return
        self.__no_evt_selection__ = True
        sys.argv.remove(arg[0])
        
    def get_no_inspect(self):
        arg = [arg for arg in sys.argv if self.NO_INSPECT_ARG in arg]
        if not arg: return
        self.__no_inspect__ = True
        sys.argv.remove(arg[0])

    def get_no_alignment(self):
        arg = [arg for arg in sys.argv if self.NO_ALIGNMENT_ARG in arg]
        if not arg: return
        self.__no_alignment__ = True
        sys.argv.remove(arg[0])



    def do_all_sensors(self):
        arg = [arg for arg in sys.argv if self.ALIGNSENSORS_ARG in arg]
        if not arg: return
        self.__print_info__("Aligning all sensors!")
        self.__all_sensors__ = True
        self.N_STK_LADDERS = self.N_STK_LADDERS * NSENSORS
        sys.argv.remove(arg[0])
         
        
    def di_version_tag(self):
        arg = [arg for arg in sys.argv if self.NAME_VERSION in arg]
        if not arg: return
        self.versiontag = arg[0].split("=")[-1]
        self.__print_info__("Using version tag: "+self.versiontag)
        self.versiontag  = "_"+self.versiontag
        sys.argv.remove(arg[0]) 
            
        
    def do_exclude_cracks(self):
        arg = [arg for arg in sys.argv if self.SELECTION_NO_CRACKS in arg]
        if not arg: return
        self.exclude_cracks = True
        self.__print_info__("Exclude crack regions: "+str(self.exclude_cracks))
        sys.argv.remove(arg[0]) 
        
    def get_initial_alignment(self):
        arg = [arg for arg in sys.argv if self.ALIGNMENT_FILE in arg]
        if not arg: return
        fname = arg[0].split("=")[-1]
        self.read_alignment_from_file(fname)
        sys.argv.remove(arg[0]) 
    
    def get_skip_events(self):
        arg = [arg for arg in sys.argv if self.SKIP_EVENTS_ARG in arg]
        if not arg: return
        self.skipevents = int(arg[0].split("=")[-1])
        sys.argv.remove(arg[0])
        
    def get_max_events(self):
        arg = [arg for arg in sys.argv if self.MAX_EVENTS_ARG in arg]
        if not arg: return
        self.MAX_ENTRIES = int(arg[0].split("=")[-1])
        sys.argv.remove(arg[0])
        
    def get_n_tracks_file(self):
        arg = [arg for arg in sys.argv if self.NTRACKS_FILE in arg]
        if not arg: return
        normalize_n_tracks_filename = arg[0].split("=")[-1]
        try:
            normalize_n_tracks_file = open(normalize_n_tracks_filename, "rb")
            self.normalize_n_tracks = cPickle.load(normalize_n_tracks_file)
            self.normalize_n_tracks = cPickle.load(normalize_n_tracks_file)
            normalize_n_tracks_file.close()
            self.__is_ntracks_file__ = True
        except:
            print "WARNING: n-tracks file not found ",normalize_n_tracks_filename
            print "Not doing normalization!"
        sys.argv.remove(arg[0])
        
    def get_track_chi_sq_cut(self):
        arg = [arg for arg in sys.argv if self.TRACKCHISQMAX_ARG in arg]
        if not arg: return
        self.chi_sq_cut = float(arg[0].split("=")[-1])
        sys.argv.remove(arg[0])

    def get_aligntrack_chi_sq_cut(self):
        arg = [arg for arg in sys.argv if self.TRACKCHISQMAX_AL_ARG in arg]
        if not arg: return
        self.aligntrack_chi_sq_cut = float(arg[0].split("=")[-1])
        sys.argv.remove(arg[0])
         

    def get_cut_unbiased_track_residuals(self):
        arg = [arg for arg in sys.argv if self.CUT_UNBIASED_RESID in arg]
        if not arg: return
        self.__unbiased_track_residuals_max__ = float(arg[0].split("=")[-1])
        sys.argv.remove(arg[0])

    def get_cut_other_projection(self):
        arg = [arg for arg in sys.argv if self.CUT_OTHERPROJECTION in arg]
        if not arg: return
        self.__unbiased_track_residuals_otherprojection_max__ = float(arg[0].split("=")[-1])
        sys.argv.remove(arg[0])
        

   
                
    def read_bad_channels(self):
        arg = [arg for arg in sys.argv if self.BAD_CHANNELS_ARG in arg]
        if not arg: return
        fname = arg[0].split("=")[-1]
        f = open(fname)
        self.__print_info__("Using bad channels file: "+fname)
        l = f.readlines()
        for line in l:
            v = line.split("\n")[0].split(",")
            if len(v)<2: continue
            ladder = int(v[0])
            chnl   = int(v[1])
            self.__print_debug__("badchannel: %3d  %3d"%(ladder,chnl))
            if ladder not in self.badchannels:
                self.badchannels[ladder] = []
            self.badchannels[ladder].append(chnl)
        f.close()
        sys.argv.remove(arg[0])
        
        
    def read_alignment_from_file(self, fname):
        f = open(fname, "r")
        lines = f.readlines()
        for line in lines:
            vals = line.strip().split("\n")[0]
            if not vals: continue
            if vals[0] == "#": continue
            vals = vals.split(",")
            if len(vals)<3: continue
            ladderid = int   (vals[0])
            dx       = float (vals[1])
            dy       = float (vals[2])
            dz       = float (vals[3])
            dthetax  = float (vals[4])
            dthetay  = float (vals[5])
            dthetaz  = float (vals[6])
            if len(vals)>7:
                alphax   = float (vals[7])
                alphay   = float (vals[8])
                betax    = float (vals[9])
                betay    = float (vals[10])
            else:
                alphax   = 0.
                alphay   = 0.
                betax    = 0.
                betay    = 0. 
               
            self.alignments[ladderid] = {"dx": dx, "dy":dy, "dz":dz, "dthetax": dthetax, "dthetay": dthetay, "dthetaz": dthetaz, "alphax":alphax, "betax":betax, "alphay":alphay, "betay":betay }
        self.__print_info__("Alignment parametes retrieved sucessfully from the file: " + fname)
        f.close()
    
    

    def get_aligned_x(self, cluster, y=None, track=None):

        if self.__running_on_aligned__:
	    for i in xrange(track.GetNPoints()):
	        if track.GetClusterX(i,self.stkclusters) != cluster: continue
	        return track.getHitMeasX(i)
	    raise Exception("Aligned x not found!")
	
	if self.__run_on_compr_data__:
            ladderid = cluster["lad"]
            x = cluster["x"]
        else:
            ladderid = self.__get_ladder_ultimate__(cluster,None,y)
            x  = cluster.GetX()

        newx = x+ self.alignments[ladderid]["dx"] * MICRON 
        if y is not None:
            newx = newx - y * self.alignments[ladderid]["dthetaz"] * MICRORAD
        return newx
    
    def get_aligned_y(self, cluster, x=None, track = None):
        if self.__running_on_aligned__:
	    for i in xrange(track.GetNPoints()):
                if track.GetClusterY(i,self.stkclusters) != cluster: continue
                return track.getHitMeasY(i)
            raise Exception("Aligned y not found!")

	if self.__run_on_compr_data__:
            ladderid = cluster["lad"]
            y = cluster["y"]
        else:
            ladderid = self.__get_ladder_ultimate__(cluster, x, None)
            y  = cluster.GetY()

        newy = y + self.alignments[ladderid]["dy"] * MICRON
        if x is not None:
            newy = newy + x * self.alignments[ladderid]["dthetaz"] * MICRORAD
        return newy 
    
    
    def get_aligned_z(self, cluster, x=None, y=None, track = None):
        if self.__running_on_aligned__:
	    for i in xrange(track.GetNPoints()):
                if track.GetClusterX(i,self.stkclusters) == cluster:
                    return track.getHitMeasX_Z(i)
                if track.GetClusterY(i,self.stkclusters) == cluster:
                    return track.getHitMeasY_Z(i)
            raise Exception("Aligned z not found!")

	if self.__run_on_compr_data__:
            ladderid = cluster["lad"]
            z = cluster["z"]
        else:
            ladderid = self.__get_ladder_ultimate__(cluster, x, y)
            z  = cluster.GetZ()

        newz = z + self.alignments[ladderid]["dz"] * MICRON
        if x is not None:
            newz = newz - x * self.alignments[ladderid]["dthetay"] * MICRORAD
            newz = newz + x**2 * self.alignments[ladderid]["alphax"] +  x**3 * self.alignments[ladderid]["betax"]
        if y is not None:
            newz = newz + y * self.alignments[ladderid]["dthetax"] * MICRORAD
            newz = newz + y**2 * self.alignments[ladderid]["alphay"] +  y**3 * self.alignments[ladderid]["betay"]

	#arg = y if cluster.isX() else x

        return newz
    
            
    def is_cluster_bad_channel(self, cluster):
        if not self.badchannels:
            return False
        ladder = cluster.getLadderHardware()
        if ladder not in self.badchannels:
            return False
        minc = cluster.getIndex1()
        maxc = cluster.getIndex1() + cluster.getNstrip() -1
        bd = self.badchannels[ladder]
        for c in bd:
            if c<minc: continue
            if c>maxc: continue
            return True
        return False 
    
    def calculate_profile(self, hiso2d):
        #h = hiso2d.ProfileX()
        result = []
        c = TCanvas(hiso2d.GetName())
        c.cd()
        hiso2d.Draw()
        step = self.RESIDUAL_2D_HISTO_COOR_NBINS / self.RESIDUAL_2D_HISTO_COORSECTION_NSTEPS
        bins = xrange(1, self.RESIDUAL_2D_HISTO_COOR_NBINS, step)
        for b in bins:
            h = hiso2d.ProfileX(hiso2d.GetName()+"_px_" + str(b) ,b, b+step -1, "")
            y = (
                 self.RESIDUAL_2D_HISTO_COOR_MAX - 
                 self.RESIDUAL_2D_HISTO_COOR_MIN 
                 ) * (b+step/2) / self.RESIDUAL_2D_HISTO_COOR_NBINS + self.RESIDUAL_2D_HISTO_COOR_MIN
            x = h.GetMean()    
            p = TMarker(x,y,21)
            p.SetMarkerSize(0.6)
            p.SetMarkerColor(kBlue)
            p.Draw("same")
            result.append(p)
        c.Write()
        return result
            
            
    
    def Finalize(self):
        # write histograms
        print 
        self.fout.cd()
        """
        plotdir = self.fout.mkdir("Plots")
        plotdir.cd()
        markers1 = [self.calculate_profile(h) for h in self.residx_vs_x.values()]
        markers2 = [self.calculate_profile(h) for h in self.residx_vs_y.values()]
        markers3 = [self.calculate_profile(h) for h in self.residy_vs_y.values()]
        markers4 = [self.calculate_profile(h) for h in self.residy_vs_x.values()]
        """

        histodir = self.fout.mkdir("Histos")
        histodir.cd()
        self.histo_distance.Write()
        self.histo_nclustersbefore.Write()
        self.histo_nclustersafter.Write()
        self.histo_angdistance.Write()
        self.histo_chargex.Write()
        self.histo_chargey.Write()
        self.histo_vertex_dz.Write()
        self.histo_vertex_xz.Write()
        self.histo_vertex_yz.Write()
        self.histo_vertex_focus_dz.Write()
        self.histo_vertex_focus_xz.Write()
        self.histo_vertex_focus_yz.Write()
        self.ngoodtracks.Write()
        self.trackchisqx.Write()
        self.trackchisqy.Write()
        folders = ["ANGLE_" + str(val) for val in ANGLES] + ["ANGLE_ALL"]
        for i in xrange(len(folders)):
            fld = histodir.mkdir(folders[i])
            fld.cd() 
            
            [h.Write() for h in self.residx[i].values()]
            [h.Write() for h in self.residy[i].values()]
            [h.Write() for h in self.residx_vs_x[i].values()]
            [h.Write() for h in self.residx_vs_x_posv[i].values()]
            [h.Write() for h in self.residx_vs_x_negv[i].values()]
            [h.Write() for h in self.residx_vs_y[i].values()]
            [h.Write() for h in self.residy_vs_y[i].values()]
            [h.Write() for h in self.residy_vs_y_posv[i].values()]
            [h.Write() for h in self.residy_vs_y_negv[i].values()]
            [h.Write() for h in self.residy_vs_x[i].values()]
            
            [h.Write() for h in self.chisqx[i].values()]
            [h.Write() for h in self.chisqy[i].values()]
            [h.Write() for h in self.covxx[i].values()]
            [h.Write() for h in self.covyy[i].values()]
            [h.Write() for h in self.covxy[i].values()]


        histodir2 = self.fout.mkdir("Histos_Ladders_x")
        histodir2.cd()
        for i in xrange(len(folders)):
            fld = histodir2.mkdir(folders[i])
            fld.cd()
            [h.Write() for h in self.resid_all_lad_x[i].values()]


        histodir3 = self.fout.mkdir("Histos_Ladders_y")
        histodir3.cd()
        for i in xrange(len(folders)):
            fld = histodir3.mkdir(folders[i])
            fld.cd()
            [h.Write() for h in self.resid_all_lad_y[i].values()]
            
        
        histodir_before = self.fout.mkdir("Histos_before_selection")
        histodir_before.cd()
        [h.Write() for h in self.nxcluster_vs_y]
        [h.Write() for h in self.nycluster_vs_x]
        #[h.Write() for h in self.centralcrack_x]
        #[h.Write() for h in self.centralcrack_y]
        self.centralcrack_x.Write()
        self.centralcrack_y.Write()
        
        histodir_der = self.fout.mkdir("Histos_derivatives")
        histodir_der.cd()
        self.histo_track_unnormalized_chisq_for_derivatives.Write()
        self.histo_alignedtrack_chisq.Write()
        self.histo_track_kalman_chisq_for_derivatives.Write()
        self.histo_track_angle_x.Write()
        self.histo_track_angle_y.Write()
        self.histo_track_angle.Write()
        self.histo_track_anglex_chisqx.Write()
        self.histo_track_angley_chisqy.Write()
        self.histo_track_angle_chisq.Write()
        if not self.__run_on_compr_data__:
            self.histo_incline_nstrips_2D_x.Write()
            self.histo_incline_nstrips_2D_y.Write()
            self.histo_incline_nstrips_2D.Write()
            self.histo_incline_charge_2D_x.Write()
            self.histo_incline_charge_2D_y.Write() 
            self.histo_incline_charge_2D.Write() 
            self.histo_incline_chargevacorrected_2D_x.Write() 
            self.histo_incline_chargevacorrected_2D_y.Write()
            self.histo_incline_chargevacorrected_2D.Write()
            [histo.Write() for histo in self.histo_incline_nstrips_x]
            [histo.Write() for histo in self.histo_incline_nstrips_y]
            [histo.Write() for histo in self.histo_incline_nstrips]
            [histo.Write() for histo in self.histo_incline_charge_x]
            [histo.Write() for histo in self.histo_incline_charge_y]
            [histo.Write() for histo in self.histo_incline_charge] 
            [histo.Write() for histo in self.histo_incline_chargevacorr_x]
            [histo.Write() for histo in self.histo_incline_chargevacorr_y]
            [histo.Write() for histo in self.histo_incline_chargevacorr]
        


        # output data file
        if self.__dump_data__:
            foutdata = open(self.OUT_DATA_FILE,"wb")
            #self.data_dump_all = np.array(self.data_dump_all)
            for key in self.data_dump_all.keys():
                self.data_dump_all[key] = np.array(self.data_dump_all[key])
            np.savez(foutdata,**self.data_dump_all)
            foutdata.close()

        #
        # dump derivatives to the file
        #
        f = open(self.OUT_DERIVATIVES_FILE, "wb")
        cPickle.dump(self.dchish_dpars,f)
        cPickle.dump(self.tracks_for_lad,f)
        f.close()

        
        
        
        #@ Get ChiSq derivatives
        #if self.dchish_dpars["ntracks"]:
        #    self.dchish_dpars["derivatives"] = np.divide(["derivatives"], self.dchish_dpars["ntracks"])
            
        
        self.fout.Close()
        outrootname = self.ROOTFILE_PREFIX + self.versiontag +".root"
        os.system("mv "+self.TMP_ROOT_FILE_NAME+" "+ outrootname)
        #self.fin.Close()

       
        
        
    def PrinChiSqDerivatives(self):


        print
        print "===========> Chi-square derivatives: " 
        print
        print "-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" 
        print " lad |          dchidx           dchidy           dchidz         dchidthx         dchidthx         dchidthz      alphax         alphay      betax           betay          |      #tracks"
        print "-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"
        for i in xrange(self.N_STK_LADDERS):
            #print "test:", self.dchish_dpars
            #print 
            print "%3d  | %15f  %15f  %15f  %15f  %15f  %15f %15f %15f %15f %15f |  %10d"%tuple([i] + list(self.dchish_dpars["derivatives"][i]) + [sum(self.tracks_for_lad[i])])
            
            # NORMALIZED - DO NOT DELETE - USEFUL FOR DEBUGGING!
            #print "%3d  | %15f  %15f  %15f  %15f  %15f  %15f |  %10d"%tuple([i] + [der for der in self.dchish_dpars["derivatives"][i][:3]] + [der/400. for der in self.dchish_dpars["derivatives"][i][3:]] + [self.tracks_for_lad[i]])
        print "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"
        print
        
    def PrintTimeStatistics(self):
        MASK  = "%20s  :  %20f  ,  %10d"
        print
        print "Time consumption statistsics: "
        #self.__time_elapsed__["Total"] =    
        for key, value in self.__time_elapsed__.items():
            #if "test" in key.lower().strip() : continue
            if key.lower().strip() == "total": continue
            print MASK%(key, value, self.__number_of_calls__[key])
        print "----------------------------------------------------------------"
        key = "TOTAL"
        print MASK%("Total",self.__time_elapsed__[key], self.__number_of_calls__[key])
        print "----------------------------------------------------------------"
        
        
        
            
        
                
                
            
a = Alignment()
#a.MAX_ENTRIES = 100



one = datetime.datetime.now()        
a.Init()
a.Run()
a.Finalize()
two = datetime.datetime.now()        



#a.PrintTimeStatistics()
a.PrinChiSqDerivatives()

print
print two - one
print 
a.PrintTimeStatistics()
print






