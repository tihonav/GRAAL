'''
Created on Jul 7, 2015

@author: andrii
'''

import sys
import os

#sys.path.append(os.environ["DMPSWSYS"]+"/python/COSM_APR2015/")
#from Align import Alignment

from glob import glob
import commands
import time
import ROOT
from datetime import datetime
import cPickle
import numpy as np


#DETECTOR_PREREQUISITE_SCRIPT = "/afs/cern.ch/work/a/andrii/public/detector/sw/setup_detector_prerequisits.sh"
DETECTOR_PREREQUISITE_SCRIPT = "/home/software/detector/setup/setup-externals.sh"
#DAMPME_INSTALL_PATH       = "/afs/cern.ch/work/a/andrii/public/detector/trunk/Install"
#DAMPME_INSTALL_PATH       = "/detector/data2/BTNOV2015/software/trunk/Install"
#DAMPME_INSTALL_PATH       = "/afs/cern.ch/work/a/andrii/public/detector/trunk_rep1_merged/Install"
DAMPME_INSTALL_PATH       = os.environ["DMPSWSYS"]
DETECTOR_LAUNCHER            = "$DMPSWSYS/script/BT2014/DetectorJobLauncher.sh"
#SCRIPT                    = "python/COSM_APR2015/Align10.py %s --bdchnl=%s/share/Calibration/ParametersTRACKER_FM/bad_chan.txt --align=%s" # --align=$DMPSWSYS/../Calibration/ParametersTRACKER_FM/tracker_alignment_step0.txt"
SCRIPT                    = "python/COSM_APR2015/Align.py %s --bdchnl=%s/share/Calibration/ParametersTRACKER_FM/bad_chan.txt --align=%s" # --align=$DMPSWSYS/../Calibration/ParametersTRACKER_FM/tracker_alignment_step0.txt"
#FILE_PATH                 = "/detector/data1/cosmics_dpnc_april_2015/REC2"
#FILE_PATH                 = "/detector/data2/FMdata/Satellite_Aging_Test/processing1/REC"
#FILE_PATH                 = "/detector/data2/FMdata/Satellite_Aging_Test/processing1/part*/REC/"
#FILE_PATH                 = "/detector/data2/FMdata/Satellite_Inflight_Simulation_2015/processing1/REC"
#FILE_PATH                 = "/detector/data2/FMdata/Satellite_Inflight_Simulation_2015/processing_testaligntracker/REC"
#FILE_PATH                 = "/detector/data1/shanghai_tests_may_2015/REC"
#FILE_PATH                 = "/detector/data2/FMdata/Shanghai_Satellite_Integration_May2015/processing2/REC/"
#FILE_PATH                 = "/detector/data2/FMdata/Satellite_Inflight_Simulation_2015/processing_testaligntracker/REC/"
#FILE_PATH                 = "/detector/data2/FMdata/Satellite_Inflight_Simulation_2015/processing1/REC/" 
"""
FILE_PATH                 = ["/detector/data2/FMdata/Orbit_Data/processing3/REC/nrm_201601",
                             "/detector/data2/FMdata/Orbit_Data/processing3/REC/nrm_20160201",
                             "/detector/data2/FMdata/Orbit_Data/processing3/REC/nrm_20160202",
                             "/detector/data2/FMdata/Orbit_Data/processing3/REC/nrm_20160203",
                             "/detector/data2/FMdata/Orbit_Data/processing3/REC/nrm_20160204",
                             ]
"""
#OUTPUT_PATH               = "/detector/data1/andrii/alignment"
#ALIGNMENT_TMP             = "/detector/data1/andrii/alignment/alignment_tmp.txt"
"""
FILE_PATH                 = [
                             "/detector/data2/FMdata/Orbit_Data/processing3/REC/nrm_20160118",
                             "/detector/data2/FMdata/Orbit_Data/processing3/REC/nrm_20160119",
                             "/detector/data2/FMdata/Orbit_Data/processing3/REC/nrm_20160120",
                             "/detector/data2/FMdata/Orbit_Data/processing3/REC/nrm_20160121",
                             "/detector/data2/FMdata/Orbit_Data/processing3/REC/nrm_20160122",
                            ]
"""
#FILE_PATH                 = ["/detector/data2/FMdata/Orbit_Data/processing3/REC/nrm_20160204003000-20160204003000__20160204B000000S031-20160204B003115S894_ext_REC"]
#FILE_PATH                 = ["/detector/data2/datageneva/FMdata/Orbit_Data/processing3/REC/nrm_20160120*root"]
#FILE_PATH                 = ["/detector/data2/datageneva/FMdata/Orbit_Data/processing3/REC/nrm_2016*root"]


"""
#@Run on compresses data - align
#----------------------------------------------------------------------------------------------
FILE_PATH                 = ["/detector/data1/andrii/alignment_input_data_25_02_2016_v2/part*npy"]
OUTPUT_PATH               = "/detector/data1/andrii/alignment2"
ALIGNMENT_TMP             = "/detector/data1/andrii/alignment2/alignment_tmp.txt"
DUMP_DATA                 = False
VERBOSITY                 = False
RUN_ON_COMPR_DATA         = True
WAIT                      = 14400  # run      (sec)
ALGSTEPFIRST              = 3.125  # test2
MAXALGITERATIONS          = 400    # test2
MAXEVENTSPERJOB           = 1000000
RERUN_FIRST               = True
NO_INSPECT                = True
#----------------------------------------------------------------------------------------------
"""

"""
#@Run on compresses data - ispect
#----------------------------------------------------------------------------------------------
FILE_PATH                 = ["/detector/data1/andrii/alignment_input_data_25_02_2016_v2/part*npy"]
OUTPUT_PATH               = "/detector/data1/andrii/alignment2"
ALIGNMENT_TMP             = "/detector/data1/andrii/alignment2/alignment_tmp.txt"
DUMP_DATA                 = False
VERBOSITY                 = False
RUN_ON_COMPR_DATA         = True
WAIT                      = 14400  # run      (sec)
ALGSTEPFIRST              = 0.0    # test2
MAXALGITERATIONS          = 1      # test2
MAXEVENTSPERJOB           = 1000000
RERUN_FIRST               = False
NO_INSPECT                = False
#----------------------------------------------------------------------------------------------
"""


#@Run on compresses data - ispect
#----------------------------------------------------------------------------------------------
"""
FILE_PATH                 = [
                             "/detector/data1/andrii/alignment_input_data_05_06_2016/part_out_000*",
                             "/detector/data1/andrii/alignment_input_data_05_06_2016/part_out_001*",
                             "/detector/data1/andrii/alignment_input_data_05_06_2016/part_out_002*",
                             "/detector/data1/andrii/alignment_input_data_05_06_2016/part_out_003*",
                             "/detector/data1/andrii/alignment_input_data_05_06_2016/part_out_004*",
                             "/detector/data1/andrii/alignment_input_data_05_06_2016/part_out_005*",
                             "/detector/data1/andrii/alignment_input_data_05_06_2016/part_out_006*",
                            ]
"""
"""
FILE_PATH                 = [
				"/detector/data1/andrii/alignment_input_data_firsthalf2016/part_out_000--20160101005224_20160101012153_0-20160101060642_20160101070947_0.npy",
				"/detector/data1/andrii/alignment_input_data_firsthalf2016/part_out_001--20160101071048_20160101074017_0-20160101122506_20160101132811_0.npy",
				"/detector/data1/andrii/alignment_input_data_firsthalf2016/part_out_002--20160101132912_20160101135841_0-20160101170854_20160101181159_0.npy",
				"/detector/data1/andrii/alignment_input_data_firsthalf2016/part_out_003--20160101181300_20160101184229_0-20160101232718_20160102003023_0.npy",
				"/detector/data1/andrii/alignment_input_data_firsthalf2016/part_out_004--20160102003124_20160102010052_0-20160102051512_20160102054441_0.npy",
				"/detector/data1/andrii/alignment_input_data_firsthalf2016/part_out_005--20160102054542_20160102064847_0-20160102102930_20160102113235_0.npy",
				"/detector/data1/andrii/alignment_input_data_firsthalf2016/part_out_006--20160102113336_20160102120305_0-20160102164754_20160102175059_0.npy",
				"/detector/data1/andrii/alignment_input_data_firsthalf2016/part_out_007--20160102175200_20160102182129_0-20160102223548_20160102230517_0.npy",
                            ]
"""
FILE_PATH                 = [
"part_out_096--20160812T052804-20160812T104750.npy",
#"part_out_097--20160812T114628-20160812T163016.npy",
"part_out_098--20160812T170614-20160812T215002.npy",
#"part_out_099--20160812T224840-20160813T095052.npy",
"part_out_100--20160813T112528-20160814T110428.npy",
#"part_out_101--20160814T114026-20160814T162414.npy",
"part_out_102--20160814T172252-20160814T224238.npy",
#"part_out_000--20160814T234116-20160815T050102.npy",
"part_out_001--20160815T055940-20160815T104328.npy",
#"part_out_002--20160815T111926-20160815T160314.npy",
"part_out_003--20160815T170152-20160815T222138.npy",
#"part_out_004--20160815T232016-20160816T040404.npy",
"part_out_005--20160816T044002-20160816T092350.npy",
#"part_out_006--20160816T102228-20160816T154214.npy",

]
FILE_PATH                 = ["/detector/data1/andrii/alignment_input_data_firsthalf2016/" + f for f in FILE_PATH]
OUTPUT_PATH               = "/detector/data1/andrii/alignment20_junecrosscheck/FULLSTEP_tmp_16_____INCLINATIONNORMALIZATION"
ALIGNMENT_TMP             = "/detector/data1/andrii/alignment20_junecrosscheck/FULLSTEP_tmp_16_____INCLINATIONNORMALIZATION/alignment_tmp.txt"
DUMP_DATA                 = False
VERBOSITY                 = True
RUN_ON_COMPR_DATA         = True
WAIT                      = 21600  # 6 hours # 14400  # run      (sec)
ALGSTEPFIRST              = 3.90625 #500.   #6.25   #100.0  #200.0  #1.0   # test2
MAXALGITERATIONS          = 400    #200    # 400      # test2
MAXEVENTSPERJOB           = 10000  #20000 #1000000
RERUN_FIRST               = False
NO_INSPECT                = False
#----------------------------------------------------------------------------------------------

ALIGNMENT_TMP_DEFAULT     = "/afs/cern.ch/work/a/andrii/public/detector/trunk/Calibration/ParametersTRACKER_FM/tracker_noalignment.txt"           
MERGED_OUT_DUMP           = OUTPUT_PATH + "/out.dump"
TRACK_CHI_SQ_CUT          = 6.  # 2.  #6. 
TRACK_UNBIAS_RES_CUT      = 0.3
#TRACK_ALIGN_CHI_SQ_CUT    = 0.7



def run_alignment(f, skipevents = None, maxevents = None, version = None, alltracksforladder = None, trackchisqcut = None, verbosity = VERBOSITY, alignsensors=False, rerun=False):

    TMP_INPUT_ROOT = "input_tmp_REC.root"

    if verbosity: print "Processong run", f

    os.environ["DETECTOR_PREREQUISITE_SCRIPT"] = DETECTOR_PREREQUISITE_SCRIPT
    os.environ["DAMPME_INSTALL_PATH"]       = DAMPME_INSTALL_PATH
    #os.environ["DETECTORSCRIPT"]               = SCRIPT%(f,DAMPME_INSTALL_PATH, ALIGNMENT_TMP) + (" --alignsensors" if alignsensors else "")
    os.environ["DETECTORSCRIPT"]               = SCRIPT%(TMP_INPUT_ROOT,DAMPME_INSTALL_PATH, ALIGNMENT_TMP) #+ (" --alignsensors" if alignsensors else "")
    os.environ["DMPOUTPATH"]                = OUTPUT_PATH 
    os.environ["INPUT_FILE_TO_COPY"]        = f
    os.environ["INPUT_FILE_TMP"]            =  TMP_INPUT_ROOT

    #if verbosity: print "DETECTOR scrips: ", os.environ["DETECTORSCRIPT"]
    if alignsensors:
        os.environ["DETECTORSCRIPT"]+= " --alignsensors" 

    if skipevents is not None:
        os.environ["DETECTORSCRIPT"]+= " --skipevents=%d"%skipevents

    if maxevents is not None:
        os.environ["DETECTORSCRIPT"]+= " --maxevents=%d"%maxevents

    #if version is not None:
    #    os.environ["DETECTORSCRIPT"]+= " --ver=%d"%version

    if version is not None:
        os.environ["DETECTORSCRIPT"]+= " --ver=%s_%d"%(f.split("/")[-1].split(".")[0],version)

    if alltracksforladder is not None:
        os.environ["DETECTORSCRIPT"]+= " --ntracksfile="+alltracksforladder

    if trackchisqcut is not None:
        os.environ["DETECTORSCRIPT"]+= " --trackchisqmax=%f"%trackchisqcut

    if DUMP_DATA:
        os.environ["DETECTORSCRIPT"]+= " --dumpdata"

    if RUN_ON_COMPR_DATA:
        os.environ["DETECTORSCRIPT"]+= " --runoncomprdata"

    if NO_INSPECT:
        os.environ["DETECTORSCRIPT"]+= " --noinspect"

    #os.environ["DETECTORSCRIPT"]+= " --aligntrackchisqmax=%f"%TRACK_ALIGN_CHI_SQ_CUT
    os.environ["DETECTORSCRIPT"]+= " --unbiased"
    os.environ["DETECTORSCRIPT"]+= " --cutunbiasres=%f"%TRACK_UNBIAS_RES_CUT
    #os.environ["DETECTORSCRIPT"]+= " --verttracks"
    os.environ["DETECTORSCRIPT"]+= " --notmorethanone"
    #os.environ["DETECTORSCRIPT"]+= " --runningonaligned"



    # check if file exsists
    if DUMP_DATA:
        test = glob(OUTPUT_PATH+"/*"+f.split("/")[-1].split(".root")[0]+"*")
        test = [i for i in test if ".np" in i]
        if version is not None: test = [i for i in test if "_"+  str(version)+"." in i]
        if test:
            print " ... dataset ",f, " is already processed ==> skipping it!"
            return 

    # check if derivatives file exist (only used in case of reruning the jobs)
    if rerun:
        test = glob(OUTPUT_PATH+"/*"+f.split("/")[-1].split(".")[0]+"*")
        test = [i for i in test if ".dump" in i]
        if version is not None: test = [i for i in test if "_"+str(version)+"." in i]
        if test:
            print " ... derivatives for ",f, " already exist ==> skipping it!"
            return None
    
    
    if verbosity: print "DETECTOR scrips: ", os.environ["DETECTORSCRIPT"]
    

    #cmd="qsub -q veryshort -v DETECTOR_PREREQUISITE_SCRIPT,DAMPME_INSTALL_PATH,DETECTORSCRIPT,DMPOUTPATH -l mem=6000mb -l vmem=6000mb " + DETECTOR_LAUNCHER
    #cmd="qsub -q short -v DETECTOR_PREREQUISITE_SCRIPT,DAMPME_INSTALL_PATH,DETECTORSCRIPT,DMPOUTPATH -l mem=6000mb -l vmem=6000mb " + DETECTOR_LAUNCHER
    cmd="qsub -q veryshort -v DETECTOR_PREREQUISITE_SCRIPT,DAMPME_INSTALL_PATH,DETECTORSCRIPT,DMPOUTPATH,INPUT_FILE_TO_COPY,INPUT_FILE_TMP -l mem=6000mb -l vmem=6000mb " + DETECTOR_LAUNCHER
    print cmd
    raise SystemExit
    #cmd="qsub -q long -v DETECTOR_PREREQUISITE_SCRIPT,DAMPME_INSTALL_PATH,DETECTORSCRIPT,DMPOUTPATH,INPUT_FILE_TO_COPY,INPUT_FILE_TMP -l mem=6000mb -l vmem=6000mb " + DETECTOR_LAUNCHER
    #if not verbosity:  cmd+=" > tmp.log"
    ou = commands.getoutput(cmd)
    if verbosity: print "Submitted job: ", ou
    
    jobid = ou.split(".")[0]
    if jobid.isdigit():
        return jobid
    else:
        return None
    

#for f in glob(DATA_PATH+"*frd"):
#    run_calibration(f)


def checkdir(name):
    if os.path.isdir(name):
        print "Directory exists: ", name
        return
    
    result =  commands.getstatusoutput("mkdir "+name)
    if not result[0]:
        print "Directory created: ", name
        return
    
    print result
    tmp = name[:name.rfind("/")]
    result =  commands.getstatusoutput("mkdir " + tmp)
    print result
    if result[0]:
        print "Failed to create directory: ", tmp
        return
    
    result =  commands.getstatusoutput("mkdir " + name)
    print result
    if result[0]:
        print "Failed to create directory: ", name
        return
    if not os.path.isdir(name):  
        print "Still did not create the directory!"
        return 

    print "Directory created: ", name
    
def checkaligntmp():
    if not os.path.isfile(ALIGNMENT_TMP):
        print "Alignment file does not exist ==> copying the default one..."
        command = "cp %s %s"%(ALIGNMENT_TMP_DEFAULT, ALIGNMENT_TMP)
        print "Executing command: ", command
        print commands.getstatusoutput(command)
    if not os.path.isfile(ALIGNMENT_TMP):
        print "ERROR! no alignment file created yet!"

def run(alltracksforladder = None, trackchisqcut=None, exitifstuck = False, alignsensors=False, rerun = False):


    #print alignsensors
    #raise SystemExit

    DERIVASTIVESMASK = "derivatives"
    ROOTMASK         = "out_ALIGN"
    EVENTSPERJOB     = MAXEVENTSPERJOB #60000 #400000 #60000 #10000    

    checkdir(OUTPUT_PATH)
    checkaligntmp()
    #files = glob(FILE_PATH + "/*root")
    #files = glob(FILE_PATH + "*root")

    files = reduce(lambda x,y: x+y, [glob(PATTERN) for PATTERN in FILE_PATH])
    #files = [f for f in files if "0009" not in f and "0010" not in f]


    print 
    print "-------------------- Cleaning ---------------------------------"
    print
    if not rerun:
        for mask in [DERIVASTIVESMASK, ROOTMASK]:
            commands.getstatusoutput("rm "+OUTPUT_PATH+"/*" + mask +"*")


    print 
    print "-------------------- Submitting jobs --------------------------"
    print   
    nparts = 0 
    jobs = []
    for f in files:
        if ".root" in f:
            rootf = ROOT.TFile(f)
            roott = rootf.Get("CollectionTree")
            try:
                nentries = roott.GetEntries()
                #print nentries
            except:
                print
                print "Failed to open ROOT file: ",f
                print
                continue
            rootf.Close()
        elif ".np" in f:
            a=np.load(f)
            nentries = len(a[a.keys()[0]])
        else:
            raise Exception("Unknown input file type")
            

        ver = 0
        for start in xrange(0, nentries, EVENTSPERJOB):
            jobid = run_alignment(f, skipevents = start, maxevents = EVENTSPERJOB, version=ver, alltracksforladder = alltracksforladder, trackchisqcut = trackchisqcut, alignsensors=alignsensors, rerun = rerun)
            print "jobid=",jobid

            jobs.append(jobid)
            ver+=1
            nparts+=1
            sys.stdout.write("\r Jobs submitted: %d"%nparts)
            sys.stdout.flush()         
            #if nparts > 2: break
        #if nparts > 2: break
        
    
    print 
    print "-------------------- Chacking job status-----------------------"
    start = datetime.now()
    print " Time: ", start
    while True:
        #time.sleep(30)
        time.sleep(1)
        processed = glob(OUTPUT_PATH+"/*" + DERIVASTIVESMASK + "*")
        percentage = 100. * len(processed) / nparts
        elapsed = datetime.now() - start
        #if (len(processed) == nparts) or (elapsed.seconds>600 and percentage>99.5):
        if len(processed) == nparts:
            print "\nAll jobs processed!"
            break
        sys.stdout.write("\r Processed:  %10d / %10d ( %2d%% )"%(len(processed),nparts, percentage))
        sys.stdout.flush()

        if exitifstuck and elapsed.seconds > WAIT: #360000: #10800:  #5400: #3600:
            print "Something has stucked ==> rerun the job!"
            print " ... killing the jobs"
            for job in jobs:
                if job is None: continue
                qdelcmd = "qdel "+job
                print " ... ", qdelcmd
                os.system(qdelcmd)
            return None, None
       
        
        #CUT = 95
        #if percentage>CUT:
        #    print "Minimum number of jobs done %d %%"%CUT
        #    break
        
        # removing job logs
    
    commands.getoutput("rm DetectorJobLauncher.sh.*")        

    print " Time: ", datetime.now()  
   

    print 
    print "-------------------- Collecting results-----------------------"
    
    all_dchish_dpars = None
    all_tracks_for_lad = None
    for f in processed:
        fin = open(f, "rb")
        dchish_dpars = cPickle.load(fin)
        tracks_for_lad = cPickle.load(fin)
        fin.close()
        
        if all_dchish_dpars is None:
            all_dchish_dpars = dchish_dpars
            all_tracks_for_lad = tracks_for_lad
            continue

        # append
        all_tracks_for_lad = np.add(all_tracks_for_lad, tracks_for_lad)
        all_dchish_dpars["derivatives"] = np.add(all_dchish_dpars["derivatives"], dchish_dpars["derivatives"])
        all_dchish_dpars["ntracks"] += dchish_dpars["ntracks"]
        all_dchish_dpars["chisq"]   += dchish_dpars["chisq"] 
   
    # put into output file
    fout = open(MERGED_OUT_DUMP,"wb")
    cPickle.dump(all_dchish_dpars, fout)
    cPickle.dump(all_tracks_for_lad, fout)
    fout.close()
    return all_dchish_dpars, all_tracks_for_lad






NSENSORS = 4
NLADDERS = 192
NLADDERSTRB = 24
NALIGNPARSLAD = 10
ALIGNMENTFILE = "align_tmp.txt"
def dchisq_dposnegx(dchish_dpars):  
    result = 0 
    for ladder in xrange(NLADDERS):
        trbid = ladder/ NLADDERSTRB
        #ladderid = ladder % NLADDERSTRB
        if trbid not in [2,3,6,7]:
            continue
        dchidx = dchish_dpars["derivatives"][ladder][0]
        if trbid in [2,3]:
            result+= dchidx
        elif trbid in [6,7]:
            result-= dchidx
    return result


def dchisq_dposnegy(dchish_dpars):  
    result = 0 
    for ladder in xrange(NLADDERS):
        trbid = ladder/ NLADDERSTRB
        #ladderid = ladder % NLADDERSTRB
        if trbid not in [0,1,4,5]:
            continue
        dchidy = dchish_dpars["derivatives"][ladder][1]
        if trbid in [0,1]:
            result+= dchidy
        elif trbid in [4,5]:
            result-= dchidy
    return result


def dchisq_dangle(dchish_dpars):
    result = {
               "l0": { "x" : [0., 0., 0.], "y" : [0., 0., 0.] },
               "l1": { "x" : [0., 0., 0.], "y" : [0., 0., 0.] },
               "l2": { "x" : [0., 0., 0.], "y" : [0., 0., 0.] },
               "l3": { "x" : [0., 0., 0.], "y" : [0., 0., 0.] },
               "l4": { "x" : [0., 0., 0.], "y" : [0., 0., 0.] },
               "l5": { "x" : [0., 0., 0.], "y" : [0., 0., 0.] },
             }
    for ladder in xrange(NLADDERS):
        trbid = ladder/ NLADDERSTRB
        ladderid = ladder % NLADDERSTRB
        layer= 5-ladderid/4
        if trbid in [0,1,4,5]:
            coor = "y"
        else:
            coor = "x"
        for i in xrange(3):
            result["l%d"%layer][coor][i]+= dchish_dpars["derivatives"][ladder][3+i] 
    return result


def dchisq_layers(dchish_dpars):
    result = {
               "l0": { "x" : [0., 0., 0., 0., 0., 0.], "y" : [0., 0., 0., 0., 0., 0.] },
               "l1": { "x" : [0., 0., 0., 0., 0., 0.], "y" : [0., 0., 0., 0., 0., 0.] },
               "l2": { "x" : [0., 0., 0., 0., 0., 0.], "y" : [0., 0., 0., 0., 0., 0.] },
               "l3": { "x" : [0., 0., 0., 0., 0., 0.], "y" : [0., 0., 0., 0., 0., 0.] },
               "l4": { "x" : [0., 0., 0., 0., 0., 0.], "y" : [0., 0., 0., 0., 0., 0.] },
               "l5": { "x" : [0., 0., 0., 0., 0., 0.], "y" : [0., 0., 0., 0., 0., 0.] },
             }
    for sensor in xrange(NLADDERS * NSENSORS):
        ladder = sensor / 4
        trbid = ladder/ NLADDERSTRB
        ladderid = ladder % NLADDERSTRB
        layer= 5-ladderid/4
        if trbid in [0,1,4,5]:
            coor = "y"
        else:
            coor = "x"
        for i in xrange(6):
            result["l%d"%layer][coor][i]+= dchish_dpars["derivatives"][sensor][i] 
    return result


def add_angle_to_alignment(alignment, angledrv, total_derivative, step):
    for ladder in xrange(NLADDERS):
        trbid = ladder/ NLADDERSTRB
        ladderid = ladder % NLADDERSTRB
        layer= 5-ladderid/4
        if trbid in [0,1,4,5]:
            coor = "y"
        else:
            coor = "x"
        alignment[ladder]["dthetax"] -= angledrv["l%d"%layer][coor][0] * step / total_derivative
        alignment[ladder]["dthetay"] -= angledrv["l%d"%layer][coor][1] * step / total_derivative
        alignment[ladder]["dthetaz"] -= angledrv["l%d"%layer][coor][2] * step / total_derivative



def add_layers_to_alignment(alignment, angledrv, total_derivative, step):
    
    for sensor in xrange(NLADDERS * NSENSORS):
        ladder = sensor / 4
        trbid = ladder/ NLADDERSTRB
        ladderid = ladder % NLADDERSTRB
        layer= 5-ladderid/4
        if trbid in [0,1,4,5]:
            coor = "y"
        else:
            coor = "x"
        alignment[sensor]["dx"] -= angledrv["l%d"%layer][coor][0] * step / total_derivative
        alignment[sensor]["dy"] -= angledrv["l%d"%layer][coor][1] * step / total_derivative
        alignment[sensor]["dz"] -= angledrv["l%d"%layer][coor][2] * step / total_derivative
        alignment[sensor]["dthetax"] -= angledrv["l%d"%layer][coor][3] * step / total_derivative
        alignment[sensor]["dthetay"] -= angledrv["l%d"%layer][coor][4] * step / total_derivative
        alignment[sensor]["dthetaz"] -= angledrv["l%d"%layer][coor][5] * step / total_derivative


def add_posnegx_to_alignment(alignment, posneg):
    for ladder in xrange(NLADDERS):
        trbid = ladder/ NLADDERSTRB
        if trbid not in [2,3,6,7]:
             continue
        if trbid in [2,3]:
            alignment[ladder]["dx"] += posneg
        elif trbid in [6,7]:
            alignment[ladder]["dx"] -= posneg


def add_posnegy_to_alignment(alignment, posneg):
    for ladder in xrange(NLADDERS):
        trbid = ladder/ NLADDERSTRB
        if trbid not in [0,1,4,5]:
             continue
        if trbid in [0,1]:
            alignment[ladder]["dy"] += posneg
        elif trbid in [4,5]:
            alignment[ladder]["dy"] -= posneg


def add_new_offsets_to_alignment(alignment, dchish_dpars, total_derivative, step, anglescalefactorXY = 1.0,anglescalefactorZ = 1.0, alphascalefactor = 1.0, betascalefactor = 1.0, forbiddenladders = [], alignsensors = False):
    n = NLADDERS
    if alignsensors: n*=NSENSORS
    for ladder in xrange(n):
        if ladder in forbiddenladders: continue
	#print "\n\n\n\n", alignment
        alignment[ladder]["dx"] -= dchish_dpars[ladder][0] * step / total_derivative
        alignment[ladder]["dy"] -= dchish_dpars[ladder][1] * step / total_derivative
        alignment[ladder]["dz"] -= dchish_dpars[ladder][2] * step / total_derivative
        alignment[ladder]["dthetax"] -= dchish_dpars[ladder][3] * step / total_derivative  * anglescalefactorXY
        alignment[ladder]["dthetay"] -= dchish_dpars[ladder][4] * step / total_derivative  * anglescalefactorXY
        alignment[ladder]["dthetaz"] -= dchish_dpars[ladder][5] * step / total_derivative  * anglescalefactorZ
        """
        alignment[ladder]["alphax"] -= dchish_dpars[ladder][6] * step / total_derivative  * alphascalefactor
        alignment[ladder]["alphay"] -= dchish_dpars[ladder][7] * step / total_derivative  * alphascalefactor
        alignment[ladder]["betax"] -= dchish_dpars[ladder][6] * step / total_derivative  * betascalefactor
        alignment[ladder]["betay"] -= dchish_dpars[ladder][7] * step / total_derivative  * betascalefactor
        """

def add_new_offsets_to_alignment_perelement(alignment, dchish_dpars, total_derivative, step, anglescalefactorXY = 1.0,anglescalefactorZ = 1.0, alphascalefactor = 1.0, betascalefactor = 1.0, forbiddenladders = []):
    for ladder in xrange(len(step)):
        if ladder in forbiddenladders: continue
        alignment[ladder]["dx"] -= dchish_dpars[ladder][0] * step [ladder] / total_derivative[ladder]
        alignment[ladder]["dy"] -= dchish_dpars[ladder][1] * step [ladder] / total_derivative[ladder]
        alignment[ladder]["dz"] -= dchish_dpars[ladder][2] * step [ladder] / total_derivative[ladder]
        alignment[ladder]["dthetax"] -= dchish_dpars[ladder][3] * step [ladder] / total_derivative[ladder]  * anglescalefactorXY
        alignment[ladder]["dthetay"] -= dchish_dpars[ladder][4] * step [ladder] / total_derivative[ladder]  * anglescalefactorXY
        alignment[ladder]["dthetaz"] -= dchish_dpars[ladder][5] * step [ladder] / total_derivative[ladder]  * anglescalefactorZ



def read_alignment_from_file(fname):
    f = open(fname, "r")
    lines = f.readlines()
    alignments = {}
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
 
        alignments[ladderid] = {"dx": dx, "dy":dy, "dz":dz, "dthetax": dthetax, "dthetay": dthetay, "dthetaz": dthetaz, "alphax":alphax, "betax":betax, "alphay":alphay, "betay":betay }
    f.close()
    return alignments

def eformat(f, prec = 9, exp_digits = 2):
    s = "%+.*e"%(prec, f)
    try:
        mantissa, exp = s.split('e')
    except:
         print "prec = ", prec, "   s=", s, "  exp_digits=",exp_digits
         raise Exception
    # add 1 to digits as 1 is taken by sign +/-
    return "%se%+0*d"%(mantissa, exp_digits+1, int(exp))


def PrinChiSqDerivatives(dchish_dpars, tracks_for_lad,alignsensors=False):
    print
    print "===========> Chi-square derivatives: " 
    print
    print "---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" 
    print " lad |          dchidx           dchidy           dchidz         dchidthx         dchidthx         dchidthz        alphax        alphay       betax       betay    |      #tracks"
    print "---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"
    n = NLADDERS
    if alignsensors: n*=NSENSORS
    for i in xrange(n):
        print "%3d  | %s  %s  %s  %s  %s  %s  %s  %s  %s  %s |  %10d"%tuple([i] + map(lambda x: eformat(x, 9, 2), list(dchish_dpars["derivatives"][i])) + [sum(sum(tracks_for_lad[i]))])
    print "---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"
    print


def PrinChiSqDerivatives4(dchish_dpars, total_derivative, step, tracks_for_lad):
    print
    print "===========> Chi-square derivatives: " 
    print
    print "------------------------------------------------------------------------------------------------------------------------------------------------------------------" 
    print " lad |          dchidx           dchidy           dchidz         dchidthx         dchidthx         dchidthz      |   totalder    |      step    |     tracks     |     "
    print "------------------------------------------------------------------------------------------------------------------------------------------------------------------"
    for i in xrange(len(step)):
        print "%3d  | %s  %s  %s  %s  %s  %s | %s | %s |  %10d"%tuple([i] + map(lambda x: eformat(x, 9, 2), list(dchish_dpars["derivatives"][i][:6]) + [total_derivative[i]] + [step[i]]) + [tracks_for_lad[i]])
    print "------------------------------------------------------------------------------------------------------------------------------------------------------------------"
    print


def write_alignment(fname, alignment, alignsensors=False):
    f=open(fname, "w")
    n = NLADDERS
    if alignsensors: n*=NSENSORS
    for ladder in xrange(n):
        al = alignment[ladder]
        line = "%5d , %9.1f , %9.1f , %9.1f , %9.1f , %9.1f , %9.1f , %s, %s, %s, %s \n"%(ladder, al["dx"], al["dy"], al["dz"], al["dthetax"], al["dthetay"], al["dthetaz"], eformat(al["alphax"]), eformat(al["alphay"]), eformat(al["betax"]), eformat(al["betay"]))
        f.write(line)
    f.close()
     

def test_x():
    step = 100.
    for i in xrange(1):
        print "\n\n\n\n\n\n"
        print "Iteration: ", i
        all_dchish_dpars, all_tracks_for_lad = run(trackchisqcut = TRACK_CHI_SQ_CUT)
        print "\n\n\n\n\n\n"
        print "================================="
        dci_dposneg = dchisq_dposnegx(all_dchish_dpars)
        print "dci_dposneg  = ", dci_dposneg
        print "chisq/tracks = ", all_dchish_dpars["chisq"] / all_dchish_dpars["ntracks"]
        print "================================="
        alignment = read_alignment_from_file(ALIGNMENT_TMP)
        add_posnegx_to_alignment(alignment, -step * np.sign(dci_dposneg)  )
        write_alignment(ALIGNMENT_TMP, alignment)


def test_y():
    step = 100.
    for i in xrange(1):
        print "\n\n\n\n\n\n"
        print "Iteration: ", i
        all_dchish_dpars, all_tracks_for_lad = run(trackchisqcut = TRACK_CHI_SQ_CUT)
        print "\n\n\n\n\n\n"
        print "================================="
        dci_dposneg = dchisq_dposnegy(all_dchish_dpars)
        print "dci_dposneg  = ", dci_dposneg
        print "chisq/tracks = ", all_dchish_dpars["chisq"] / all_dchish_dpars["ntracks"]
        print "================================="
        alignment = read_alignment_from_file(ALIGNMENT_TMP)
        add_posnegy_to_alignment(alignment, -step * np.sign(dci_dposneg)  )
        write_alignment(ALIGNMENT_TMP, alignment)


def test2():
    alignsensors = True
    #step = 0. # 25. #3.125   #100. 
    step = ALGSTEPFIRST  #3.125   #100. 

    previouschisq = 99999999999.
    total_derivative_previous = 0.0000000001
    all_dchish_dpars_old = None

    minimalstep = 0.01
    anglescalefactorXY = 5.  #2.  #0.2  # 100.0
    anglescalefactorZ  = 5.  #2.  #0.2  # 100.0
    alphascalefactor   = 10.* 1e-8
    betascalefactor  = 2. * 1e-10
    #while(True):
    all_steps = []
    previosder = None
    scalar = None


    for i in xrange(MAXALGITERATIONS):
    #for i in xrange(1):
        #print info
        print "\n\n\n\n\n\n"
        print "Iteration: ", i
        


        # do analysis 
        rerun = False
        if not i and RERUN_FIRST: rerun = True
        while(True):
            #all_dchish_dpars, all_tracks_for_lad = run( exitifstuck = True, alignsensors=alignsensors)
            all_dchish_dpars, all_tracks_for_lad = run(alltracksforladder = MERGED_OUT_DUMP, trackchisqcut = TRACK_CHI_SQ_CUT, exitifstuck = True, alignsensors=alignsensors, rerun=rerun)
            rerun = True
            #all_dchish_dpars, all_tracks_for_lad = run(alltracksforladder = MERGED_OUT_DUMP, trackchisqcut = TRACK_CHI_SQ_CUT, exitifstuck = True, alignsensors=alignsensors)
            if all_dchish_dpars is not None: break
        
        n = NLADDERS
        if alignsensors: n*=NSENSORS
        for lad in xrange(n): # <--- normalize to um and microns
            all_dchish_dpars["derivatives"][lad][0] /= 1000.
            all_dchish_dpars["derivatives"][lad][1] /= 1000.
            all_dchish_dpars["derivatives"][lad][2] /= 1000.
            all_dchish_dpars["derivatives"][lad][3] /= (1000000. / anglescalefactorXY )
            all_dchish_dpars["derivatives"][lad][4] /= (1000000. / anglescalefactorXY )
            all_dchish_dpars["derivatives"][lad][5] /= (1000000. / anglescalefactorZ  )
            all_dchish_dpars["derivatives"][lad][6] *= alphascalefactor
            all_dchish_dpars["derivatives"][lad][7] *= alphascalefactor
            all_dchish_dpars["derivatives"][lad][8] *= betascalefactor
            all_dchish_dpars["derivatives"][lad][9] *= betascalefactor
        

        #total_derivative = sum([ np.power( all_dchish_dpars["derivatives"][lad][par],2) for lad in xrange(n) for par in [2,3,4]  ])
        total_derivative = sum([ np.power( all_dchish_dpars["derivatives"][lad][par],2) for lad in xrange(n) for par in xrange(6)  ])
        #total_derivative = sum([ np.power( all_dchish_dpars["derivatives"][lad][par],2) for lad in xrange(n) for par in xrange(NALIGNPARSLAD)  ])
        #total_derivative = sum([ np.power( all_dchish_dpars["derivatives"][lad][par],2) for lad in xrange(NLADDERS) for par in [0,1,2]  ]) # coordinates
        #total_derivative = sum([ np.power( all_dchish_dpars["derivatives"][lad][par],2) for lad in xrange(NLADDERS) for par in [3,4,5]  ]) # angles
        total_derivative = np.power( total_derivative, 0.5)  


        #print "test3"
        #print "\n\n\n\n all_dchish_dpars['chisq'] = ", all_dchish_dpars['chisq']
        #print " total_derivative_previous = ", total_derivative_previous
        #print "\n\n\n\n"


        # update step
        #if previouschisq < all_dchish_dpars['chisq']: step/=2.
        divided_step = False
        #if (previouschisq - all_dchish_dpars['chisq']) < step * total_derivative_previous: # * 0.66:  #0.5:  

        if previosder is not None:
            scalar = sum([  all_dchish_dpars["derivatives"][lad][par] * previosder[lad][par] for lad in xrange(n) for par in xrange(6)  ])/ (total_derivative_previous * total_derivative)

        if previosder is not None and scalar  < 0.5:
            step/=2.
            divided_step = True
        elif len(filter(lambda x: x==all_steps[-1], all_steps[-10:-1] )) ==9:
        #elif len(all_steps)>3 and all_steps[-4] == all_steps[-3] and all_steps[-3] == all_steps[-2] and  all_steps[-2] == all_steps[-1]: 
            step*=2.


        # update previous chisq
        previouschisq = all_dchish_dpars['chisq']
        all_steps.append(step)  
 

        # read alignment
        print "Applying step: ", step
        alignment = read_alignment_from_file(ALIGNMENT_TMP)
        add_new_offsets_to_alignment(alignment, all_dchish_dpars["derivatives"], total_derivative, step , 
                                      anglescalefactorXY = anglescalefactorXY, 
                                      anglescalefactorZ  = anglescalefactorZ, 
                                      alphascalefactor = alphascalefactor, 
                                      betascalefactor  = betascalefactor,   
                                      alignsensors=alignsensors)
        

	"""
        # ------------------------  TEST --------------------------
        # check for overshoot
        forbiddenladders = []
        if all_dchish_dpars_old is not None and not divided_step:
            n = NLADDERS
            if alignsensors: n*=NSENSORS
            for lad in xrange(n):
                overshoot = True
                for par in xrange(NALIGNPARSLAD):
                    old = all_dchish_dpars_old["derivatives"][lad][par]
                    new = all_dchish_dpars["derivatives"][lad][par]
                    if old==0 and new==0: continue
                    if old*new < 0 and abs(new)> abs(old) * 0.5: continue
                    overshoot = False
                    break  
                if overshoot: forbiddenladders.append(lad)
        if forbiddenladders:
            total_derivative_notcorrected = total_derivative
            total_derivative = sum([ np.power( all_dchish_dpars["derivatives"][lad][par],2) for lad in xrange(n) if lad not in forbiddenladders for par in xrange(NALIGNPARSLAD)  ])
            total_derivative = np.power( total_derivative, 0.5)  
            alignment = read_alignment_from_file(ALIGNMENT_TMP)
            add_new_offsets_to_alignment(alignment, all_dchish_dpars["derivatives"], total_derivative, step , 
                                           anglescalefactor = anglescalefactor, 
                                           alphascalefactor = alphascalefactor, 
                                           betascalefactor  = betascalefactor,
                                           forbiddenladders =forbiddenladders, alignsensors=alignsensors)
            print
            print "=========================================================="
            print "using corrected step!"
            print "total_derivative = ", total_derivative
            print "total_derivative_notcorrected = ", total_derivative_notcorrected
            print "=========================================================="
            print 
        print
        print "forbiddenladders = ", forbiddenladders
        print 
        all_dchish_dpars_old = all_dchish_dpars
        # ------------------------  TEST --------------------------
	"""


        # write alignment to file
        write_alignment(ALIGNMENT_TMP, alignment, alignsensors=alignsensors)
        total_derivative_previous = total_derivative
        previosder = all_dchish_dpars["derivatives"]
        
        
        # print info        
        print "\n\n\n\n\n\n"
        PrinChiSqDerivatives(all_dchish_dpars, all_tracks_for_lad, alignsensors=alignsensors)
        print "================================="
        print "current step                = ", step
        print "all_dchish_dpars['chisq']   = ", all_dchish_dpars['chisq']
        print "all_dchish_dpars['ntracks'] = ", all_dchish_dpars['ntracks']
        print "total_derivative            = ", total_derivative
        print "scalar                      = ", scalar
        print "================================="

 
        #
        print "changing permissions"
        ft = glob(OUTPUT_PATH+"/*")
        for f in ft:
            os.system("chmod 777 "+f)
       
       
        

        # check if limit reached 
        if minimalstep>step:
            break



def test3():
    alignsensors = True
    step = ALGSTEPFIRST #100.0 #10.0
    previouschisq = 99999999999.
    total_derivative_previous = 0.0000000001
    minimalstep = 1.
    all_steps = []
    for i in xrange(MAXALGITERATIONS):
        print "\n\n\n\n\n\n"
        print "Iteration: ", i
        
        # do analysis 
        rerun = False
        if not i and RERUN_FIRST: rerun = True
        while(True):
            all_dchish_dpars, all_tracks_for_lad = run(trackchisqcut = TRACK_CHI_SQ_CUT, exitifstuck = True, alignsensors=alignsensors, rerun=rerun)
            rerun = True
            if all_dchish_dpars is not None: break

        """
        # do analysis 
        while(True):
            all_dchish_dpars, all_tracks_for_lad = run(trackchisqcut = TRACK_CHI_SQ_CUT, exitifstuck = True, alignsensors=True)
            if all_dchish_dpars is not None: break
        """

        for lad in xrange(NLADDERS*NSENSORS): # <--- normalize to um and microns
            all_dchish_dpars["derivatives"][lad][0] /= 1000.
            all_dchish_dpars["derivatives"][lad][1] /= 1000.
            all_dchish_dpars["derivatives"][lad][2] /= 1000.
            all_dchish_dpars["derivatives"][lad][3] /= 1000000. # / anglescalefactor )
            all_dchish_dpars["derivatives"][lad][4] /= 1000000. # / anglescalefactor )
            all_dchish_dpars["derivatives"][lad][5] /= 1000000. # / anglescalefactor )



        layerdrv = dchisq_layers(all_dchish_dpars)
         
        total_derivative = sum([np.power(d,2) for l in layerdrv.values() for c in l.values() for d in c])
        total_derivative = np.power( total_derivative, 0.5)  


        # update step
        #if previouschisq < all_dchish_dpars['chisq']: step/=2.
        if (previouschisq - all_dchish_dpars['chisq']) < step * total_derivative_previous * 0.66: #* 0.5:
            step/=2.
        elif len(filter(lambda x: x==all_steps[-1], all_steps[-10:-1] )) ==9:
            step*=2.

        # update previous chisq
        previouschisq = all_dchish_dpars['chisq']
        all_steps.append(step)  
 





        # create new alignment file
        alignment = read_alignment_from_file(ALIGNMENT_TMP)
        add_layers_to_alignment(alignment, layerdrv, total_derivative, step)


        write_alignment(ALIGNMENT_TMP, alignment, alignsensors=True)
        #write_alignment(ALIGNMENT_TMP, alignment)
        total_derivative_previous = total_derivative
        
        # print info        
        print "\n\n\n\n\n\n"
        #PrinChiSqDerivatives(all_dchish_dpars, all_tracks_for_lad)
        print "angledrv="
        for k in sorted(layerdrv.keys()):
            print k, "   :   ", layerdrv[k]
        print "================================="
        print "previous step               = ", step
        print "all_dchish_dpars['chisq']   = ", all_dchish_dpars['chisq']
        print "all_dchish_dpars['ntracks'] = ", all_dchish_dpars['ntracks']
        print "total_derivative            = ", total_derivative
        print "================================="

       
        

        # check if limit reached 
        #if minimalstep>step:
        #    break



def test4():
    alignsensors = True

    nelements = NLADDERS 
    if alignsensors:
        nelements = nelements * NSENSORS

    step = [ALGSTEPFIRST]*nelements
    total_derivative_previous = [99999999999.] * nelements

    minimalstep = 0.1
    anglescalefactorXY = 5.  #2.  #0.2  # 100.0
    anglescalefactorZ  = 5.  #2.  #0.2  # 100.0
    alphascalefactor   = 10.* 1e-8
    betascalefactor  = 2. * 1e-10
    all_steps = []


    for i in xrange(MAXALGITERATIONS):
        print "\n\n\n\n\n\n"
        print "Iteration: ", i
        


        # do analysis 
        rerun = False
        if not i and RERUN_FIRST: rerun = True
        while(True):
            all_dchish_dpars, all_tracks_for_lad = run(trackchisqcut = TRACK_CHI_SQ_CUT, exitifstuck = True, alignsensors=alignsensors, rerun=rerun)
            rerun = True
            if all_dchish_dpars is not None: break
        
        for lad in xrange(nelements): # <--- normalize to um and microns
            all_dchish_dpars["derivatives"][lad][0] /= 1000.
            all_dchish_dpars["derivatives"][lad][1] /= 1000.
            all_dchish_dpars["derivatives"][lad][2] /= 1000.
            all_dchish_dpars["derivatives"][lad][3] /= (1000000. / anglescalefactorXY )
            all_dchish_dpars["derivatives"][lad][4] /= (1000000. / anglescalefactorXY )
            all_dchish_dpars["derivatives"][lad][5] /= (1000000. / anglescalefactorZ  )
            all_dchish_dpars["derivatives"][lad][6] *= alphascalefactor
            all_dchish_dpars["derivatives"][lad][7] *= alphascalefactor
            all_dchish_dpars["derivatives"][lad][8] *= betascalefactor
            all_dchish_dpars["derivatives"][lad][9] *= betascalefactor
        

        total_derivative = [sum([ np.power( all_dchish_dpars["derivatives"][lad][par],2) for par in xrange(6)  ]) for lad in xrange(nelements)]
        total_derivative = [np.power( x, 0.5) for x in total_derivative] 




        # update step
        #if previouschisq < all_dchish_dpars['chisq']: step/=2.
        """
        for el in xrange(len(total_derivative)):
            divided_step = False
            if total_derivative_previous[el]<total_derivative[el]: # * 0.66:  #0.5:
                step[el]/=2.
                divided_step = True
            elif len(filter(lambda x: x==all_steps[-1][el], [stp[el] for stp in all_steps[-10:-1]]   )) ==9:
                step[el]*=2.
        """


        # update previous chisq
        #previouschisq = all_dchish_dpars['chisq']
        all_steps.append(list(step))  
 

        # read alignment
        #print "\nApplying step: ["
        #for el in step: print el
        #print "]\n"
        
        alignment = read_alignment_from_file(ALIGNMENT_TMP)
        add_new_offsets_to_alignment_perelement(alignment, all_dchish_dpars["derivatives"], total_derivative, step , 
                                      anglescalefactorXY = anglescalefactorXY, 
                                      anglescalefactorZ  = anglescalefactorZ, 
                                      alphascalefactor = alphascalefactor, 
                                      betascalefactor  = betascalefactor)
                                      
        



        # write alignment to file
        write_alignment(ALIGNMENT_TMP, alignment, alignsensors=alignsensors)
        total_derivative_previous = list(total_derivative)
        
        # print info        
        print "\n\n\n\n\n\n"
        PrinChiSqDerivatives4(all_dchish_dpars, total_derivative, step, all_tracks_for_lad)
        print "================================="
        #print "current step                = ", step
        #print "current step: ["
        #for el in step: print el
        #print "]"
        print "all_dchish_dpars['chisq']   = ", all_dchish_dpars['chisq']
        print "all_dchish_dpars['ntracks'] = ", all_dchish_dpars['ntracks']
        #print "total_derivative            = ", total_derivative
        print "================================="


        #
        print "changing permissions"
        ft = glob(OUTPUT_PATH+"/*")
        for f in ft:
            os.system("chmod 777 "+f)
       

       
        

        # check if limit reached 
        if not [x for x in step if x>minimalstep]:
            break



test2()
