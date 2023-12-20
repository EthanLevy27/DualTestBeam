from array import *
import argparse


# python massjobs2.py -g fSCEPonly

argParser = argparse.ArgumentParser()
argParser.add_argument("-g", "--geometry", help="geometry code")

args = argParser.parse_args()
print("args=%s" % args)

print("args.name=%s" % args.geometry)


hostarea="/data/users/eno/CalVision/dd4hep/DD4hep/examples/DualTestBeam/compact/jobs/"
outputarea="/data/users/eno/CalVision/dd4hep/DD4hep/examples/DualTestBeam/compact/output/"




nenergy=2
energies=[10,30]
name="condor-executable-"+args.geometry+"-"

# create the .sh files for electrons
i=0
while (i<nenergy):
    print(i)
    shfile = open(name+str(energies[i])+'_GeV-e.sh',"w")

    shfile.write('#!/bin/bash'+'\n')
    shfile.write('cd /data/users/eno/CalVision/dd4hep/DD4hep/examples/DualTestBeam/compact/'+'\n')
    shfile.write('START_TIME=`/bin/date`'+'\n')
    shfile.write('echo "started at $START_TIME"'+'\n')
    shfile.write('echo "started at $START_TIME on ${HOSTNAME}"'+'\n')
    shfile.write('source /cvmfs/sft.cern.ch/lcg/views/LCG_102b/x86_64-centos7-gcc11-opt/setup.sh'+'\n')
    shfile.write('echo "ran setup"'+'\n')
    shfile.write('source  /data/users/eno/CalVision/dd4hep/DD4hep/bin/thisdd4hep.sh'+'\n')
    shfile.write('echo "ran thisdd4hep"'+'\n')
    shfile.write('ddsim --compactFile=/home/eno/CalVision/dd4hep/DD4hep/examples/DualTestBeam/compact/DRFSCEPonly.xml --runType=batch -G --steeringFile /home/eno/CalVision/dd4hep/DD4hep/examples/DualTestBeam/compact/SCEPCALsteering.py --outputFile='+outputarea+'out_'+args.geometry+'_20GeV_e-_100.root --part.userParticleHandler='' -G --gun.position="0.,0.,-210*cm" --gun.direction "0 0.05 0.99875" --gun.energy "'+str(energies[i])+'*GeV" --gun.particle="e-" -N 100 >& '+outputarea+'sce_e_fSCEPonly_'+str(energies[i])+'.log'+'\n')
    shfile.write('exitcode=$?'+'\n')
    shfile.write('echo ""'+'\n')
    shfile.write('END_TIME=`/bin/date`'+'\n')
    shfile.write('echo "finished at $END_TIME"'+'\n')
    shfile.write('exit $exitcode'+'\n')

    shfile.close()
    i=i+1
    print("file closed")

# create the .sh files for pions
i=0
while (i<nenergy):
    print(i)
    shfile = open(name+str(energies[i])+'_GeV-pi.sh',"w")

    shfile.write('#!/bin/bash'+'\n')
    shfile.write('cd /data/users/eno/CalVision/dd4hep/DD4hep/examples/DualTestBeam/compact/'+'\n')
    shfile.write('START_TIME=`/bin/date`'+'\n')
    shfile.write('echo "started at $START_TIME"'+'\n')
    shfile.write('echo "started at $START_TIME on ${HOSTNAME}"'+'\n')
    shfile.write('source /cvmfs/sft.cern.ch/lcg/views/LCG_102b/x86_64-centos7-gcc11-opt/setup.sh'+'\n')
    shfile.write('echo "ran setup"'+'\n')
    shfile.write('source  /data/users/eno/CalVision/dd4hep/DD4hep/bin/thisdd4hep.sh'+'\n')
    shfile.write('echo "ran thisdd4hep"'+'\n')
    shfile.write('ddsim --compactFile=/home/eno/CalVision/dd4hep/DD4hep/examples/DualTestBeam/compact/DRFSCEPonly.xml --runType=batch -G --steeringFile /home/eno/CalVision/dd4hep/DD4hep/examples/DualTestBeam/compact/SCEPCALsteering.py --outputFile='+outputarea+'out_fSCEPonly_20GeV_pi-_100.root --part.userParticleHandler='' -G --gun.position="0.,0.,-210*cm" --gun.direction "0 0.05 0.99875" --gun.energy "'+str(energies[i])+'*GeV" --gun.particle="pi-" -N 100 >& '+outputarea+'sce_pi_fSCEPonly_'+str(energies[i])+'.log'+'\n')
    shfile.write('exitcode=$?'+'\n')
    shfile.write('echo ""'+'\n')
    shfile.write('END_TIME=`/bin/date`'+'\n')
    shfile.write('echo "finished at $END_TIME"'+'\n')
    shfile.write('exit $exitcode'+'\n')

    shfile.close()
    i=i+1
    print("file closed")


# create the .jdl files for electrons
i=0
while (i<nenergy):
    print(i)
    jdlfile = open(name+str(energies[i])+'-e.jdl',"w")
    jdlfile.write("universe = vanilla"+'\n')
    jdlfile.write("Executable ="+name+str(energies[i])+"_GeV-e.sh"+'\n')
    jdlfile.write("should_transfer_files = NO"+'\n')
    jdlfile.write("Requirements = TARGET.FileSystemDomain == \"privnet\""+'\n')
    jdlfile.write("Output = "+hostarea+name+str(energies[i])+"-e_sce_$(cluster)_$(process).stdout"+'\n')
    jdlfile.write("Error = "+hostarea+name+str(energies[i])+"-e_sce_$(cluster)_$(process).stderr"+'\n')
    jdlfile.write("Log = "+hostarea+name+str(energies[i])+"-e_sce_$(cluster)_$(process).condor"+'\n')
    jdlfile.write("Arguments = SCE"+'\n')
    jdlfile.write("Queue 1"+'\n')
    jdlfile.close()
    i=i+1
    print("file closed")

# create the .jdl files for pins
i=0
while (i<nenergy):
    print(i)
    jdlfile = open(name+str(energies[i])+'-pi.jdl',"w")
    jdlfile.write("universe = vanilla"+'\n')
    jdlfile.write("Executable ="+name+str(energies[i])+"_GeV-pi.sh"+'\n')
    jdlfile.write("should_transfer_files = NO"+'\n')
    jdlfile.write("Requirements = TARGET.FileSystemDomain == \"privnet\""+'\n')
    jdlfile.write("Output = "+hostarea+name+str(energies[i])+"-pi_sce_$(cluster)_$(process).stdout"+'\n')
    jdlfile.write("Error = "+hostarea+name+str(energies[i])+"-pi_sce_$(cluster)_$(process).stderr"+'\n')
    jdlfile.write("Log = "+hostarea+name+str(energies[i])+"-pi_sce_$(cluster)_$(process).condor"+'\n')
    jdlfile.write("Arguments = SCE"+'\n')
    jdlfile.write("Queue 1"+'\n')
    jdlfile.close()
    i=i+1
    print("file closed")



# create the submitter file
f = open("massjobs.sh",'w')
i=0
while (i<nenergy):
    f.write("condor_submit "+name+str(energies[i])+'-e.jdl'+'\n')
    f.write("condor_submit "+name+str(energies[i])+'-pi.jdl'+'\n')
    i=i+1
f.write("condor_q")
f.close()



