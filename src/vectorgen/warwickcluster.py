import os
import tempfile
import time
import socket

###############################################################################

class ClusterDefs:
    #see http://www2.warwick.ac.uk/fac/sci/physics/research/epp/internal/computing/hardware/cluster/batch/
    queueNames = ["express", # 10 minutes
                  "short", # 30 minutes
                  "medium", # 2 hours
                  "long", # 12 hours
                  "xlong", # 48 hours
                  "xxl", # 96 hours
                  ]

###############################################################################

class ClusterJobException(Exception):
    pass

###############################################################################

class ClusterJob:

    #a globally accessible time string that all jobs created this session may use to label themselves.
    initialiseTime = time.strftime("%Y-%m-%d.%H%M")

    def __init__(self, name, queue, cmd, workingDirectory="./", isTest=False, isVerbose=False, isQuiet=False,
                 memoryLimit=4000,
                 swapLimit=4000,
    ):
        self.cmd = cmd
        self.jobName = name
        self.queue = queue
        self.workingDirectory = self.getAbsolutePath(workingDirectory)
        self.isTest = isTest
        self.isQuiet = isQuiet
        self.isVerbose = isVerbose
        self.scriptPath = None
        self.timeStr = ClusterJob.initialiseTime
        #set memory limits
        self.memoryLimit = memoryLimit # in MB
        self.swapLimit = swapLimit # in MB
        self.sanitiseInputs()

    def __str__(self):
        optStr = ",".join([self.jobName,self.queue])
        return "ClusterJob("+optStr+")"

    def getAbsolutePath(self,path):
        return os.path.abspath(os.path.expanduser(path))

    def sanitiseInputs(self):
        #ensure that working directory exists
        if not ( os.path.exists(self.workingDirectory) and os.path.isdir(self.workingDirectory) ):
            raise ClusterJobException("working directory does not exist",self.workingDirectory)
        #ensure that queue name is ok
        if self.queue not in ClusterDefs.queueNames:
            raise ClusterJobException("given unknown queue",self.queue, ClusterDefs.queueNames)
        #ensure that job name is not null
        if not self.jobName:
            raise ClusterJobException("job name not correctly set: it should be a none null string",self.jobName)

    def submit(self):
        self.generateScript()
        cmd = "bsub < "+self.scriptPath
        if self.isTest or not self.isQuiet:
            print str(self),"executing:",cmd
        if self._isLoginNode():
            os.system(cmd)
        elif not self.isTest:
            raise ClusterJobException("host not recognised as cluster login node",socket.getfqdn())

    def isLoginNode(self):
        return "epp-ui01" in socket.getfqdn()

    def generateScript(self):
        #create a unique file to write the script to.
        #fileHandle = tempfile.NamedTemporaryFile(prefix=self.jobName, suffix=".sh", dir=self.workingDirectory, delete=False)
        #fileName = fileHandle.name
        fileHandle = None
        counter = 0
        while not fileHandle:
            fileName = "%(pwd)s/%(name)s.%(time)s.%(count)03d.sh" % {"pwd":self.workingDirectory,
                                                                  "name":self.jobName,
                                                                  "time":self.timeStr,
                                                                  "count":counter,
                                                                  }
            counter += 1
            if not os.path.exists(fileName):
                fileHandle = open(fileName,"w")
        #make sure we are using bash
        print >>fileHandle,"#!/bin/bash"
        #set batch system options
        self.stdoutFile = os.path.splitext(fileName)[0]+".out"
        options = ["-q "+self.queue,
                   "-J "+self.jobName,
                   " -o "+self.stdoutFile,
                   "-M "+str(self.memoryLimit),
                   "-v "+str(self.swapLimit),
                   ]
        for opt in options:
            print >>fileHandle,"#BSUB",opt
        #change working directory to current working directory
        if self.workingDirectory is None:
            self.workingDirectory = os.getcwd()
        if not self.workingDirectory == "":
            print >>fileHandle,"cd "+self.workingDirectory
        #finally, add the command to be run
        print >>fileHandle,self.cmd
        #clean up resources and record the script filename
        fileHandle.close()
        #return the absolute path of the file we created
        self.scriptPath = fileName
        return fileName

###############################################################################

def _unitTestClusterJob():
    cmd = "echo Hello World && sleep 5"
    job = ClusterJob("unitTestClusterJob", "express",cmd)
    job.submit()

###############################################################################

def main():
    _unitTestClusterJob()
    return

if __name__ == "__main__":
    main()

