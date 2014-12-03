import fcntl
import time

###############################################################################

class FileLock:
    def __init__(self, name=".lockFile_vectorgen", verbose=False, timeout=60*5):
        self.fname = name
        self.lockfile = None
        self.verbose = verbose
        self.waitTime = 5
        self.waitAttempts = max(1,timeout/5)
    def lock(self):
        hasLock = False
        #remove null character as this sometimes causes problems
        self.fname = "".join(self.fname.split("\x00"))
        self.lockFile = open(self.fname, "w")
        for attempt in xrange(self.waitAttempts):
            try:
                fcntl.lockf(self.lockFile, fcntl.LOCK_EX | fcntl.LOCK_NB)
                #successfully got lock
                hasLock = True
                break
            except IOError:
                if self.verbose:
                    print "waiting for lock on ",self.fname,attempt,"/",self.waitAttempts
                time.sleep(self.waitTime)
        if not hasLock and self.verbose:
            print "warning : failed to get lock on ",self.fname
        return hasLock

    def release(self):
        if self.lockFile:
            self.lockFile.close()
            self.lockFile = None
            
###############################################################################

def main():
    lf = FileLock(verbose=True)
    lf.lock()
    for i in xrange(10):
        print "have lock, waiting",i
        time.sleep(1)
    lf.release()
    return

###############################################################################

if __name__ == "__main__":
    main()
