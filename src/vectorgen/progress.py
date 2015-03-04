"""Provides utilities to print the status of a loop as it is processed.
"""

import datetime
import StringIO
import sys

###############################################################################

def printprogress(name, num_entries, iterable, update=True):
    """Prints to stdout every 1 or 10% way through the loop.
    
    nEntries should be the size of the inputer iterable.
    When you loop over this generator the generator will return elements in the 
    input iterable and print information along the way.
    It prints out the progress at several milestones. 
    The milestones are 0, 1, 2 ... 10% complete and then 20, 30 ... 100%.
    
    If the update option is set the previous message is deleted and updated with
    new information. In this case updates are given every 1% of progress.
    
    The information printed is:
    * The input name of the generator.
    * N entries processed and total entries.
    * Time elapsed and estimated time until completion.
    * Current date and time.
    
    Returns a generator.
    """
    nEntries = num_entries
    startTime = _getNow()
    #print a header before any iteration is done
    print "----",name,"(entries="+str(nEntries),"begin="+_getTimeStampString(startTime,date=True)+")"
    header = " ".join([name,"  %"]) + " (" + ", ".join(["entries","elapsed/remaining","timestamp"]) + ")"
    print header
    i = 0
    nextUpdatePercentage = 0
    lastMessageSize = 0
    lastForcedNewLine = startTime
    
    for entry in iterable:
        #calculate percentage complete (as an integer)
        done = (100*i)/nEntries
        if done >= nextUpdatePercentage:
            #print an update of progress
            #fraction of job complete
            doneFracStr = str("%2.0f"%(done))+"%"
            try:
                #Add other information
                #get string for current time
                nowTime = _getNow()
                #nowStr = nowTime.strftime("%Y.%m.%d-%H:%M")
                nowStr = _getTimeStampString(nowTime)
                #get string for elapsed time
                elapsedTime = nowTime - startTime
                elapsedSeconds = _getTotalSecondsFromTimeDelta(elapsedTime)
                elapsedStr = _formatDeltaSeconds(elapsedSeconds)
                #get string for estimated time remaining
                try:
                    doneFloat = (100.0*float(i))/float(nEntries)
                    remainingSeconds = int(abs(float(elapsedSeconds) * (100.0-doneFloat)/doneFloat))
                    remainingStr = _formatDeltaSeconds(remainingSeconds)
                except ZeroDivisionError:
                    remainingStr = "???"
                timeStr = "/".join([elapsedStr,remainingStr])
                #get a string showing how many entries have been processed
                unit = ""
                scale = 1
                #for very large numbers shorten them by adding SI prefixes
                if nEntries>=100000:
                    unit = "k"
                    scale = 1000
                if nEntries>=100000000:
                    unit = "M"
                    scale = 1000*1000
                if nEntries>=100000000000:
                    unit = "G"
                    scale = 1000*1000*1000
                n = str(max(1,len("%d"%(nEntries/scale))))
                formatStr = "{0:"+n+"d}"+unit+"/{1:d}"+unit
                entriesStr = formatStr.format(i/scale,nEntries/scale)
                #finally join all strings to form the printed message
                message = " ".join([name,doneFracStr]) + " (" + ", ".join([entriesStr,timeStr,nowStr]) + ")"
            except Exception:
                raise
                #this should never happen, but if there is a rare bug in this 
                # code we don't want it end a long job.
                # if any exception is thrown just print the job name and 
                # fraction of job complete.
                message = " ".join([name,doneFracStr])
            #never print first entry as there is no time remaining estimate   
            
            if update:
                #delete the last message
                if lastMessageSize > 0:
                    delCommand = "\b"*(lastMessageSize)
                    sys.stdout.write(delCommand)
                    sys.stdout.flush()
                #print the new message
                sys.stdout.write(message)
                #store the new message size to print the next line
                lastMessageSize = len(message)
                #force a new line every 10 minutes to keep a record in the 
                #log file for long jobs
                if _getTotalSecondsFromTimeDelta(nowTime - lastForcedNewLine) > 60*10:
                    lastForcedNewLine = nowTime
                    sys.stdout.write("\n")
                    lastMessageSize = 0
                if i>0:#this condition forces it to print a line at 0% and first entry
                    #set next update step 1 percent
                    nextUpdatePercentage += 1
                #force python to flush it's buffer
                sys.stdout.flush()
            else:
                #print the message and start a new line 
                sys.stdout.write(message)
                sys.stdout.write("\n")
                #force python to flush it's buffer
                sys.stdout.flush()
                #determine next entry to print an update
                if i>0:#this condition forces it to print a line at 0% and first entry
                    if nextUpdatePercentage < 10:
                        #first 10%, print every 1% processed
                        nextUpdatePercentage += 1
                    else:
                        #print every 10% processed
                        nextUpdatePercentage += 10
        yield entry
        i += 1
    #if in update, there is usually a need to add a new line
    if update:
        sys.stdout.write("\n")
    #finally write a message to indicate that the loop finished
    print "-----",name,"(end="+_getTimeStampString(_getNow(),date=True)+")"
    #print name,"finished at",_getTimeStampString(_getNow(),date=True)
    return

def _getTotalSecondsFromTimeDelta(timeDelta):
    r = -1
    try:
        r = timeDelta.total_seconds()
    except AttributeError:
        r = timeDelta.seconds
    return r

        
def _formatDeltaSeconds(seconds):
    sign = ""
    if seconds < 0:
        sign = "-"
    seconds = abs(seconds)
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    r = "%s%dh%dm%ds" % (sign,hours,minutes,seconds)
    if hours>0:
        r = "%s%dh%02dm" % (sign,hours,minutes)
    elif minutes>0:
        r = "%s%dm%02ds" % (sign,minutes,seconds)
    else:
        r = "%s%02ds" % (sign,seconds)
    return r

def _getNow():
    return datetime.datetime.now()

def _getTimeStampString(time, date=False):
    seconds = _systemTimeOffsetFromUtc
    sign = "+"
    if seconds < 0:
        sign = "-"
    seconds = abs(seconds)
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    utcOffset = ("{0:s}{1:02d}:{2:02d}").format(sign,hours,minutes)
    timeStampFormat = "%H:%M:%S"+utcOffset
    if date:
        timeStampFormat += " %Y.%m.%d"
    return time.strftime(timeStampFormat)

###############################################################################

def _testPrintProgress():
    import time
    nEvents = 100000
    def dummyIterable():
        total = 5.0
        for i in xrange(nEvents):
            time.sleep(total/float(nEvents))
            yield i
        return
    #try update mode
    iterable = dummyIterable()
    for i in printprogress("_testPrintProgress", nEvents, iterable, update=True):
        #do nothing
        pass 
    #try without update mode
    iterable = dummyIterable()    
    for i in printprogress("_testPrintProgress", nEvents, iterable, update=False):
        #do nothing
        pass 
    return

###############################################################################

_systemTimeOffsetFromUtc = int(round( _getTotalSecondsFromTimeDelta(_getNow()-datetime.datetime.utcnow()) ))

###############################################################################

def main():
    _testPrintProgress()

###############################################################################

if __name__ == "__main__":
    main()
