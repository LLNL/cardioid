#!/usr/bin/env python

import os
import sys

testroot = ""

class Profile:
    def __init__(self,name):
        self.name = name
        #get the tags from the file

class Test:
    def __init__(self,name):
        self.name = name

    def directory(self):
        return os.join("tests",name)
        
    def execute(self, profile, build, action):
        myscript = """
        source %(testroot)s/loadTest %(build)s %(profile)s %(testname)s;
        %(action)s;
        """ % {
            "profile" : profile.name,
            "build" : build,
            "action" : action,
            "testname" : self.name,
            "testroot" : testroot,
            }
        bash_pipe = os.popen("bash", "w", len(myscript))
        print >>bash_pipe, myscript

def nameFromTestDir(dirname):
    return os.path.relpath(dirname,os.path.join(testroot, "tests"))
        
def main():
    global testroot
    testroot = os.path.dirname(os.path.realpath(__file__))
    print testroot
    
    import sys
    import argparse
    ap = argparse.ArgumentParser(description="Run all the tests for the Cardioid framework.")
    ap.add_argument("--profile", "-p",
                    help="Runtime profile to use for the tests.",
                    type=str,
                    default="seq",
                    )
    ap.add_argument("--build", "-b",
                    help="What build of the code we should use.",
                    type=str,
                    default="osx",
                    )
    ap.add_argument("--action", "-a",
                    help="What action should be performed on the tests",
                    type=str,
                    default="run"
                    )
    #ap.add_argument("--runtime", "-r",
    #                help="Keep the total running time near the time specified.",
    #                type=str,
    #                default="10m",
    #                )
    #ap.add_argument("--tags", "-t",
    #                help="Only run tests that match (all) the following tags.",
    #                action="append",
    #                default=[],
    #                )
    #ap.add_argument("--update-runtimes", "-u",
    #                help="Update the running times for all the tests we run",
    #                type=str,
    #                default=""
    #                )
                    
    ap.add_argument("tests", nargs=argparse.REMAINDER,
                    help="list of all the manually specified tests we want to run",
                    )

    options = ap.parse_args()

    #create the profile
    profile = Profile(options.profile)

    manualTests = []
    for name in options.tests:
        manualTests.append(Test(name))
    
    #scan the directory for candidate tests
    if manualTests:
        runnableTests = manualTests
    else:
        runnableTests = []
        allsubdirs = [x[0] for x in os.walk(os.path.join(testroot, "tests"))]
        for subdir in allsubdirs:
            definitionFilename=os.path.join(subdir,"definition")
            if os.path.isfile(definitionFilename):
                #scan the test looking for tags.
                #if the tags match add it to the list of tests.
                runnableTests.append(Test(nameFromTestDir(subdir)))

    #sort the runnableTests by runtime.

    if options.action == "run":
        if len(runnableTests):
            print "1.."+str(len(runnableTests))
        testCounter=1
        for test in runnableTests:
            resultFile = os.path.join(testroot, "scratch", options.build, profile.name, test.name, "result")
            if os.path.isfile(resultFile):
                os.remove(resultFile)
            test.execute(profile, options.build, options.action)
            

            if not os.path.isfile(resultFile):
                print "not ok %d %s" % (testCounter, test.name)
                print "# missing result file %s" % resultFile
            elif os.stat(resultFile).st_size != 0:
                print "not ok %d %s" % (testCounter, test.name)
                for line in open(resultFile, "r"):
                    print "# "+line.rstrip()
            else:
                print "ok %d %s" % (testCounter, test.name)
            testCounter += 1
    else:
        for test in runnableTests:
            test.execute(profile, options.build, options.action)
        
if __name__=='__main__':
    main()
