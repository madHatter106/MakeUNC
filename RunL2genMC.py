#! /usr/bin/env python3

from subprocess import Popen, DEVNULL
import glob
import re
import os
import logging
import pickle
import sys
import multiprocessing as mp
import argparse
from itertools import islice
from datetime import datetime as dt

__author__ = "Erdem K."
__version__ = "0.5"


class L2genRunner():

    def __init__(self, workers, debug=False):
        maxProcs = (mp.cpu_count() - 1) * 2
        if workers > maxProcs:
            self.workers = maxProcs
        else:
            self.workers = workers
        self.dbg = debug
        self.__set_logger()

    def __set_logger(self):
        self.logger = logging.getLogger(__name__)
        formatter = logging.Formatter(' % (asctime)s [ % (levelname)s]:\
                                      % (message)s')
        if self.dbg:
            self.logger.setLevel(logging.DEBUG)
            formatter = logging.Formatter(' % (asctime)s [ % (levelname)s]\
                                          % (name)s, % (funcName)s,\
                                          % (lineno)d: % (message)s')
            ch = logging.StreamHandler()  # debug statements will appear in the console
            ch.setLevel(logging.DEBUG)
            ch.setFormatter(formatter)
            self.logger.addHandler(ch)
        else:
            self.logger.setLevel(logging.INFO)
        fh = logging.FileHandler('l2genRunner.log')
        fh.setLevel(logging.INFO)
        fh.setFormatter(formatter)
        self.logger.addHandler(fh)

    def Runner(self, cmdList):
        '''
        Creates a generator for processes then slices through the iterator
        by the number ofconcurrent processes allowed.
        cmdList is a generator yielding l2gen command lines for each process.
        '''
        status = False
        # create process generator
        processes = (Popen(cmd, shell=True, stdout=DEVNULL)
                     for cmd in cmdList)
        runningProcs = list(islice(processes, self.workers))  # start new ps
        while runningProcs:
            for i, process in enumerate(runningProcs):
                if process.poll() is not None:  # process has finished
                    runningProcs[i] = next(processes, None)  # start new ps
                    if runningProcs[i] is None:  # no new processes
                        del runningProcs[i]
                        status = True
                        break
        return status

    def GetCmdList(self):
        raise NotImplementedError


class CMCRunner(L2genRunner):
    '''
    Class to run l2gen monte carlo process, by default in parallel.
    Creates silent/noisy files  in the appropriate directories for later use
    by the uncertainty computation script.
    '''
    def __init__(self, pArgs):
        '''
        Takes pre-parsed command line arguments.
        '''

        self.l1path = pArgs.ifile
        self.l2MainPath = pArgs.opath
        self.silParFi = pArgs.prsil
        self.noiParFi = pArgs.prnoi
        self.itNum = pArgs.mcrns
        self.debug = pArgs.debug
        self.filesProcessed = 0
        self.l2SilFname = None
        self.l2NoiPath = None
        self.basename = None
        self.logfname = None
        self.logMeta = None
        super(CMCRunner, self).__init__(pArgs.workers, pArgs.debug)
        self._GetL2FilePath()
        self.logger.info("L1 file: %s" % self.l1path)
        self.logger.info("L2 main path %s" % self.l2MainPath)
        self.logger.info("silent ParFile %s" % self.silParFi)
        self.logger.info("noisy ParFile %s" % self.noiParFi)
        self.logger.info("number of iterations %d" % self.itNum)
        self.logger.info("number of concurrent processes %d" % self.workers)
        self.logger.info("silent L2 file: %s" % self.l2SilFname)
        self.logger.info("noisy L2 path: %s" % self.l2NoiPath)

    def _GetL2FilePath(self):
        '''
        Path handling and where necessary directory creation.
        '''
        pattern = '(S[0-9]+).L1A'
        basename = re.findall(pattern, self.l1path)[0]
        l2path = os.path.join(self.l2MainPath, basename)
        if not os.path.exists(l2path):
            os.makedirs(l2path)
        self.l2SilFname = os.path.join(l2path, basename+'_silent.L2')
        self.l2NoiPath = os.path.join(l2path, 'Noisy/')
        if not os.path.exists(self.l2NoiPath):
            os.makedirs(self.l2NoiPath)
        self.basename = basename

    def GetCmdList(self):
        '''Generates cmdList for subprocess calls'''
        cmdBase = 'l2gen ifile=%s ofile=' % self.l1path
        if os.path.exists(self.l2SilFname):
            if self.verbose:
                with open(self.logfname, 'a') as lf:
                    print('skipping silent L2', file=lf)
        else:
            # silent L2 does not exist, add it to the tasklist
            cmd = cmdBase + '%s par=%s' % (self.l2SilFname, self.silParFi)
            yield cmd

        for it in range(self.itNum):
            l2f = '%s_noisy_%d.L2' % (self.basename, it+1)
            ofile = os.path.join(self.l2NoiPath, l2f)
            if os.path.exists(ofile):
                if self.verbose:
                    with open(self.logfname, 'a') as lf:
                        print('skipping noisy file %s' % l2f, file=lf)
                continue
            cmd = cmdBase + '%s par=%s' % (ofile, self.noiParFi)
            yield cmd

    def Runner(self, cmdList):
        '''
        Creates a generator for processes then slices through the iterator
        by the number ofconcurrent processes allowed.
        cmdList is a generator yielding l2gen command lines for each process.
        '''
        status = False
        # create process generator
        processes = (Popen(cmd, shell=True, stdout=DEVNULL)
                     for cmd in cmdList)
        runningProcs = list(islice(processes, self.workers))  # start new ps
        while runningProcs:
            for i, process in enumerate(runningProcs):
                if process.poll() is not None:  # process has finished
                    runningProcs[i] = next(processes, None)  # start new ps
                    if runningProcs[i] is None:  # no new processes
                        del runningProcs[i]
                        status = True
                        break
        return status


class CNamespace():
    '''
    Class to replace command line argument parser for IPython calls.
    Usage: args=Namespace(ifile='',opath='',prsil='',prnoi='')
    '''

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)
        return None


class CBatchManager():
    '''
    Class to manage batch processing of multiple L1A files by the MCRunner.
    '''

    def __init__(self, bArgs, isdir=True):
        '''Takes a directory containing L1A or a text file listing
        L1Apaths on each line.'''
        if isdir:
            matchPattern = os.path.join(bArgs.ifile, '*.L1A*')
            self.ifileGen = glob.iglob(matchPattern)  # L1AGenerator
            self.pArgs = bArgs
            self.verbose = self.pArgs.verbose
            self.l2MainPath = self.pArgs.opath

            if self.verbose:
                self.logMeta = os.path.join(self.l2MainPath, 'Meta.log')

    def ProcessL1A(self):
        '''Calls L1AGenerator to get next file to process'''
        for ifile in self.ifileGen:
            self.pArgs.ifile = ifile
            mcr = CMCRunner(self.pArgs)
            pickle.dump(mcr, open(os.path.join(mcr.l2MainPath, 'mcr_%s.pkl'
                                               % mcr.basename), 'wb'))
            cmdGen = mcr.GetCmdList()
            status = mcr.Runner(cmdGen)
            if status:
                if self.verbose:
                    print('\r%s: Finished processing %s' % (dt.now(), ifile),
                          end='', flush=True)
                    with open(self.logMeta, 'a') as fmeta:
                        print('Finished processing %s' % ifile, file=fmeta)
                del mcr  # make room for the next mc set
        return None

    def CreateCmdLineArgs(**kwargs):
        pArgs = CNamespace(**kwargs)
        return pArgs


# TODO def ConsolidateParfile()

def ParseCommandLine(args):
    '''
    Returns argparse object with parsed arguments as attributes
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--ifile', help='l1a input file',
                        type=str, required='True')
    parser.add_argument('-o', '--opath', help='l2 output main path',
                        type=str, required='True')
    parser.add_argument('-s', '--prsil', help='silent param. file',
                        type=str, required='True')
    parser.add_argument('-n', '--prnoi', help='noisy param. file',
                        type=str, required='True')
    parser.add_argument('-m', '--mcrns', help='number of MC iterations',
                        type=int, default=1000)
    parser.add_argument('-w', '--workers', help='process # to allocate',
                        type=int, default=1)
    parser.add_argument('-v', '--verbose', help='increase output verbosity',
                        action='store_true')
    parser.add_argument('-b', '--batch', help='batch processing',
                        action='store_true')
    parsedArgs = parser.parse_args(args)
    # TODO parsedArgs = ConsolidateParfile(parsedArgs)
    return parsedArgs


def Main(args):

    pArgs = ParseCommandLine(args)

    if pArgs.batch:
        bcr = CBatchManager(pArgs)
        bcr.ProcessL1A()
    else:
        # Init MCRUnner Object, passing the args
        mcr = CMCRunner(pArgs)
        # Run MC process; includes creating silent file and noisy files
        taskList = mcr.GetCmdList()
        mcr.Runner(taskList)


if __name__ == '__main__':
    Main(sys.argv[1:])
