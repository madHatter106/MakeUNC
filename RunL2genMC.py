#! /usr/bin/env python3

from subprocess import Popen, DEVNULL, PIPE
import glob
import re
import os
import logging
import pickle
import sys
import multiprocessing as mp
import argparse
from itertools import islice
from datetime import datetime as DT

__author__ = "Erdem K."
__version__ = "0.5"


class L2genRunner():

    def __init__(self, workers):
        maxProcs = (mp.cpu_count() - 1) * 2
        if workers > maxProcs:
            self.workers = maxProcs
        else:
            self.workers = workers

    def Runner(self, cmdList):
        '''
        Creates a generator for processes then slices through the iterator
        by the number ofconcurrent processes allowed.
        cmdList is a generator yielding l2gen command lines for each process.
        '''
        status = False
        # create process generator
        processes = (Popen(cmd, shell=True, stderr=DEVNULL, stdout=DEVNULL)
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
    def __init__(self, pArgs, parent_logger_name):
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
        super(CMCRunner, self).__init__(pArgs.workers)
        self._GetL2FilePath()
        self._SetLogger(parent_logger_name)
        self.logger.info("L1 file: %s" % self.l1path)
        self.logger.info("L2 main path %s" % self.l2MainPath)
        self.logger.info("silent ParFile %s" % self.silParFi)
        self.logger.info("noisy ParFile %s" % self.noiParFi)
        self.logger.info("number of iterations %d" % self.itNum)
        self.logger.info("number of concurrent processes %d" % self.workers)
        self.logger.info("silent L2 file: %s" % self.l2SilFname)
        self.logger.info("noisy L2 path: %s" % self.l2NoiPath)

    def _SetLogger(self, pln):
        '''
        The user is expected to set the logger in the 'main' module and
        pass on its name to this child logger. Otherwise no logging will occurr.
        '''
        self.logger = logging.getLogger('%s.RunL2genMC.CMCRunner' % pln)
        self.logger.info('%s initialized' % self.logger.name)

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

    def BuildCmdGen(self):
        '''Generator: generates cmdList for subprocess calls'''
        cmdBase = 'l2gen ifile=%s ofile=' % self.l1path
        if os.path.exists(self.l2SilFname):
            self.logger.info('skipping silent L2')
        else:
            # silent L2 does not exist, add it to the tasklist
            cmd = cmdBase + '%s par=%s' % (self.l2SilFname, self.silParFi)
            yield cmd

        for it in range(self.itNum):
            l2f = '%s_noisy_%d.L2' % (self.basename, it+1)
            ofile = os.path.join(self.l2NoiPath, l2f)
            if os.path.exists(ofile):
                self.logger.info('skipping noisy file %s' % l2f)
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
        processes = (Popen(cmd, shell=True, stdout=DEVNULL, stderr=PIPE)
                     for cmd in cmdList)
        runningProcs = list(islice(processes, self.workers))  # start new ps
        while runningProcs:
            for i, process in enumerate(runningProcs):
                if self.debug:
                    if process.stderr:
                        for line in process.stderr.readlines():
                            self.logger.debug('%d %s' % (i, line))
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

    def __init__(self, bArgs, isdir=True, parent_logger_name=None):
        '''Takes a directory containing L1A or a text file listing
        L1Apaths on each line.'''
        if isdir:
            matchPattern = os.path.join(bArgs.ifile, '*.L1A*')
            self.ifileGen = glob.iglob(matchPattern)  # L1AGenerator
            self.pArgs = bArgs
            self.verbose = self.pArgs.verbose
            self.l2MainPath = self.pArgs.opath
            if parent_logger_name is not None:
                self._SetLogger(parent_logger_name)

    def _SetLogger(self, parentloggername):
        '''
        The user is expected to set the logger in the 'main' module and
        pass on its name to this child logger. Otherwise no logging will occurr.
        '''
        self.logger = logging.getLogger('%s.RunL2genMC.CMCRunner' % parentloggername)
        self.logger.info('logger initialized')

    def ProcessL1A(self):
        '''Calls L1AGenerator to get next file to process'''
        for ifile in self.ifileGen:
            self.pArgs.ifile = ifile
            mcr = CMCRunner(self.pArgs)
            pickle.dump(mcr, open(os.path.join(mcr.l2MainPath, 'mcr_%s.pkl'
                                               % mcr.basename), 'wb'))
            cmdGen = mcr.BuildCmdGen()
            status = mcr.Runner(cmdGen)
            if status:
                if self.verbose:
                    print('\r%s: Finished processing %s' % (DT.now(), ifile),
                          end='', flush=True)
                    with open(self.logMeta, 'a') as fmeta:
                        print('Finished processing %s' % ifile, file=fmeta)
                del mcr  # make room for the next mc set
        return None

    def CreateCmdLineArgs(**kwargs):
        pArgs = CNamespace(**kwargs)
        return pArgs


# TODO def ConsolidateParfile()

def SetLogger(dbg_lvl=False):
    '''

    '''
    logger_name = 'RL2MC_%s' % DT.strftime(DT.now(), '%Y-%m-%dT%H:%M:%S')
    logfn = '%s.log' % logger_name
    logger = logging.getLogger(logger_name)
    if dbg_lvl:
        logger.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s -'
                                      + ' [%(module)s..%(funcName)s..%(lineno)d]'
                                      + ' - %(message)s')
    else:
        logger.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh = logging.FileHandler(logfn)
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    ch = logging.StreamHandler()
    ch.setLevel(logging.ERROR)
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    logger.addHandler(fh)
    logger.debug('logging')
    return logger


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
    parser.add_argument('-d', '--debug', help='increase output verbosity',
                        action='store_true', default=False)
    parser.add_argument('-b', '--batch', help='batch processing',
                        action='store_true')
    parsedArgs = parser.parse_args(args)
    # TODO parsedArgs = ConsolidateParfile(parsedArgs)
    return parsedArgs


def Main(args):

    pArgs = ParseCommandLine(args)
    mainLogger = SetLogger(dbg_lvl=pArgs.debug)

    if not os.path.exists(pArgs.ifile):
        sys.exit('\n %s not found!\n exiting...' % pArgs.ifile)
    if pArgs.batch:
        mainLogger.info('Initializing batch processor')
        bcr = CBatchManager(pArgs)
        bcr.ProcessL1A()
    else:
        mainLogger.info('Init MCRUnner Object w/ pArgs')
        mcr = CMCRunner(pArgs, mainLogger.name)
        mainLogger.info('Creating task list')
        taskList = mcr.BuildCmdGen()
        mainLogger.info('Feeding tasklist l2gen runner')
        mcr.Runner(taskList)


if __name__ == '__main__':
    Main(sys.argv[1:])
