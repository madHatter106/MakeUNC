#!/usr/bin/env python3
from RunL2genMC import CMCRunner, CBatchManager
# from MakeUnc import MakeSwfUnc, MakeHMA
import argparse
import logging
import sys


def SetLogger(logger_name, dbg_lvl=False):
    '''

    '''
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
    parser.add_argument('--logName', help='logger name',
                        type=str, default='CallMCLog')
    parser.add_argument('-b', '--batch', help='batch processing',
                        action='store_true')
    parsedArgs = parser.parse_args(args)
    # TODO parsedArgs = ConsolidateParfile(parsedArgs)
    return parsedArgs


def Main(args):
    pArgs = ParseCommandLine(args)
    mainLogger = SetLogger(logger_name=pArgs.logName, dbg_lvl=pArgs.debug)
    mainLogger.info('Get this?')
    if pArgs.batch:
        mainLogger.info('Initializing batch processor')
        bcr = CBatchManager(pArgs)
        bcr.ProcessL1A()
    else:
        mainLogger.info('Init MCRUnner Object w/ pArgs')
        mcr = CMCRunner(pArgs, mainLogger.name)
        mainLogger.info('Creating task list')
        taskList = mcr.GetCmdList()
        mainLogger.info('Feeding tasklist l2gen runner')
        mcr.Runner(taskList)

if __name__ == '__main__':
    Main(sys.argv[1:])
