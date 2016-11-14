
from RunL2genMC import CMCRunner
import logging


def SetLogger(dbg=False, logger_name=__name__, logfn='CallMC.log'):
    # create logger with 'spam_application'
    logger = logging.getLogger(logger_name)
    if dbg:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    # create file handler which logs even debug messages
    fh = logging.FileHandler(logfn)
    fh.setLevel(logging.DEBUG)
    # create console handler with a higher log level
    ch = logging.StreamHandler()
    ch.setLevel(logging.ERROR)
    # create formatter and add it to the handlers
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    # add the handlers to the logger
    logger.addHandler(fh)
    logger.addHandler(ch)
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
                        action='store_true')
    parser.add_argument('-b', '--batch', help='batch processing',
                        action='store_true')
    parsedArgs = parser.parse_args(args)
    # TODO parsedArgs = ConsolidateParfile(parsedArgs)
    return parsedArgs


def Main(args):
    pArgs = ParseCommandLine(args)
    parent_logger = SetLogger(pArgs.debug)
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
