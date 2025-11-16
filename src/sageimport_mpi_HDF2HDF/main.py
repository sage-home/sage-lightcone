"""
Adapted from sageimport_mpi_HDF
Main File
Execute the importing code
JUST SERIAL FOR NOW
@author: Ray Seikel
"""
## Import Helper modules
import sys  # for listing directory contents
from mpi4py import MPI  # MPI Implementation
import time
import logging
import os

import SAGEReader  # Read the SAGE Files into memory
import settingReader  # Read the XML settings

mydtype = []

def SetupLogFile(CommRank):
    FilePath = f'log/logfile{CommRank}.log'
    if os.path.exists(FilePath):
        os.remove(FilePath)
    logging.basicConfig(filename=FilePath, level=logging.DEBUG, format='%(asctime)s %(message)s')

if __name__ == '__main__':

    if len(sys.argv) < 2:
        print("Error Not Enough Arguments")
        exit()
    SettingFile = sys.argv[1]
    try:
        ResumeProcessing = (sys.argv[2] == 'True')
    except IndexError:
        ResumeProcessing = False
    
    ## MPI already initiated in the import statement
    ## Get The current Process Rank and the total number of processes
    start = time.perf_counter()  # Changed from time.clock() which is deprecated
    comm = MPI.COMM_WORLD
    CommRank = comm.Get_rank()
    
    CommSize = comm.Get_size()

    SetupLogFile(CommRank)
    
    if CommRank == 0:
        logging.info('SAGE Data Importing Tool ( MPI version)')
    
    ### Read Running Settings
    CurrentSAGEStruct, Options, sageFieldsNode = settingReader.ParseParams(SettingFile)

    comm.Barrier()
    logging.info("Start together")

    ## Init files reader
    Reader = SAGEReader.SAGEDataReader(CurrentSAGEStruct, Options, sageFieldsNode, CommSize, CommRank)
    ## Start Processing the files
    Reader.ProcessAllTrees()

    comm.Barrier()
    logging.info("Finish together")

    logging.info(f'{CommRank}:Processing Done')
    sys.exit(0)
