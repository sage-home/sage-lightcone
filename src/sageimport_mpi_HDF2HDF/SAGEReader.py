"""
Created on 28/09/2012

@author: Amr Hassan
"""
import math
import numpy
import logging
import h5py
import time
import ProcessTreeTraversal
from mpi4py import MPI  # MPI Implementation


def print_attrs(name, obj):
    if name == 'galaxies' or name == 'Data':
        global mydtype
        mydtype = obj.dtype
        print(len(obj), len(obj.dtype))
        print(obj.dtype)

class SAGEDataReader:    
    # The Module handles the data reading from SAGE output to a memory data structure.
    
    CurrentInputFilePath = ""
    CurrentOutputFilePath = ""
    CurrentGlobalTreeID = 0
    FormatMapping = {'int': 'i', 'float': 'f', 'long long': 'q'}
    
    
    
    def __init__(self, CurrentSAGEStruct, Options, sageFieldsNode, CommSize, CommRank):
        
        
        # Initialize the Class to handle a specific file path        
        self.CurrentInputFilePath = Options['RunningSettings:InputFile']
        self.CurrentOutputFilePath = Options['RunningSettings:OutputFile']
        self.CurrentSAGEStruct = CurrentSAGEStruct
        self.Options = Options
        self.sageFieldsNode = sageFieldsNode
        self.CommSize = CommSize
        self.CommRank = CommRank

        # logging
        self.localgalaxyid_status_needs_printing = True

    def copyGroup(self, name):
        group_path = self.InputFile[name].parent.name
        group_id = self.OutputFile.require_group(group_path)
        self.InputFile.copy(name, group_id, name=name)

    def createGalaxies(self):
        # create new composite datatype for output galaxies - adding some fields
        self.InputFile.visititems(print_attrs)
        dynamic_dtype_description = mydtype.descr
        dynamic_dtype_description.append(('globaltreeid', '<i8'))
        dynamic_dtype_description.append(('breadthfirst_traversalorder', '<i8'))
        dynamic_dtype_description.append(('depthfirst_traversalorder', '<i8'))
        dynamic_dtype_description.append(('subtree_count', '<i8'))
        dynamic_dtype_description.append(('localgalaxyid', '<i4'))
        self.odtype = numpy.dtype(dynamic_dtype_description)
        self.OutputFile.create_dataset('galaxies', self.InputFile['galaxies'].shape, dtype=self.odtype)

        # Create the sidecar xml
        filename_sidecar = self.CurrentOutputFilePath[0:-3] + ".xml"
        # filename_sidecar = "abc.txt"
        with open(filename_sidecar, mode="w+") as xml_sidecar:
            try:
                xml_sidecar.write("  <sageinput>\n")
                count = 1
                for field in self.sageFieldsNode:
                    xml_sidecar.write(f"    <Field Type=\"{field.attrib['Type']}\"\n")
                    xml_sidecar.write(f"      label=\"{field.attrib['label']}\"\n")
                    xml_sidecar.write(f"      description=\"{field.attrib['description']}\"\n")
                    xml_sidecar.write(f"      order=\"{field.attrib['order']}\"\n")
                    xml_sidecar.write(f"      units=\"{field.attrib['units']}\"\n")
                    xml_sidecar.write(f"      group=\"{field.attrib['group']}\">{field.text}</Field>\n")
                    count = count + 1
                xml_sidecar.write(f"    <Field Type=\"{'long long'}\"\n")
                xml_sidecar.write(f"      label=\"{'globaltreeid'}\"\n")
                xml_sidecar.write(f"      description=\"{'globaltreeid'}\"\n")
                xml_sidecar.write(f"      order=\"{count}\"\n")
                xml_sidecar.write(f"      units=\"{''}\"\n")
                xml_sidecar.write(f"      group=\"{'internal'}\">{'globaltreeid'}</Field>\n")
                count = count + 1
                xml_sidecar.write(f"    <Field Type=\"{'long long'}\"\n")
                xml_sidecar.write(f"      label=\"{'breadthfirst_traversalorder'}\"\n")
                xml_sidecar.write(f"      description=\"{'breadthfirst_traversalorder'}\"\n")
                xml_sidecar.write(f"      order=\"{count}\"\n")
                xml_sidecar.write(f"      units=\"{''}\"\n")
                xml_sidecar.write(f"      group=\"{'internal'}\">{'breadthfirst_traversalorder'}</Field>\n")
                count = count + 1
                xml_sidecar.write(f"    <Field Type=\"{'long long'}\"\n")
                xml_sidecar.write(f"      label=\"{'depthfirst_traversalorder'}\"\n")
                xml_sidecar.write(f"      description=\"{'depthfirst_traversalorder'}\"\n")
                xml_sidecar.write(f"      order=\"{count}\"\n")
                xml_sidecar.write(f"      units=\"{''}\"\n")
                xml_sidecar.write(f"      group=\"{'internal'}\">{'depthfirst_traversalorder'}</Field>\n")
                count = count + 1
                xml_sidecar.write(f"    <Field Type=\"{'long long'}\"\n")
                xml_sidecar.write(f"      label=\"{'subtree_count'}\"\n")
                xml_sidecar.write(f"      description=\"{'subtree_count'}\"\n")
                xml_sidecar.write(f"      order=\"{count}\"\n")
                xml_sidecar.write(f"      units=\"{''}\"\n")
                xml_sidecar.write(f"      group=\"{'internal'}\">{'subtree_count'}</Field>\n")
                count = count + 1
                xml_sidecar.write(f"    <Field Type=\"{'int'}\"\n")
                xml_sidecar.write(f"      label=\"{'localgalaxyid'}\"\n")
                xml_sidecar.write(f"      description=\"{'localgalaxyid'}\"\n")
                xml_sidecar.write(f"      order=\"{count}\"\n")
                xml_sidecar.write(f"      units=\"{''}\"\n")
                xml_sidecar.write(f"      group=\"{'internal'}\">{'localgalaxyid'}</Field>\n")
                xml_sidecar.write("  </sageinput>\n")
                xml_sidecar.close()
            except Exception as e:
                print(e)

    def ProcessAllTrees(self):

        
        print("ProcessAllTrees Rank=", self.CommRank)
        # Local change until I have hdf5 working with mpi
        # self.InputFile = h5py.File(self.CurrentInputFilePath, 'r', driver='mpio', comm=MPI.COMM_WORLD)
        self.InputFile = h5py.File(self.CurrentInputFilePath, 'r')


        print("\t Opened InputFile Rank=", self.CommRank)

        # Local change until I have hdf5 working with mpi
        # self.OutputFile = h5py.File(self.CurrentOutputFilePath, 'w', driver='mpio', comm=MPI.COMM_WORLD)
        self.OutputFile = h5py.File(self.CurrentOutputFilePath, 'w')
        print("\t Opened Output Rank=", self.CommRank)

        print("\t Copy Groups Rank=", self.CommRank)
        self.copyGroup('cosmology')
        self.createGalaxies()
        self.copyGroup('snapshot_redshifts')
        self.copyGroup('tree_counts')
        self.copyGroup('tree_displs')

        # Everyone wait for master to copy the header to output
        print("\t Waiting Rank=", self.CommRank)
        MPI.COMM_WORLD.Barrier()
        print("\t ...proceed Rank=", self.CommRank)

        # Process All the Trees

        tree_counts = self.InputFile['tree_counts']
        tree_displs = self.InputFile['tree_displs']
        
        TotalNumberofUnprocessedTrees = len(tree_counts)
        MyNumberOfUnprocessedTrees = TotalNumberofUnprocessedTrees // self.CommSize + 1
        MyFirstTree = self.CommRank * MyNumberOfUnprocessedTrees
        MyLastTree = min(MyFirstTree + MyNumberOfUnprocessedTrees, TotalNumberofUnprocessedTrees)

        offset = 0
        for i in range(0, MyFirstTree):
            offset += tree_counts[i]
        print(f"CommRank[{self.CommRank}]:First={MyFirstTree} Last={MyLastTree-1} Offset={offset} TotalTrees={TotalNumberofUnprocessedTrees}")
        percentage = 0
        for UnProcessedTree in range(MyFirstTree, MyLastTree):
            # Updating the user with what is going on
            # logging.info(f"{self.CommRank}:Processing Tree ({0}-{TotalNumberofUnprocessedTrees-1}):{UnProcessedTree}")
            offset = self.ProcessTree(UnProcessedTree, tree_counts, tree_displs, offset)
            sofar = int((UnProcessedTree - MyFirstTree + 1) * 100 / (MyLastTree - MyFirstTree + 1))
            if sofar != percentage:
                percentage = sofar
                if self.CommRank == 0:
                    print(f"Progress:{sofar}")
                    logging.info(f"Progress:{sofar}")

        # Wait till everyone is done
        print("\t Waiting Rank=", self.CommRank)
        MPI.COMM_WORLD.Barrier()
        self.InputFile.close()
        self.OutputFile.close()
        print("\t ...done  Rank=", self.CommRank)


    def GenerateDictFromFields(self, TreeLoadingID, TreeData, offset):
        TreeDict = []
        
        pgcopy_dtype = [('num_fields', '>i2')]
        FieldsList = []
        FieldsIndex = 0
        for field, dtype in TreeData.dtype.descr:
            FieldsList += [self.CurrentSAGEStruct[FieldsIndex][0]]
            FieldName = self.CurrentSAGEStruct[FieldsIndex][0]
            #print(f"ExistingField={FieldName}:\t{dtype}")

            pgcopy_dtype += [(FieldName + '_length', '>i4'), (FieldName, dtype.replace('<', '>'))]
            FieldsIndex = FieldsIndex + 1
        
        
        
        ####### Add Generated Fields (Computed) ###############################
        
        FieldName = 'globaltreeid'
        pgcopy_dtype += [(FieldName + '_length', '>i4'), (FieldName, '>i8')]
        
        # FieldName='CentralGalaxyGlobalID'        
        # pgcopy_dtype += [(FieldName + '_length', '>i4'),(FieldName, '>i8')]
        
        
        
        FieldName = 'breadthfirst_traversalorder'
        pgcopy_dtype += [(FieldName + '_length', '>i4'), (FieldName, '>i8')]
        
        FieldName = 'depthfirst_traversalorder'
        pgcopy_dtype += [(FieldName + '_length', '>i4'), (FieldName, '>i8')]
        
        FieldName = 'subtree_count'
        pgcopy_dtype += [(FieldName + '_length', '>i4'), (FieldName, '>i8')]
        
        
        FieldsList += ['globaltreeid']
        # FieldsList+=['centralgalaxyglobalid']       
        FieldsList += ['breadthfirst_traversalorder']
        FieldsList += ['depthfirst_traversalorder']
        FieldsList += ['subtree_count']
        #########################################################################
        
        # TODO: only print this once
        if self.localgalaxyid_status_needs_printing:
            self.localgalaxyid_status_needs_printing = False
            if FieldsList.count('localgalaxyid') > 0:
                logging.info("### localgalaxyid already Exists. No Data Will be generated")
            else:
                logging.info("### localgalaxyid  is Missing. Regenerate Local GalaxyID")
        
        if FieldsList.count('localgalaxyid') == 0:
            FieldName = 'localgalaxyid'
            pgcopy_dtype += [(FieldName + '_length', '>i4'), (FieldName, '>i4')]


        pgcopy = numpy.empty(TreeData.shape, pgcopy_dtype)
        
        
        
        pgcopy['globaltreeid_length'] = numpy.dtype('>i8').alignment
        # pgcopy['CentralGalaxyGlobalID_length'] = numpy.dtype('>i8').alignment
                       
        pgcopy['breadthfirst_traversalorder_length'] = numpy.dtype('>i8').alignment
        pgcopy['depthfirst_traversalorder_length'] = numpy.dtype('>i8').alignment
        pgcopy['subtree_count_length'] = numpy.dtype('>i8').alignment

        
        GeneratedFields = 0
        if FieldsList.count('localgalaxyid') == 0:
            GeneratedFields = 5
        else:
            GeneratedFields = 4
        
        
        pgcopy['num_fields'] = len(TreeData.dtype) + GeneratedFields
        #print(f" position_x={TreeData['position_x']}")
        #print(f" position_y={TreeData['position_y']}")
        #print(f" position_z={TreeData['position_z']}")

        for i in range(0, len(TreeData.dtype)):
            field = self.CurrentSAGEStruct[i][0]                            
            pgcopy[field + '_length'] = TreeData.dtype[i].alignment
            pgcopy[field] = TreeData[TreeData.dtype.names[i]]


        pgcopy['globaltreeid'].fill(TreeLoadingID)
        # Turn globaltreeid into a globaltreeindex (An index into the galaxyTree array)
        pgcopy['globaltreeid'] = range(offset, offset + len(TreeData))


        
        if FieldsList.count('localgalaxyid') == 0:
            pgcopy['localgalaxyid'] = range(0, len(TreeData))
            pgcopy['localgalaxyid_length'] = numpy.dtype('>i4').alignment




        
        ManageTreeIndexObj = ProcessTreeTraversal.ManageTreeIndex(self.Options)
        
        
        ManageTreeIndexObj.BuildTree(TreeData)
        
        ManageTreeIndexObj.BreadthFirst(ManageTreeIndexObj.ParentNode)
        ManageTreeIndexObj.DepthFirst_PreOrder(ManageTreeIndexObj.ParentNode)
        ManageTreeIndexObj.CountChildNodes(ManageTreeIndexObj.ParentNode)
        
        NodesList = {}
        ManageTreeIndexObj.TreeToList(ManageTreeIndexObj.ParentNode, NodesList)      
        
        logging.info(f'NodesList Length={len(NodesList)}')
        logging.info(f'Tree Loading ID={TreeLoadingID}')

        if not (len(NodesList) == len(TreeData)):
            logging.info('NodesList.Length!= TreeData.Length')
            logging.info(f"NodeList Length={len(NodesList)}")
            logging.info(f"TreeData Length={len(TreeData)}")
            for i in range(0, len(TreeData)):
                if TreeData[i][self.Options['TreeMapping_0']] not in NodesList:
                    logging.info(f"{TreeData[i][self.Options['TreeMapping_0']]} is missing!")
        

        for TreeField in pgcopy:         

            GlobalIndex = TreeField[self.Options['TreeMapping_0']]  
            
            TreeField['breadthfirst_traversalorder'] = NodesList[GlobalIndex]['BreadthFirstIndex']
            TreeField['depthfirst_traversalorder'] = NodesList[GlobalIndex]['DepthFirstIndex']
            snapnum = NodesList[GlobalIndex]['SnapNum']
            subsize = NodesList[GlobalIndex]['SubTreeSize']
            #print(f'all={snapnum}\tGlobalIndex={GlobalIndex}\t SubTreeSize={subsize}')  
            TreeField['subtree_count'] = NodesList[GlobalIndex]['SubTreeSize']
            #print(f"  +LocalIndex={i}\t SubtreeCount={TreeField['subtree_count']}  x={TreeField['posx']} xlength={TreeField['posx_length']}")
        ###### This Part Work for Validation only ############
        ###### Please un-comment it only when you have a data issue #######
        ###### It might significantly reduce the importing speed ########
        # SampleTreeItem=pgcopy[0]
        # ActualFieldsCount=len(SampleTreeItem)
        # logging.info(f"Actual Number of Fields:{ActualFieldsCount}")
        # logging.info(f"Number of Data Fields:{SampleTreeItem[0]}\t Number of actual fields should be= {(SampleTreeItem[0]*2)+1}")
        # for i in range(1,ActualFieldsCount,2):
        #     logging.info(f"{SampleTreeItem[i]}\t{type(SampleTreeItem[i+1])}\t{numpy.dtype(SampleTreeItem[i+1]).itemsize}")
        #     if(numpy.dtype(SampleTreeItem[i+1]).itemsize != SampleTreeItem[i]):
        #         logging.error(f"Data Size Does not match at :{i}")

        return pgcopy
            
    def ProcessTree(self, UnProcessedTree, tree_counts, tree_displs, offset):
        
        
        LoadingTreeID = UnProcessedTree
        StartIndex = tree_displs[UnProcessedTree]
        GalaxiesCount = tree_counts[UnProcessedTree]

        logging.info('\t '+str(self.CommRank)+': Number of Galaxies in Tree ('+str(LoadingTreeID)+')='+str(GalaxiesCount))
        if GalaxiesCount > 0:
            start_time = time.time()
            TreeData=self.InputFile['galaxies'][StartIndex:StartIndex+GalaxiesCount] 
            logging.info("Reading Data="+str( time.time() - start_time)+ " seconds")
            #print(f" READING TreeID={LoadingTreeID}\t StartIndex={StartIndex}\t GalaxiesCount={GalaxiesCount} Data={TreeData[0]}")
            start_time = time.time()      
            TreeData=self.GenerateDictFromFields(LoadingTreeID, TreeData, offset)
            self.OutputFile['galaxies'][offset:offset+len(TreeData)] = TreeData
            #print(f" WROTE TreeID={LoadingTreeID}\t StartIndex={StartIndex}\t GalaxiesCount={GalaxiesCount} Data={TreeData}")
            for i in range(0, len(TreeData)):
                if TreeData[i]['subtree_count'] == 0:
                    print(f"Data Issue in TreeID={LoadingTreeID}\t StartIndex={StartIndex}\t GalaxiesCount={GalaxiesCount} at LocalIndex={i} SubtreeCount={TreeData[i]['subtree_count']}")
                #print(f"  LocalIndex={i}\t SubtreeCount={TreeData[i]['subtree_count']}  x={TreeData[i]['posx']} xlength={TreeData[i]['posx_length']}")
            offset = offset + len(TreeData)
            start_time = time.time()
            #self.ComputeFields(TreeData)
            #logging.info("Compute Fields="+str( time.time() - start_time)+ " seconds")
            start_time = time.time()
         
        return offset











