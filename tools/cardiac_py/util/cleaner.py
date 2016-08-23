'''
@author: butler
'''
import os
import file_manipulation
import glob
class Cleaner():
    """ Class to crawl directories and clean cardioid simulations
    
    Crawl across a directory structure to check for cardioid sims and do things
    to these directories. 'os.walk' is used for simplicity. However, once a 
    cardioid simulation directory is hit things are managed manually.
    
    This means that it is unnessarily searching through the snapshot.0*
    directories. We might be able to avoid this we
    
    Assumption: cardioid directories are not recursive
    
    Parameters:
    root_directory: root directory to crawl

    """
    def __init__(self, root_directory):
        self.root = root_directory
        self.real = False
        ''' Checkpoint clean mode
        mode == 0:
        Print what mode 1 and mode 2 will do
        mode == 1:
        Delete all state files except if they occur in the final directory.
        Will ignore sshot directories which are sym-links
        mode == 2:
        Delete all state files contained in the cardioid case directory.
        Will ignore sshot directories which are sym-links
        '''
        self.mode = 1
        '''
        snapshot remove match: Provide an optional list of globs to delete from
        every snapshot directory. Warning .. this does not check whether files
        are sim-linked in.
        '''
        self.sshot_remove_matchs = []
        #self.sshot_remove_match = ["anatomy#*","coarsened_anatomy#*"]
        self.cardioid_remove_matchs = ['torusMap#000000','func_*','*.core',
                                       'fit.data','vis*.vtu','vis*.vtk']
        
    def _clean_cardioid_dir(self, directory):
        if self.real:
            file_manipulation.remove_checkpoints(directory,self.mode)
            if len(self.sshot_remove_matchs) > 0:
                for sshot_rm_match in self.sshot_remove_match:
                    globber = directory + os.sep + 'snapshot.0*' + os.sep + \
                    sshot_rm_match
                    rm_list = glob.glob(globber)
                    for rm_file in rm_list:
                        os.remove(rm_file)
            if (len(self.cardioid_remove_matchs) > 0):
                for cardioid_remove_match in self.cardioid_remove_matchs:
                    globber = directory + os.sep + cardioid_remove_match
                    rm_list = glob.glob(globber)
                    for rm_file in rm_list:
                        os.remove(rm_file)
        else:
            print 'cardioid_dir', directory
            # Remove state files

            
            
    def _is_cardioid_dir_run(self, contents):
        ''' Check function for whether a directory is a cardioid simulation dir.
        '''
        pass
        
        if ('object.data' in contents):
            # Okay probably a cardioid directory check if run 
            for path in contents:
                if 'snapshot.0' in path:
                    # WE have a valid directory return true
                    return True
            # Directory has not run
            return False
        # not a cardioid directory
        return False
    
    def run(self):
        self.sweep(self.root)
        
    def sweep(self, directory):
        contents = os.listdir(directory)
        if self._is_cardioid_dir_run(contents):
            self._clean_cardioid_dir(directory)
        else:
            for fs_obj in contents:
                fs_path = directory + os.sep + fs_obj
                if os.path.isdir(fs_path):
                    self.sweep(fs_path)
                
            
        
        
        
        
                    
        
            