import sys,os,time
import subprocess
from .tools import loadCmsProcess
from .DTWorkflow import DTWorkflow
from .DTTtrigWorkflow import DTttrigWorkflow
from .DTVdriftWorkflow import DTvdriftWorkflow
from .DTT0WireWorkflow import DTT0WireWorkflow
import logging
# setup logging
log = logging.getLogger(__name__)


class DTCalibrationWorker(object):
    """ This class serves as a top level helper to perform all available
        workflows. Additional workflow classes should use the naming scheme
        DT${WORKFLOWNAME}Workflow and implement a classmethod function add_parser_options.
    """
    available_workflows = ["ttrig","vdrift","T0Wire"]
    def __init__(self, options):
        self.options = options
        if not self.has_crab3_env:
            self.setup_crab_env()

    def run(self):
        # get class object dependent on workflow
        class_name = "DT" + self.options.workflow + "Workflow"
        workflow_class = eval(class_name)
        workflow_class_instance = workflow_class(self.options)
        workflow_class_instance.run()
        return workflow_class_instance.local_path

    @property
    def has_crab3_env(self):
        if not "/crabclient/3" in os.environ["PATH"]:
            return False
        return True

    def setup_crab_env(self):
        log.info('Setting up crab environment')
        # following
        # http://stackoverflow.com/questions/3503719/emulating-bash-source-in-python
        command = ['bash', '-c', 'unset module;source /cvmfs/cms.cern.ch/crab3/crab.sh && env']
        proc = subprocess.Popen(command, executable = '/bin/bash', stdout = subprocess.PIPE)

        for line in proc.stdout:
          (key, _, value) = line.partition(b"=")
          os.environ[key.decode()] = value.decode().replace("\n","")
        for path in os.environ['PYTHONPATH'].split(':'):
            sys.path.append(path)
        proc.communicate()

    @classmethod
    def add_arguments(cls, parser):
        workflow_parser = DTWorkflow.add_parser_options(parser)
        for workflow in cls.available_workflows:
            class_name = "DT" + workflow + "Workflow"
            try:
                workflow_class = eval( class_name )
                workflow_class.add_parser_options(workflow_parser)
            except:
                log.error("No class with name: %s exists but workflow exists in %s" %
                            (class_name, DTCalibrationWorker)
                         )

