from KratosBRepApplication import *
application = KratosBRepApplication()
application_name = "KratosBRepApplication"
application_folder = "brep_application"

# The following lines are common for all applications
from . import application_importer
import inspect
caller = inspect.stack()[1] # Information about the file that imported this, to check for unexpected imports
application_importer.ImportApplication(application,application_name,application_folder,caller)
