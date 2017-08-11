"""
A set of general functions for use in the GATK package.
"""

###############################################################################
### Function definitions
###############################################################################

def cleanup_files(files):
    """
    A method for cleaning up intermediate files.

    USAGE:
        file_cleanup(files)

    INPUT:
        * files: an array of files to cleanup

    OUTPUT:
        Returns a list containing a command for file cleanup.
    """
    program = 'rm -f'
    files_string = ' '.join(files)
    cmd = ' '.join([program, files_string])
    return {"command":cmd}
