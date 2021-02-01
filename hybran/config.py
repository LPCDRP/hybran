import tempfile

# Setting up the temporary files directory
def init():
    global hybran_tmp_dir
    hybran_tmp_dir = tempfile.mkdtemp()

