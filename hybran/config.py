import tempfile

def init():
    global hybran_tmp_dir
    # Make the temporary files directory
    hybran_tmp_dir = tempfile.mkdtemp()

