import pkg_resources
version = pkg_resources.get_distribution('dcore').version
#version = '3.0.0b1'

def show_version():
  print("\nYou are using the DCore version %s\n"%version)

def show_git_hash():
  print("\nYou are using the DCore git hash %s based on triqs git hash %s\n"%(dcore_hash, triqs_hash))
