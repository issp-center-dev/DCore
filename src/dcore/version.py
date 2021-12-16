from dcore._version import get_versions

version = get_versions()['version']

def show_version():
  print("\nYou are using the DCore version %s\n"%version)

def show_git_hash():
  print("\nYou are using the DCore git hash %s based on triqs git hash %s\n"%(dcore_hash, triqs_hash))

def print_header():
  print(f"#### DCore {version} ####")