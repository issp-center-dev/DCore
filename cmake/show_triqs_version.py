try:
    from triqs import version
except ImportError:
    print('NOT_FOUND', end='')
else:
    print(version.version, end='')