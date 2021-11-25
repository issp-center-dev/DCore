class Intracomm:
    def __init__(self):
        pass

    def Get_rank(self):
        return 0

    def Get_size(self):
        return 1
    
    def bcast(self, x, root):
        assert root == 0
        return x

    def allreduce(self, x):
        return x

    def send(self, x, dest):
        pass

    def recv(self, source):
        pass

    def barrier(self):
        pass

world = Intracomm()

rank = world.Get_rank()
size = world.Get_size()