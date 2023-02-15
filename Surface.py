#File detailing the Surface class
class Surface:

    def __init__(self,subdivisions):

        self.N = subdivisions

        self.verticies = [
            [(i,j) for i in range(self.N)]
            for j in range(self.N)
        ]
        
        self.edges = [[

        ]]

