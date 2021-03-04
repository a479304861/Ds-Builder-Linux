class LCS():
  
    def input(self, x, y):
        self.x = x
        self.y = y

    def Compute_LCS(self):
        xlength = len(self.x)
        ylength = len(self.y)
        self.direction_list = [None] * xlength 
        for i in range(xlength):
            self.direction_list[i] = [None] * ylength
        self.lcslength_list = [None] * (xlength + 1)
        for j in range(xlength + 1):
            self.lcslength_list[j] = [None] * (ylength + 1)
        for i in range(0, xlength + 1):
            self.lcslength_list[i][0] = 0
        for j in range(0, ylength + 1):
            self.lcslength_list[0][j] = 0
        for i in range(1, xlength + 1):
            for j in range(1, ylength + 1):
                if self.x[i - 1] == self.y[j - 1]:
                    self.lcslength_list[i][j] = self.lcslength_list[i - 1][j - 1] + 1
                    self.direction_list[i - 1][j - 1] = 0  # left up
                elif self.lcslength_list[i - 1][j] > self.lcslength_list[i][j - 1]:
                    self.lcslength_list[i][j] = self.lcslength_list[i - 1][j]
                    self.direction_list[i - 1][j - 1] = 1  # up
                else:
                    self.lcslength_list[i][j] = self.lcslength_list[i][j - 1]
                    self.direction_list[i - 1][j - 1] = -1  # left
        self.lcslength = self.lcslength_list[-1][-1]
        return self.direction_list, self.lcslength_list

    def printOneLCS(self):
        starti = len(self.x)-1
        startj = len(self.y)-1
        final_result = []
        while startj != -1 and starti != -1:
            if self.direction_list[starti][startj] == -1:
                starti = starti
                startj = startj - 1
            elif self.direction_list[starti][startj] == 1:
                starti = starti - 1
                startj = startj
            else:
                # 
                final_result.append([starti,startj])
                starti = starti - 1
                startj = startj - 1
        return final_result[::-1]
