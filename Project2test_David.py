import numpy as np
import enum
import matplotlib.pyplot as plt

class Polymer:

    def __init__(self, size, init_pos, color):
        self.size = size
        self.init_pos = init_pos
        self.body = [init_pos]
        self.step = 0   
        self.color = color          
        WorldMap[init_pos[0], init_pos[1]] = 1

    # Updates the coordinates in the given direction and returns the new coordinates
    def walk(self, direction):
        x_old = self.body[self.step][0]
        y_old = self.body[self.step][1]
        x_new = x_old
        y_new = y_old

        if(direction == 1):
            y_new += 1
        if(direction == 2):
            y_new -= 1
        if(direction == 3):
            x_new += 1
        if(direction == 4):
            x_new -= 1

        return [x_new, y_new]

    # Updates the polymer. It chooses a random direction and checks if the new spot is empty.
    def update(self):   

        max_it = 0  # Max number of iterations, to stop the loop if the polymer gets stuck
        while(max_it < 10):
            i = np.random.randint(1, 5)         # Choose a random direction
            new_coords = self.walk(i)           # Walk in that direction

            if(is_empty(new_coords)):           # Check is that new spot is empty. If yes, then add the new position to body 
                self.body.append(new_coords)    # and update worldmap. If no, choose new direction and repeat. 
                WorldMap[new_coords[0], new_coords[1]] = 1
                self.step += 1
                break
            max_it += 1
        pass

    # Plot polymer
    def plot(self):
        body = self.body
        x = [body[i][0] for i in range(len(body))]
        y = [body[i][1] for i in range(len(body))]
        plt.plot(x, y, color = self.color, marker = '.')

# Checks if the coordinate is occupied on the worldmap
def is_empty(coord):
    x = coord[0]
    y = coord[1]

    if(WorldMap[x, y] == 0):
        return True
    if(WorldMap[x, y] == 1):
        return False

worldsize = 5000
WorldMap = np.zeros(shape=(worldsize, worldsize))

polymers = []           # List to store the polymers
polymer_size = 1000     # Max length of polymer
num_polymer = 20        # Number of polymers

# Initialize polymers
for i in range(num_polymer):
    polymers.append(Polymer(polymer_size, [450 + 5 * i, 500], np.random.uniform(0, 1, size= 3)))

# Update polymers
for i in range(polymer_size):
    for polymer in polymers:
        polymer.update()

# Plot polymers
for polymer in polymers:
    polymer.plot()

plt.axis('equal')
plt.show()