import math
from matplotlib import pyplot
import numpy as np

x = np.linspace(0, 2 * 3.14, 100)
y = np.sin(x)

pyplot.plt(x, y)
pyplot.show()