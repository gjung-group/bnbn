from .base import *
from .space import *
from .hamilton import *
from .solver import *
from .post import *

try:    # Cupy check
    import cupy
except:
    print("QUEST_Warning: Cupy module can't be imported, You cannot use GPU solver")
    pass