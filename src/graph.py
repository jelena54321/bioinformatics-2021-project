from abc import ABC, abstractmethod

class Graph(ABC):

    @abstractmethod
    def __init__(self, args):
        pass

    @abstractmethod
    def construct(self):
        pass

    @abstractmethod
    def generate_paths(self):
        pass

    @abstractmethod
    def generate_sequences(self):
        pass
