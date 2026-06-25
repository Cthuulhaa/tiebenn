from abc import ABC, abstractmethod


class BaseClient(ABC):

    @abstractmethod
    def get_input_params(self) -> dict:
        pass

    @abstractmethod
    def run(self):
        pass
