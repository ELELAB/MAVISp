class MAVISpError(Exception):
    pass

class MAVISpCriticalError(Exception):
    pass

class MAVISpWarningError(Exception):
    pass

class MAVISpMultipleError(Exception):
    def __init__(self, message="", warning=list(), critical=list()):
        super().__init__(message)

        self.warning = warning
        self.critical = critical

class MAVISpEmptySystemsError(Exception):
    def __init__(self, empty_systems, message=""):
        super().__init__(message)

        self.empty_systems = empty_systems


