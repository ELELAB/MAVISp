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


