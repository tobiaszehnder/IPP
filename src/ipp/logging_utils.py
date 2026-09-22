"""Logging that coexists with a tqdm progress bar (pauses/refreshes around prints)."""

LOG_LEVEL_QUIET = 0
LOG_LEVEL_LOG = 1
LOG_LEVEL_DEBUG = 2


class ProjectionLogger:
    """Level-gated print helper; `.pbar` set to a tqdm bar pauses around prints."""

    def __init__(self, level: int = LOG_LEVEL_LOG) -> None:
        self.level = level
        self.pbar = None

    def is_log(self) -> bool:
        return self.level >= LOG_LEVEL_LOG

    def is_debug(self) -> bool:
        return self.level >= LOG_LEVEL_DEBUG

    def do_print(self, *args, **kwargs) -> None:
        """Print unconditionally, pausing/refreshing the progress bar if one is set."""
        if self.pbar:
            self.pbar.clear()
        print(*args, **kwargs)
        if self.pbar:
            self.pbar.refresh()

    def log(self, *args, **kwargs) -> None:
        if self.is_log():
            self.do_print(*args, **kwargs)

    def debug(self, *args, **kwargs) -> None:
        if self.is_debug():
            self.do_print(*args, **kwargs)
