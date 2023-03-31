import collections
import time
import timeit

import utils


class ProgressMonitor:
    def __init__(self, total_iterations, average_n_iterations=1):
        self.total_iterations = total_iterations
        self.iteration_count = 0
        self.start_time = None
        self.iteration_durations = collections.deque([], average_n_iterations)
        self.eta = None

    def start_iteration(self):
        self.iteration_count += 1
        self.start_time = timeit.default_timer()

    def end_iteration(self):
        if not self.start_time: return None
        self.iteration_durations.append(timeit.default_timer() - self.start_time)
        average_iteration_duration = utils.mean(self.iteration_durations)
        elapsed_time = float(self.iteration_count) * average_iteration_duration
        total_time = float(self.total_iterations) * average_iteration_duration
        self.eta = time.localtime(timeit.default_timer() + (total_time - elapsed_time))

    @property
    def arrival_time_str(self):
        return time.strftime("%d/%m/%Y %H:%M", self.eta)

    @property
    def iterations_str(self):
        return str(self.iteration_count) + "/" + str(self.total_iterations)

    @property
    def progress_str(self):
        if self.eta:
            return self.iterations_str + ", ETA=" + self.arrival_time_str
        else:
            return self.iterations_str
