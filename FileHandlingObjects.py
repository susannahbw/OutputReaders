class OpenFile:

    def __init__(self, open_file, count=0):
        self.f = open_file
        self.count = count
        self.current_line = None

    def read_to_line_containing(self, target_str, include_current_line=False):
        found = False

        if include_current_line:
            if self.current_line is not None:
                if target_str in self.current_line:
                    return self.current_line

        while not found:
            self.count += 1
            self.current_line = self.f.readline()

            if target_str in self.current_line:
                return self.current_line
            if not self.current_line:
                break

        return 'not found'

    def read_line(self):
        self.count += 1
        self.current_line = self.f.readline()
        return self.current_line


def find_next_nondigit_character(string):
    ix = [i for i, j in enumerate(string) if not j.isdigit()]
    return ix[0]
