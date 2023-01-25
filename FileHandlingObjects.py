class OpenFile:

    def __init__(self, open_file, count=0):
        self.f = open_file
        self.count = count
        self.current_line = None

    def read_to_line_containing(self, target_str):
        found = False
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
