class OpenFile:

    def __init__(self, open_file, count=0):
        self.f = open_file
        self.count = count

    def read_to_line_containing(self, target_str):
        found = False
        while not found:
            self.count += 1
            line = self.f.readline()

            if target_str in line:
                return line
            if not line:
                break

        return 'not found'

    def read_line(self):
        self.count += 1
        return self.f.readline()
